"""
MultiModalBPNet: BPNet variant taking DNA sequence + base-pair accessibility as input.

Input: (N, 5, L) tensor
    Channels 0-3: one-hot DNA sequence
    Channel 4:    accessibility signal (ATAC or DNase, log1p-normalized)

Middle fusion architecture:
    seq (4ch) → Conv1d → (N, n_filters, L)
                                             → cat → (N, n_filters+n_acc_filters, L)
    acc (1ch) → Conv1d → (N, n_acc_filters, L)
                                             → dilated residual layers
                                             → profile head + counts head

DeepLIFT/SHAP attributions on the single 5-channel input give:
    attr[:, :4, :] = sequence importance
    attr[:, 4:, :]  = accessibility importance
"""

import time
import numpy
import torch

from bpnetlite.losses import MNLLLoss, log1pMSELoss, _mixture_loss
from bpnetlite.performance import calculate_performance_measures
from bpnetlite.logging import Logger

from tangermeme.predict import predict


class MultiModalBPNet(torch.nn.Module):
    """BPNet model with DNA sequence and base-pair accessibility inputs.

    Parameters
    ----------
    n_filters: int
        Filters for the sequence initial convolution. Default 64.
    n_acc_filters: int
        Filters for the accessibility initial convolution. Default 8.
    n_layers: int
        Number of dilated residual layers. Default 8.
    n_outputs: int
        Number of profile output tracks (1 for unstranded). Default 1.
    count_loss_weight: float
        Weight on the counts loss term. Default 1.
    profile_output_bias: bool
        Whether to include bias in the final profile convolution. Default True.
    count_output_bias: bool
        Whether to include bias in the counts linear layer. Default True.
    name: str or None
        Prefix for saved model files. Default None (auto-generated).
    trimming: int or None
        Bases trimmed from each side of the input to produce the output window.
        Default None (uses BPNet formula: 47 + sum(2^i for i in 1..n_layers)).
    verbose: bool
        Whether to print training statistics. Default True.
    """

    def __init__(self, n_filters=64, n_acc_filters=8, n_layers=8, n_outputs=2,
                 count_loss_weight=1, profile_output_bias=True,
                 count_output_bias=True, name=None, trimming=None, verbose=True):
        super().__init__()
        self.n_filters = n_filters
        self.n_acc_filters = n_acc_filters
        self.n_layers = n_layers
        self.n_outputs = n_outputs
        self.count_loss_weight = count_loss_weight
        self.name = name or f"multimodal_bpnet.{n_filters}.{n_acc_filters}.{n_layers}"
        self.trimming = trimming or 47 + sum(2**i for i in range(1, n_layers + 1))

        n_merged = n_filters + n_acc_filters

        # Separate initial convolutions for each modality
        self.seq_conv = torch.nn.Conv1d(4, n_filters, kernel_size=21, padding=10)
        self.seq_relu = torch.nn.ReLU()
        self.acc_conv = torch.nn.Conv1d(1, n_acc_filters, kernel_size=21, padding=10)
        self.acc_relu = torch.nn.ReLU()

        # Dilated residual layers on merged representation
        self.rconvs = torch.nn.ModuleList([
            torch.nn.Conv1d(n_merged, n_merged, kernel_size=3,
                            padding=2**i, dilation=2**i)
            for i in range(1, n_layers + 1)
        ])
        self.rrelus = torch.nn.ModuleList([
            torch.nn.ReLU() for _ in range(n_layers)
        ])

        # Profile head
        self.fconv = torch.nn.Conv1d(n_merged, n_outputs, kernel_size=75,
                                     padding=37, bias=profile_output_bias)

        # Counts head
        self.linear = torch.nn.Linear(n_merged, 1, bias=count_output_bias)

        self.logger = Logger(
            ["Epoch", "Iteration", "Training Time", "Validation Time",
             "Training MNLL", "Training Count MSE",
             "Validation MNLL", "Validation Profile Pearson",
             "Validation Count Pearson", "Validation Count MSE", "Saved?"],
            verbose=verbose
        )

    def forward(self, X):
        """Forward pass.

        Parameters
        ----------
        X: torch.Tensor, shape (N, 5, L)
            Channels 0-3: one-hot sequence; channel 4: accessibility signal.

        Returns
        -------
        y_profile: torch.Tensor, shape (N, n_outputs, out_length)
        y_counts:  torch.Tensor, shape (N, 1)
        """
        start, end = self.trimming, X.shape[2] - self.trimming

        X_seq = self.seq_relu(self.seq_conv(X[:, :4, :]))
        X_acc = self.acc_relu(self.acc_conv(X[:, 4:, :]))
        X_merged = torch.cat([X_seq, X_acc], dim=1)

        for i in range(self.n_layers):
            X_conv = self.rrelus[i](self.rconvs[i](X_merged))
            X_merged = torch.add(X_merged, X_conv)

        y_profile = self.fconv(X_merged)[:, :, start:end]

        X_pooled = torch.mean(X_merged[:, :, start - 37:end + 37], dim=2)
        y_counts = self.linear(X_pooled).reshape(X.shape[0], 1)

        return y_profile, y_counts

    def fit(self, training_data, optimizer, scheduler=None,
            X_valid=None, y_valid=None, max_epochs=100, batch_size=64,
            dtype='float32', device='cuda', early_stopping=None):
        """Train the model.

        Parameters
        ----------
        training_data: DataLoader
            Yields (X, y, labels) tuples where X has shape (N, 5, in_length).
        optimizer: torch.optim.Optimizer
        scheduler: lr_scheduler or None
        X_valid: torch.Tensor or None, shape (n, 5, in_length)
        y_valid: torch.Tensor or None, shape (n, n_outputs, out_length)
        max_epochs: int
        batch_size: int
        dtype: str or torch.dtype
        device: str
        early_stopping: int or None
        """
        dtype = getattr(torch, dtype) if isinstance(dtype, str) else dtype
        device_type = device.split(':')[0]
        self.to(device)

        iteration = 0
        early_stop_count = 0
        best_loss = float("inf")
        self.logger.start()

        for epoch in range(max_epochs):
            tic = time.time()

            for data in training_data:
                X, y, labels = data[0], data[-2], data[-1]
                X = X.to(device).float()
                y = y.to(device)
                labels = labels.to(device)

                optimizer.zero_grad()
                self.train()

                with torch.autocast(device_type=device_type, dtype=dtype):
                    y_hat_logits, y_hat_logcounts = self(X)
                    train_profile_loss, train_count_loss, loss = _mixture_loss(
                        y, y_hat_logits, y_hat_logcounts,
                        self.count_loss_weight, labels
                    )

                loss.backward()
                torch.nn.utils.clip_grad_norm_(self.parameters(), 1)
                optimizer.step()
                iteration += 1

            train_time = time.time() - tic

            with torch.no_grad():
                self.eval()
                tic = time.time()

                y_hat_logits, y_hat_logcounts = predict(
                    self, X_valid, batch_size=batch_size,
                    dtype=dtype, device=device
                )

                valid_profile_loss, valid_count_loss, valid_loss = _mixture_loss(
                    y_valid, y_hat_logits, y_hat_logcounts, self.count_loss_weight
                )

                measures = calculate_performance_measures(
                    y_hat_logits, y_valid, y_hat_logcounts,
                    kernel_sigma=7, kernel_width=81,
                    measures=['profile_pearson', 'count_pearson']
                )

                valid_profile_corr = numpy.nan_to_num(measures['profile_pearson'])
                valid_count_corr = numpy.nan_to_num(measures['count_pearson'])
                valid_time = time.time() - tic

                self.logger.add([
                    epoch, iteration, train_time, valid_time,
                    train_profile_loss.item(), train_count_loss.item(),
                    valid_profile_loss.item(),
                    valid_profile_corr.mean(), valid_count_corr.mean(),
                    valid_count_loss.item(),
                    (valid_loss < best_loss).item()
                ])
                self.logger.save(f"{self.name}.log")

                if valid_loss < best_loss:
                    torch.save(self, f"{self.name}.torch")
                    best_loss = valid_loss
                    early_stop_count = -1

            if scheduler is not None:
                scheduler.step(valid_loss)

            early_stop_count += 1
            if early_stopping is not None and early_stop_count >= early_stopping:
                break

        torch.save(self, f"{self.name}.final.torch")
