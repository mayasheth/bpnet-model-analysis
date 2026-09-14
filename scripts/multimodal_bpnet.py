"""
MultiModalBPNet: BPNet variant supporting three input modes.

Modes:
  'multimodal'  Input (N, 5, L): channels 0-3 one-hot DNA, channel 4 accessibility.
                Middle fusion: seq (4→n_filters) + acc (1→n_acc_filters) → concat →
                dilated stack (n_filters+n_acc_filters channels).
  'sequence'    Input (N, 4, L): one-hot DNA only.
                seq (4→n_filters) → dilated stack (n_filters channels).
  'atac'        Input (N, 1, L): base-pair accessibility profile only.
                acc (1→n_filters) → dilated stack (n_filters channels).
                Uses n_filters (not n_acc_filters) for fair capacity comparison.

DeepLIFT/SHAP attributions on the multimodal input give:
    attr[:, :4, :] = sequence importance
    attr[:, 4:, :]  = accessibility importance
"""

import time
import numpy
import torch

from bpnetlite.losses import MNLLLoss, log1pMSELoss, _mixture_loss


def _asymmetric_mixture_loss(y, y_hat_logits, y_hat_logcounts, count_loss_weight,
                             overprediction_weight, labels=None, y_profile=None,
                             profile_loss_weight=1.0):
    """_mixture_loss with the count term made asymmetric.

    Identical to bpnetlite's `_mixture_loss` when overprediction_weight == 1.0 and
    y_profile is None -- verified numerically by scripts/test_asymmetric_loss.py, which is
    the regression gate for this file being shared with the p300 models.

    Above 1.0, squared log-count error is multiplied by that factor wherever the prediction
    EXCEEDS the truth. The endpoint (ABC ranking) is hurt far more by inventing signal at
    accessible-but-inactive elements than by missing a strong enhancer, and the symmetric
    loss encodes the opposite priority.

    y_profile MAKES THE TWO HEADS READ DIFFERENT TARGETS. Normally both terms score `y`:
    the profile head against its shape and the counts head against its sum. Given
    y_profile, the profile head is scored against that instead while the counts head still
    scores `y`. The point is an auxiliary task: H3K27ac's 1 bp inter-replicate ceiling is
    0.21, so the profile head trains against mostly noise and `profile_pearson` sits at
    0.063-0.064 whatever it is fed (F-013), whereas K562 DNase's is 0.848. Nothing about
    this changes what the counts head predicts, which is the quantity of interest.

    profile_loss_weight EXISTS BECAUSE MNLL SCALES WITH READ DEPTH. MNLL is
    -sum(y * log_softmax(logits)) plus a y-only term, so a deeper profile target multiplies
    the profile term and, since `loss = profile + count_loss_weight * count`, silently
    DOWN-weights the counts term. K562 DNase carries 17.8x the reads of K562 H3K27ac over
    the same 1 kb windows (910.2 against 51.1 mean), so swapping the target in unweighted
    would divide the effective count weight by about 18 and a loss on counts would be
    uninterpretable. Setting this to the depth ratio makes the profile term START at the
    magnitude the baseline arm had, so count_loss_weight keeps its meaning and the only
    thing that changed is what the profile head is pointed at.
    """
    y_hat_logits = y_hat_logits.reshape(y_hat_logits.shape[0], -1)
    y_hat_logits = torch.nn.functional.log_softmax(y_hat_logits, dim=-1)

    y = y.reshape(y.shape[0], -1)
    y_ = y.sum(dim=-1).reshape(y.shape[0], 1)

    # The profile term's target: y_profile when supplied, otherwise y itself.
    yp = y if y_profile is None else y_profile.reshape(y_profile.shape[0], -1)
    if y_profile is not None and yp.shape[-1] != y_hat_logits.shape[-1]:
        raise ValueError(
            f"profile target has {yp.shape[-1]} flattened positions but the profile head "
            f"emits {y_hat_logits.shape[-1]}; n_outputs must match the PROFILE target's "
            f"channel count, not the counts target's")

    if labels is not None:
        profile_loss = MNLLLoss(y_hat_logits[labels == 1], yp[labels == 1]).mean()
    else:
        profile_loss = MNLLLoss(y_hat_logits, yp).mean()

    # log1pMSELoss is mean((log_pred - log1p(true))^2); reproduce it elementwise so the
    # over-prediction cases can be weighted.
    err = y_hat_logcounts - torch.log(y_ + 1)
    w = torch.where(err > 0,
                    torch.full_like(err, float(overprediction_weight)),
                    torch.ones_like(err))
    count_loss = (w * err.pow(2)).mean()

    loss = profile_loss_weight * profile_loss + count_loss_weight * count_loss
    # profile_loss is returned UNWEIGHTED so the logged and compared value stays on the
    # same scale as every existing run's; only the gradient is reweighted.
    return profile_loss, count_loss, loss
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

    def __init__(self, n_filters=64, n_acc_filters=8, n_acc_channels=1,
                 n_layers=8, n_outputs=2,
                 mode='multimodal', count_loss_weight=1, profile_output_bias=True,
                 count_output_bias=True, name=None, trimming=None, verbose=True,
                 gate_accessibility=False, overprediction_weight=1.0,
                 profile_loss_weight=1.0):
        super().__init__()
        assert mode in ('multimodal', 'sequence', 'atac'), \
            f"mode must be 'multimodal', 'sequence', or 'atac', got '{mode}'"
        self.mode = mode
        self.n_filters = n_filters
        self.n_acc_filters = n_acc_filters
        self.n_layers = n_layers
        self.n_outputs = n_outputs
        # Number of accessibility INPUT channels (1 = a single flat track; >1 for
        # fragment-size-stratified channels). n_acc_filters is the output width.
        self.n_acc_channels = n_acc_channels
        self.gate_accessibility = gate_accessibility
        # >1 makes over-prediction cost more than under-prediction. log1pMSE is symmetric, so
        # predicting 0.80 where truth is 0.59 costs ~0.015 while missing an 18.6 enhancer
        # costs ~6.1 -- roughly 400x less sensitive to inventing signal than to missing it.
        # A gate has no gradient pressure to close without this.
        self.overprediction_weight = overprediction_weight
        self.count_loss_weight = count_loss_weight
        # Scales the profile term's gradient. Needed when the profile head reads a target
        # of a different read depth from the counts head; see _asymmetric_mixture_loss.
        self.profile_loss_weight = profile_loss_weight
        self.name = name or f"multimodal_bpnet.{mode}.{n_filters}.{n_layers}"
        self.trimming = trimming or 47 + sum(2**i for i in range(1, n_layers + 1))

        if mode == 'multimodal':
            n_merged = n_filters + n_acc_filters
            self.seq_conv = torch.nn.Conv1d(4, n_filters, kernel_size=21, padding=10)
            self.seq_relu = torch.nn.ReLU()
            self.acc_conv = torch.nn.Conv1d(n_acc_channels, n_acc_filters,
                                            kernel_size=21, padding=10)
            self.acc_relu = torch.nn.ReLU()
            if gate_accessibility:
                # Sequence decides, per position, how much accessibility to let through.
                # The trunk otherwise mixes the two branches ADDITIVELY, so accessibility
                # contributes equally everywhere -- which is why every model that sees ATAC
                # over-predicts H3K27ac at accessible-but-unacetylated elements (CpG-island
                # promoters, CTCF sites) by 7-8x while the sequence-only model elevates them
                # only 1.6x. The signal is already in sequence; it has no way to veto ATAC.
                #
                # Initialised OPEN: zero weights and bias +4 give sigmoid ~ 0.982 everywhere,
                # so at initialisation this model is the ungated one and can only learn to
                # close the gate where closing helps. That keeps the comparison clean and
                # avoids destabilising early training.
                self.gate_conv = torch.nn.Conv1d(n_filters, n_acc_filters,
                                                 kernel_size=21, padding=10)
                torch.nn.init.zeros_(self.gate_conv.weight)
                torch.nn.init.constant_(self.gate_conv.bias, 4.0)
        elif mode == 'sequence':
            n_merged = n_filters
            self.seq_conv = torch.nn.Conv1d(4, n_filters, kernel_size=21, padding=10)
            self.seq_relu = torch.nn.ReLU()
        elif mode == 'atac':
            n_merged = n_filters
            self.acc_conv = torch.nn.Conv1d(n_acc_channels, n_filters, kernel_size=21, padding=10)
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
        X: torch.Tensor
            Shape depends on mode:
              'multimodal': (N, 5, L) — channels 0-3 one-hot DNA, channel 4 accessibility
              'sequence':   (N, 4, L) — one-hot DNA only
              'atac':       (N, 1, L) — base-pair accessibility profile only

        Returns
        -------
        y_profile: torch.Tensor, shape (N, n_outputs, out_length)
        y_counts:  torch.Tensor, shape (N, 1)
        """
        start, end = self.trimming, X.shape[2] - self.trimming

        if self.mode == 'multimodal':
            X_seq = self.seq_relu(self.seq_conv(X[:, :4, :]))
            X_acc = self.acc_relu(self.acc_conv(X[:, 4:, :]))
            # getattr keeps checkpoints pickled before the gate existed loadable
            if getattr(self, "gate_accessibility", False):
                X_acc = X_acc * torch.sigmoid(self.gate_conv(X_seq))
            X_merged = torch.cat([X_seq, X_acc], dim=1)
        elif self.mode == 'sequence':
            X_merged = self.seq_relu(self.seq_conv(X))
        elif self.mode == 'atac':
            X_merged = self.acc_relu(self.acc_conv(X))

        for i in range(self.n_layers):
            X_conv = self.rrelus[i](self.rconvs[i](X_merged))
            X_merged = torch.add(X_merged, X_conv)

        y_profile = self.fconv(X_merged)[:, :, start:end]

        X_pooled = torch.mean(X_merged[:, :, start - 37:end + 37], dim=2)
        y_counts = self.linear(X_pooled).reshape(X.shape[0], 1)

        return y_profile, y_counts

    def fit(self, training_data, optimizer, scheduler=None, offset_valid=None,
            X_valid=None, y_valid=None, max_epochs=100, batch_size=64,
            dtype='float32', device='cuda', early_stopping=None,
            profile_target=False, y_profile_valid=None):
        """Train the model.

        Parameters
        ----------
        training_data: DataLoader
            Yields (X, y, labels) tuples where X has shape (N, 5, in_length).
            With profile_target, yields (X, y, labels, y_profile), and with residual
            offsets as well, (X, offset, y, labels, y_profile).
        optimizer: torch.optim.Optimizer
        scheduler: lr_scheduler or None
        X_valid: torch.Tensor or None, shape (n, 5, in_length)
        y_valid: torch.Tensor or None, shape (n, n_outputs, out_length)
        max_epochs: int
        batch_size: int
        dtype: str or torch.dtype
        device: str
        early_stopping: int or None
        profile_target: bool
            The loader supplies a separate target for the profile head. Passed
            EXPLICITLY rather than inferred from the batch, because the residual-offset
            layout already uses arity to signal itself and (X, y, labels, y_profile) and
            (X, offset, y, labels) are both length 4. Sniffing would silently read the
            label column as a target.
        y_profile_valid: torch.Tensor or None
            Validation profile target. Required when profile_target is set, since the
            validation loss must score the same objective as training.
        """
        if profile_target and y_profile_valid is None:
            raise ValueError("profile_target is set but y_profile_valid is None; the "
                             "validation loss would score a different objective from "
                             "training and early stopping would be meaningless")
        if y_profile_valid is not None and not profile_target:
            raise ValueError("y_profile_valid was given but profile_target is False")
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
                # Residual training: the dataset may insert a per-region count offset
                # at data[1]. It is ADDED to the predicted logcounts before the loss, so
                # the model learns (observed - offset) while the loss still scores the
                # real target.
                if profile_target:
                    X, y, labels, y_prof = data[0], data[-3], data[-2], data[-1]
                    offset = data[1].to(device).float() if len(data) > 4 else None
                    y_prof = y_prof.to(device)
                else:
                    X, y, labels = data[0], data[-2], data[-1]
                    offset = data[1].to(device).float() if len(data) > 3 else None
                    y_prof = None
                X = X.to(device).float()
                y = y.to(device)
                labels = labels.to(device)

                optimizer.zero_grad()
                self.train()

                with torch.autocast(device_type=device_type, dtype=dtype):
                    y_hat_logits, y_hat_logcounts = self(X)
                    if offset is not None:
                        y_hat_logcounts = y_hat_logcounts + offset.reshape(-1, 1)
                    train_profile_loss, train_count_loss, loss = _asymmetric_mixture_loss(
                        y, y_hat_logits, y_hat_logcounts,
                        self.count_loss_weight,
                        getattr(self, "overprediction_weight", 1.0), labels,
                        y_profile=y_prof,
                        profile_loss_weight=getattr(self, "profile_loss_weight", 1.0)
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
                if offset_valid is not None:
                    y_hat_logcounts = y_hat_logcounts + \
                        offset_valid.to(y_hat_logcounts.device).reshape(-1, 1)

                valid_profile_loss, valid_count_loss, valid_loss = _asymmetric_mixture_loss(
                    y_valid, y_hat_logits, y_hat_logcounts, self.count_loss_weight,
                    getattr(self, "overprediction_weight", 1.0),
                    y_profile=y_profile_valid,
                    profile_loss_weight=getattr(self, "profile_loss_weight", 1.0)
                )

                # Two calls when the heads read different targets: one y cannot serve both,
                # and logging profile_pearson against the counts target would report the
                # correlation of a head that was never trained on it.
                measures = calculate_performance_measures(
                    y_hat_logits,
                    y_valid if y_profile_valid is None else y_profile_valid,
                    y_hat_logcounts,
                    kernel_sigma=7, kernel_width=81,
                    measures=['profile_pearson', 'count_pearson']
                )
                valid_profile_corr = numpy.nan_to_num(measures['profile_pearson'])
                if y_profile_valid is None:
                    valid_count_corr = numpy.nan_to_num(measures['count_pearson'])
                else:
                    count_measures = calculate_performance_measures(
                        y_hat_logits, y_valid, y_hat_logcounts,
                        kernel_sigma=7, kernel_width=81,
                        measures=['count_pearson']
                    )
                    valid_count_corr = numpy.nan_to_num(count_measures['count_pearson'])
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
