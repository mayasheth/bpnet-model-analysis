import tensorflow as tf
import bpnet.model.arch
from bpnet.model.losses import multinomial_nll
import time, datetime, json, os, sys
import pandas as pd
import numpy as np
import pysam
import random
from scipy.special import logsumexp

def get_model(model_path):
    if model_path.endswith(".h5"):
        # ChromBPNet model
        custom_objects={"multinomial_nll":multinomial_nll, "tf": tf}    
        tf.keras.utils.get_custom_objects().update(custom_objects)    
        model = tf.keras.models.load_model(model_path, compile=False)
        return "ChromBPNet", model
    else:
        # assume BPNet model
        model = tf.keras.models.load_model(model_path)
        return "BPNet", model

def get_shuffled_peak_sequences(peak_path, fasta_path, input_seq_len = 2114, num_of_simulated_sequences=1000):
    peaks_df = pd.read_csv(peak_path, sep='\t',header=None,
        names=['chrom', 'start', 'end', 'name', 'score', 'strand', 'signalValue', 'p', 'q', 'summit'])

    peaks_df = peaks_df.sample(frac=1).reset_index(drop=True)
    
    valid_chroms = ['chr2', 'chr11', 'chr6', 'chr19', 'chr5', 'chr18',
       'chrX', 'chr1', 'chr10', 'chr12', 'chr15', 'chr22', 'chr17',
       'chr3', 'chr4', 'chr20', 'chr21', 'chr14', 'chr16', 'chr9', 'chr7',
       'chr8', 'chr13', 'chrY']
    
    peaks_df = peaks_df.loc[peaks_df['chrom'].isin(valid_chroms),:].reset_index(drop=True)

    fasta_ref = pysam.FastaFile(fasta_path)

    sequences = []

    input_start = peaks_df['start'] + peaks_df['summit'] - (input_seq_len//2)
    input_start = input_start.mask(input_start < 0,0)
    input_end = peaks_df['start'] + peaks_df['summit'] + (input_seq_len//2)


    for i in range(num_of_simulated_sequences):
        actual_sequence = fasta_ref.fetch(peaks_df['chrom'][i], input_start[i] , input_end[i]).upper()
        if (input_seq_len-len(actual_sequence))>0:
            print(peaks_df['chrom'][i], input_start[i] , input_end[i])
        padded_sequence = actual_sequence + (random_seq(input_seq_len-len(actual_sequence)))
        shuffled_seq = dinuc_shuffle(padded_sequence)
        sequences.append(shuffled_seq)
    
    return sequences

def reverse_complement(seq):
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement[base] for base in reversed(seq))

def dinuc_shuffle(seq, num_shufs=None, rng=None):
    """
    Creates shuffles of the given sequence, in which dinucleotide frequencies
    are preserved.
    Arguments:
        `seq`: either a string of length L, or an L x D NumPy array of one-hot
            encodings
        `num_shufs`: the number of shuffles to create, N; if unspecified, only
            one shuffle will be created
        `rng`: a NumPy RandomState object, to use for performing shuffles
    If `seq` is a string, returns a list of N strings of length L, each one
    being a shuffled version of `seq`. If `seq` is a 2D NumPy array, then the
    result is an N x L x D NumPy array of shuffled versions of `seq`, also
    one-hot encoded. If `num_shufs` is not specified, then the first dimension
    of N will not be present (i.e. a single string will be returned, or an L x D
    array).
    """
    if type(seq) is str:
        arr = string_to_char_array(seq)
    elif type(seq) is np.ndarray and len(seq.shape) == 2:
        seq_len, one_hot_dim = seq.shape
        arr = one_hot_to_tokens(seq)
    else:
        raise ValueError("Expected string or one-hot encoded array")

    if not rng:
        rng = np.random.RandomState()
   
    # Get the set of all characters, and a mapping of which positions have which
    # characters; use `tokens`, which are integer representations of the
    # original characters
    chars, tokens = np.unique(arr, return_inverse=True)

    # For each token, get a list of indices of all the tokens that come after it
    shuf_next_inds = []
    for t in range(len(chars)):
        mask = tokens[:-1] == t  # Excluding last char
        inds = np.where(mask)[0]
        shuf_next_inds.append(inds + 1)  # Add 1 for next token
 
    if type(seq) is str:
        all_results = []
    else:
        all_results = np.empty(
            (num_shufs if num_shufs else 1, seq_len, one_hot_dim),
            dtype=seq.dtype
        )

    for i in range(num_shufs if num_shufs else 1):
        # Shuffle the next indices
        for t in range(len(chars)):
            inds = np.arange(len(shuf_next_inds[t]))
            inds[:-1] = rng.permutation(len(inds) - 1)  # Keep last index same
            shuf_next_inds[t] = shuf_next_inds[t][inds]

        counters = [0] * len(chars)
       
        # Build the resulting array
        ind = 0
        result = np.empty_like(tokens)
        result[0] = tokens[ind]
        for j in range(1, len(tokens)):
            t = tokens[ind]
            ind = shuf_next_inds[t][counters[t]]
            counters[t] += 1
            result[j] = tokens[ind]

        if type(seq) is str:
            all_results.append(char_array_to_string(chars[result]))
        else:
            all_results[i] = tokens_to_one_hot(chars[result], one_hot_dim)
    return all_results if num_shufs else all_results[0]

def string_to_char_array(seq):
    """
    Converts an ASCII string to a NumPy array of byte-long ASCII codes.
    e.g. "ACGT" becomes [65, 67, 71, 84].
    """
    return np.frombuffer(bytearray(seq, "utf8"), dtype=np.int8)

def char_array_to_string(arr):
    """
    Converts a NumPy array of byte-long ASCII codes into an ASCII string.
    e.g. [65, 67, 71, 84] becomes "ACGT".
    """
    return arr.tobytes().decode("ascii")

def one_hot_to_tokens(one_hot):
    """
    Converts an L x D one-hot encoding into an L-vector of integers in the range
    [0, D], where the token D is used when the one-hot encoding is all 0. This
    assumes that the one-hot encoding is well-formed, with at most one 1 in each
    column (and 0s elsewhere).
    """
    tokens = np.tile(one_hot.shape[1], one_hot.shape[0])  # Vector of all D
    seq_inds, dim_inds = np.where(one_hot)
    tokens[seq_inds] = dim_inds
    return tokens

def tokens_to_one_hot(tokens, one_hot_dim):
    """
    Converts an L-vector of integers in the range [0, D] to an L x D one-hot
    encoding. The value `D` must be provided as `one_hot_dim`. A token of D
    means the one-hot encoding is all 0s.
    """
    identity = np.identity(one_hot_dim + 1)[:, :-1]  # Last row is all 0s
    return identity[tokens]

def random_seq(seqlen):
    return ''.join(random.choices("ACGT", k=seqlen))

def fix_sequence_length(sequence, length):
    """
        Function to check if length of sequence matches specified
        length and then return a sequence that's either padded or
        truncated to match the given length
        Args:
            sequence (str): the input sequence
            length (int): expected length
        Returns:
            str: string of length 'length'
    """

    # check if the sequence is smaller than expected length
    if len(sequence) < length:
        # pad the sequence with 'N's
        sequence += 'N' * (length - len(sequence))
    # check if the sequence is larger than expected length
    elif len(sequence) > length:
        # truncate to expected length
        sequence = sequence[:length]

    return sequence

def one_hot_encode(sequences, seq_length=2114):

    """
    
       One hot encoding of a list of DNA sequences 
       
       Args:
           sequences (list): python list of strings of equal length
           seq_length (int): expected length of each sequence in the 
               list
           
       Returns:
           numpy.ndarray: 
               3-dimension numpy array with shape 
               (len(sequences), len(list_item), 4)
    """
    
    if len(sequences) == 0:
        raise ValueError("'sequences; is empty!")
    
    # First, let's make sure all sequences are of equal length
    sequences = list(map(
        fix_sequence_length, sequences, [seq_length] * len(sequences)))

    # Step 1. convert sequence list into a single string
    _sequences = ''.join(sequences)
    
    # Step 2. translate the alphabet to a string of digits
    transtab = str.maketrans('ACGTNYRMSWKU', '012344444444')    
    sequences_trans = _sequences.translate(transtab)
    
    # Step 3. convert to list of ints
    int_sequences = list(map(int, sequences_trans))
    
    # Step 4. one hot encode using int_sequences to index 
    # into an 'encoder' array
    encoder = np.vstack([np.eye(4), np.zeros(4)])
    X = encoder[int_sequences]

    # Step 5. reshape 
    return X.reshape(len(sequences), len(sequences[0]), 4)

def one_hot_decode(one_hot_array):
    """
    Convert one-hot encoded DNA sequences back to string sequences.

    Args:
        one_hot_array (np.ndarray): shape (N, L, 4) — N sequences, each length L

    Returns:
        list of str: DNA sequences as strings
    """
    # Index of max value across the 4 channels (A, C, G, T)
    base_indices = np.argmax(one_hot_array, axis=-1)  # shape: (N, L)

    # Map index to base: 0 → 'A', 1 → 'C', 2 → 'G', 3 → 'T'
    index_to_base = np.array(['A', 'C', 'G', 'T'])

    # Use advanced indexing to convert indices to characters
    base_chars = index_to_base[base_indices]  # shape: (N, L)

    # Join across axis 1 to make strings: ['ACGT...', 'TGC...']
    sequences = np.char.join('', base_chars)

    return sequences.tolist()

def vectorized_prediction_to_profile(predictions,output_len=1000):
        logits_arr = predictions[0]
        counts_arr = predictions[1]
        pred_profile_logits = np.reshape(logits_arr,[-1,1,output_len*2])
        probVals_array = np.exp(pred_profile_logits-logsumexp(pred_profile_logits,axis=2).reshape([len(logits_arr),1,1]))
        profile_predictions = np.multiply((np.exp(counts_arr)).reshape([len(counts_arr),1,1]),probVals_array)
        plus = np.reshape(profile_predictions,[len(counts_arr),output_len,2])[:,:,0]
        minus = np.reshape(profile_predictions,[len(counts_arr),output_len,2])[:,:,1]
        return(plus,minus)

def make_model_prediction(mod_encoded, model_info, output_seq_len=1000, number_of_strands=2, batch_size=1024):
    """
    Universal wrapper that calls the appropriate prediction function based on model type
    
    Args:
        mod_encoded: input sequences
        model_info: tuple of (model_type, model) from get_model()
        other args: same as original function
    """
    model_type, model = model_info
    
    if model_type == "BPNet":
        return make_bpnet_model_prediction(mod_encoded, model, output_seq_len, number_of_strands, batch_size)
    elif model_type == "ChromBPNet":
        return make_chrombpnet_model_prediction(mod_encoded, model, output_seq_len, number_of_strands, batch_size)
    else:
        raise ValueError(f"Unknown model type: {model_type}")
    
def make_bpnet_model_prediction(mod_encoded, model, output_seq_len=1000, number_of_strands=2, batch_size=1024):
    n = mod_encoded.shape[0]
    preds_all = []

    # Preallocate max-size bias arrays once; slice per batch
    prof_bias = np.zeros((batch_size, output_seq_len, number_of_strands), dtype='float32')
    cnt_bias  = np.zeros((batch_size, number_of_strands), dtype='float32')

    serve = model.signatures['serving_default']

    for i in range(0, n, batch_size):
        bs = min(batch_size, n - i)
        batch = mod_encoded[i:i + bs].astype('float32')
        prof_b = prof_bias[:bs]
        cnt_b  = cnt_bias[:bs]

        pred_fwd = serve(**{'profile_bias_input_0': prof_b, 'counts_bias_input_0': cnt_b, 'sequence': batch})
        pred_rev = serve(**{'profile_bias_input_0': prof_b, 'counts_bias_input_0': cnt_b, 'sequence': batch[:, ::-1, ::-1]})

        logcounts_avg = np.mean([
            pred_fwd['logcounts_predictions'].numpy().squeeze(),
            pred_rev['logcounts_predictions'].numpy().squeeze()
        ], axis=0)

        preds_all.append(logcounts_avg)

    return np.concatenate(preds_all, axis=0)

def make_chrombpnet_model_prediction(mod_encoded, model, output_seq_len=1000, number_of_strands=2, batch_size=1024):
    n = mod_encoded.shape[0]
    preds_all = []

    for i in range(0, n, batch_size):
        bs = min(batch_size, n - i)
        batch = mod_encoded[i:i + bs].astype('float32')
        
        # Call model using signature with sequence keyword
        pred_profile, pred_logcounts = model.predict(batch) # shapes (N, 1000) and (N, 1)
        pred_logcounts = pred_logcounts.squeeze()

        preds_all.append(pred_logcounts)
    
    return np.concatenate(preds_all, axis=0)


# finds motif pairs in actual peaks?
def motif_pair_dfi(dfi_filtered, motif_pair):
    """Construct the matrix of motif pairs

    Args:
      dfi_filtered: dfi filtered to the desired property
      motif_pair: tuple of two pattern_name's
    Returns:
      pd.DataFrame with columns from dfi_filtered with _x and _y suffix
    """

    ## Merge motif instances: Finds all pairs of the specified motifs in the same sequence
    dfa = dfi_filtered[dfi_filtered.pattern_name == motif_pair[0]]
    dfb = dfi_filtered[dfi_filtered.pattern_name == motif_pair[1]]

    dfab = pd.merge(dfa, dfb, on='example_idx', how='outer')
    dfab = dfab[~dfab[['pattern_center_x', 'pattern_center_y']].isnull().any(1)]

    ## Calculate spacing
    dfab['center_diff'] = dfab.pattern_center_y - dfab.pattern_center_x

    # Handles directionality
    if "pattern_center_aln_x" in dfab:
        dfab['center_diff_aln'] = dfab.pattern_center_aln_y - dfab.pattern_center_aln_x
    dfab['strand_combination'] = dfab.strand_x + dfab.strand_y
    # Assure the right strand combination
    dfab.loc[dfab.center_diff < 0, 'strand_combination'] = dfab[dfab.center_diff < 0]['strand_combination'].map(comp_strand_combination).values

    # Special handling for identical motifs
    if motif_pair[0] == motif_pair[1]:
        dfab.loc[dfab['strand_combination'] == "--", 'strand_combination'] = "++"
        dfab = dfab[dfab.center_diff > 0]
    else:
        dfab.center_diff = np.abs(dfab.center_diff)
        if "center_diff_aln" in dfab:
            dfab.center_diff_aln = np.abs(dfab.center_diff_aln)
    if "center_diff_aln" in dfab:
        dfab = dfab[dfab.center_diff_aln != 0]  # exclude perfect matches
    return dfab

### --- FROM CLAUDE --- ###
def generate_orientation_patterns(n_motifs):
    """
    Generate all possible orientation patterns for n motifs
    Returns list of tuples like ('+', '+', '-') for n=3
    """
    if n_motifs == 1:
        return [('+',)]  # Single motif, only forward orientation needed
    
    # For multiple motifs, generate all combinations
    from itertools import product
    return list(product(['+', '-'], repeat=n_motifs))

def pattern_to_string(pattern):
    """Convert pattern tuple to string representation"""
    return ''.join(pattern)

def get_representative_patterns(n_motifs):
    """
    Get biologically meaningful orientation patterns to test
    This reduces the combinatorial explosion while keeping key patterns
    """
    if n_motifs == 1:
        return [('+',)]
    elif n_motifs == 2:
        return [('+', '+'), ('+', '-'), ('-', '+')]  # Skip '--' as it's equivalent to '++'
    elif n_motifs == 3:
        return [
            ('+', '+', '+'),    # All same orientation
            ('+', '-', '+'),    # Alternating
            ('+', '+', '-')    # Two same, one different (same as +,-,-)
        ]
    elif n_motifs == 4:
        return [
            ('+', '+', '+', '+'),  # All same
            ('+', '-', '+', '-'),  # Alternating
            ('+', '+', '-', '-'),  # Two blocks
            ('+', '-', '-', '+'),  # Symmetric
        ]
    else:
        # For n > 4, use a subset of patterns to avoid explosion
        return [
            tuple(['+'] * n_motifs),  # All forward
            tuple(['+', '-'] * (n_motifs // 2) + ['+'] * (n_motifs % 2)),  # Alternating
            tuple(['+'] * (n_motifs // 2) + ['-'] * (n_motifs - n_motifs // 2)),  # Half and half
        ]
    
# get all motif patterns (more intensive, esp for larger n)
def get_all_patterns(n_motifs):
    """Get all 2^n orientation patterns"""
    from itertools import product
    return list(product(['+', '-'], repeat=n_motifs))

# to test specific pattern
def get_hypothesis_patterns(n_motifs):
    """Test specific biological hypotheses"""
    patterns = []
    patterns.append(tuple(['+'] * n_motifs))  # All cooperative
    if n_motifs >= 2:
        patterns.append(tuple(['+', '-'] * (n_motifs // 2) + ['+'] * (n_motifs % 2)))  # Alternating
    return patterns

def insert_motifs_with_orientation_general(bg_encoded, n, spacing, motif_seq, orientation_pattern, seq_len=2114):
    """
    Insert n motifs with specified spacing and orientation pattern
    
    Args:
        bg_encoded: Background sequences (N, seq_len, 4)
        n: Number of motifs to insert
        spacing: Spacing between motifs (bp)
        motif_seq: List of motif sequence strings, same length as orientation_pattern (motifs can have different lengths)
        orientation_pattern: Tuple like ('+', '-', '+') specifying orientation for each motif
        seq_len: Sequence length
    
    Returns:
        Modified sequences with motifs inserted
    """
    num_bg, seq_len, _ = bg_encoded.shape
    
    # Validate inputs
    if len(motif_seq) != len(orientation_pattern):
        raise ValueError(f"motif_seq length ({len(motif_seq)}) must match orientation_pattern length ({len(orientation_pattern)})")
    
    if len(orientation_pattern) != n:
        raise ValueError(f"orientation_pattern length ({len(orientation_pattern)}) must match n ({n})")
    
    # Get motif lengths (can be different)
    motif_lengths = [len(seq) for seq in motif_seq]
    
    # Calculate total insertion length
    total_motif_len = sum(motif_lengths)
    insert_len = total_motif_len + ((n - 1) * spacing)
    insert_center = seq_len // 2
    start = insert_center - (insert_len // 2)
    
    mod_encoded = bg_encoded.copy()
    
    # Generate motif sequences based on orientation pattern and handle N's
    motif_sequences = []
    
    for i, orientation in enumerate(orientation_pattern):
        base_seq = motif_seq[i]
        if orientation == '+':
            oriented_seq = base_seq
        else:  # orientation == '-'
            oriented_seq = reverse_complement(base_seq)
        
        # Handle N's in the sequence by replacing with background sequence
        if 'N' in oriented_seq:
            # We'll handle N replacement per background sequence during insertion
            motif_sequences.append(oriented_seq)
        else:
            motif_sequences.append(oriented_seq)
    
    # Insert each motif
    current_pos = start
    for i in range(n):
        motif_len = motif_lengths[i]
        if 0 <= current_pos < seq_len - motif_len:
            motif_seq_to_insert = motif_sequences[i]
            
            # Handle N's by replacing with background sequence at that position
            if 'N' in motif_seq_to_insert:
                # For each background sequence, replace N's with the original background nucleotide
                for bg_idx in range(num_bg):
                    # Convert background one-hot to sequence for this position range
                    bg_slice = bg_encoded[bg_idx, current_pos:current_pos+motif_len, :][np.newaxis, ...]
                    bg_seq = one_hot_decode(bg_slice)[0]

                    # Replace N's with background nucleotides
                    final_seq = ""
                    for pos, nucleotide in enumerate(motif_seq_to_insert):
                        if nucleotide == 'N':
                            final_seq += bg_seq[pos]
                        else:
                            final_seq += nucleotide

                    # One-hot encode the final sequence and insert
                    motif_oh = one_hot_encode([final_seq], motif_len)[0]
                    mod_encoded[bg_idx, current_pos:current_pos+motif_len, :] = motif_oh
            else:
                # No N's, insert the same motif for all background sequences
                motif_oh = one_hot_encode([motif_seq_to_insert], motif_len)[0]
                mod_encoded[:, current_pos:current_pos+motif_len, :] = motif_oh
        
        # Move to next position (current motif length + spacing)
        current_pos += motif_len + spacing
            
    return mod_encoded

# for motif pairs
def generate_motif_pairs(motif_dict, include_self_pairs=True):
    """
    Generate all unique pairs of motifs from motif_dict
    
    Args:
        motif_dict: Dictionary of {motif_name: motif_seq}
        include_self_pairs: Whether to include (A, A) type pairs
    
    Returns:
        List of tuples: [(name1, name2, seq1, seq2), ...]
    """
    from itertools import combinations_with_replacement, combinations
    
    motif_items = list(motif_dict.items())
    
    if include_self_pairs:
        # Use combinations_with_replacement for (A,A), (A,B), (B,B), etc.
        pairs = combinations_with_replacement(motif_items, 2)
    else:
        # Use combinations for only (A,B) pairs
        pairs = combinations(motif_items, 2)
    
    result = []
    for (name1, seq1), (name2, seq2) in pairs:
        result.append((name1, name2, seq1, seq2))
    
    return result