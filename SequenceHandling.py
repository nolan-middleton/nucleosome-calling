##############################################################################
# This file contains library functions for genome manipulation and finding   #
# sequences within the genome.                                               #
##############################################################################

#%% Setup

# Imports
import numpy as np

#%% Miscellaneous Sequence Handling Functions

def reverse_complement(seq):
    '''
    This function will return the reverse complement of a DNA sequence.
    
    Arguments:
        seq: The sequence to take the reverse complement of.
    
    Note that seq must only contain UPPERCASE letters 'A', 'C', 'G', and 'T'.
    '''
    revSeq = seq[::-1]
    
    return revSeq.replace(
        "A",
        "t"
    ).replace(
        "C",
        "g"
    ).replace(
        "G",
        "c"
    ).replace(
        "T",
        "a"
    ).upper()

def reverse_complements(seqs):
    '''
    This function will return the reverse complement of an array of DNA
    sequences.
    
    Arguments:
        seqs: An array of sequences to take the reverse complement of.
    
    Note that the sequences must only contain UPPERCASE letters 'A', 'C', 'G',
    and 'T'. The sequences must also be of uniform length.
    '''
    if (len(np.unique(np.char.str_len(seqs))) > 1):
        raise Exception("seqs must all be uniform length!")
    
    returnSeqs = np.chararray.replace(seqs, "A", "t")
    returnSeqs = np.chararray.replace(returnSeqs, "C", "g")
    returnSeqs = np.chararray.replace(returnSeqs, "G", "c")
    returnSeqs = np.chararray.replace(returnSeqs, "T", "a")
    returnSeqs = np.chararray.replace(returnSeqs, "a", "A")
    returnSeqs = np.chararray.replace(returnSeqs, "c", "C")
    returnSeqs = np.chararray.replace(returnSeqs, "g", "G")
    returnSeqs = np.chararray.replace(returnSeqs, "t", "T")
    
    return np.array([seq[::-1] for seq in returnSeqs])

def reverse_dict_complement(seq_dict):
    '''
    This function will return the reverse complement of a dictionary of DNA
    sequences.
    
    Arguments:
        seq_dict: A dictionary of sequences to take the reverse complements of.
    
    Note that the sequences in seq_dict must only contain UPPERCASE letters
    'A', 'C', 'G', and 'T'.
    '''
    return_dict = {}
    for key in seq_dict:
        return_dict[key] = reverse_complement(seq_dict[key])
    
    return return_dict

def reverse_dict_complements(seqs_dict):
    '''
    This function will return the reverse complement of a dictionary of arrays
    of DNA sequences.
    
    Arguments:
        seqs_dict: A dictionary of numpy arrays of sequences to take the
        reverse complements of.
    
    Note that the sequences in seqs_dict must only contain UPPERCASE letters
    'A', 'C', 'G', and 'T'.
    '''
    return_dict = {}
    for key in seqs_dict:
        return_dict[key] = reverse_complements(seqs_dict[key])
    
    return return_dict

#%% Sequence Retrieval

def construct_n_mers(n, genome):
    '''
    This function will return a dictionary of n-mers from the genome.
    
    Arguments:
        n: The length of the n-mer to construct
        genome: The genome array.
    '''
    print("~~~ Constructing " + str(n) + "-mers ~~~")
    stretches = {}
    
    for chrom in genome:
        print("- " + chrom + "...")
        genLen = len(genome[chrom])
        stretches[chrom] = genome[chrom][:(genLen - n + 1)]
        for i in range(1, n):
            stretches[chrom] = np.char.add(
                stretches[chrom],
                genome[chrom][i:(genLen - n + 1 + i)]
            )
    
    print("~~~ Constructed " + str(n) + "-mers ~~~")
    return stretches

def find_sequences_occurrences(seqs, genome):
    '''
    This function will find all occurrences of seq in the genome, returning
    a dictionary where each key is a chromosome and each entry is a table
    containing the inclusive start and end locations of the sequence.
    
    Arguments:
        seqs: A numpy array of sequences to find in the genome.
        genome: The genome array.
    
    Note that this function will return two dictionaries, one for the minus
    strand, one for the plus strand. They're returned as (minus, plus).
    
    Also note that the input sequences must have a uniform length.
    '''
    if (len(seqs) == 0):
        raise Exception("You must supply some sequences!")
    
    print("~~~ Finding Sequence Occurrences ~~~")
    returnDictPlus = {}
    returnDictMinus = {}
    minusSeqs = reverse_complements(seqs)
    seqLen = len(seqs[0])
    
    for seq in seqs:
        if (len(seq) != seqLen):
            raise Exception("Sequences must be of uniform length!")
    
    stretches = construct_n_mers(seqLen, genome)
    
    for chrom in genome:
        starts = np.arange(0, len(stretches[chrom]))[
            np.isin(stretches[chrom], seqs)
        ]
        ends = starts + seqLen - 1
        returnDictPlus[chrom] = np.column_stack((starts, ends))
        
        starts = np.arange(0,len(stretches[chrom]))[
            np.isin(stretches[chrom], minusSeqs)
        ]
        ends = starts + seqLen - 1
        returnDictMinus[chrom] = np.column_stack((starts, ends))
    
    print("~~~ Found Sequence Occurrences ~~~")
    return (returnDictMinus, returnDictPlus)

def get_sequences_from_locations(starts, n, genome):
    '''
    This function will return the genomic sequence of n-mers originating from
    the starts.
    
    Arguments:
        starts: A dictionary of start locations
        n: The length of the n-mers to return
        genome: The genome array
    '''
    seqs = {}
    
    for chrom in starts:
        spans = np.repeat(
            starts[chrom][:,np.newaxis],
            n,
            axis = 1
        ) + np.repeat(
            np.arange(0, n)[np.newaxis, :],
            len(starts[chrom]),
            axis = 0
        )
        
        plusLetters = genome[chrom][spans]
        
        plusSeqs = plusLetters[:,0]
        for i in range(1,n):
            plusSeqs = np.char.add(plusSeqs, plusLetters[:,i])
        
        seqs[chrom] = plusSeqs
    
    return seqs

def get_sequences_from_intervals(intervals, genome):
    '''
    This function will return the genomic sequence corresponding to a given
    set of (potentially variable-length) intervals. This function is slow. If
    your intervals are all the same length, use get_sequences_from_locations
    instead.
    
    Arguments:
        intervals: A dictionary of intervals. Each key must be a chromosome
            and each value must be a numpy array with two columns: start
            (inclusive) and stop (exclusive), of course both zero-indexed.
        genome: The genome array.
    '''
    seqs = {}
    for chrom in intervals:
        N = np.shape(intervals[chrom])[0]
        seqs[chrom] = np.zeros(N, dtype = "<U32")
        for i in range(N):
            seqs[chrom][i] = "".join(
                genome[chrom][intervals[chrom][i,0]:intervals[chrom][i,1]]
            )
    
    return seqs

#%% Sequence Alignment

def global_alignment_score(seq1,seq2,score_function,gap_penalty,dtype=float):
    '''
    This function will perform the global sequence alignment algorithm on two
    sequences and return the score.
    
    Arguments:
        seq1: The first sequence to align.
        seq2: The second sequence to align.
        score_function: The score function to judge mismatches. Must take in
            two characters as its arguments.
        gap_penalty: The penalty for starting and extending a gap. Must be a
            POSITIVE number. The quantity used for this argument is SUBTRACTED.
        dtype: The data type to use during the alignment. score_function must
            return the same data type as specified here, and gap_penalty must
            the same data type as specified here. Defaults to float.
    '''
    H = np.zeros((len(seq1) + 1, len(seq2) + 1), dtype = dtype)
    H[:,0] = -gap_penalty*np.arange(0, np.shape(H)[0], dtype = dtype)
    H[0,:] = -gap_penalty*np.arange(0, np.shape(H)[1], dtype = dtype)
    
    for i in range(1, np.shape(H)[0]):
        for j in range(1, np.shape(H)[1]):
            H[i][j] = max(
                H[i-1][j-1] + score_function(seq1[i-1], seq2[j-1]),
                H[i-1][j] - gap_penalty,
                H[i][j-1] - gap_penalty
            )
    
    return H[-1,-1]

def local_alignment_score(seq1,seq2,score_function,gap_penalty,dtype=float):
    '''
    This function will perform the local sequence alignment algorithm on two
    sequences and return the score.
    
    Arguments:
        seq1: The first sequence to align.
        seq2: The second sequence to align.
        score_function: The score function to judge mismatches. Must take in
            two characters as its arguments.
        gap_penalty: The penalty for starting and extending a gap. Must be a
            POSITIVE number. The quantity used for this argument is SUBTRACTED.
        dtype: The data type to use during the alignment. score_function must
            return the same data type as specified here, and gap_penalty must
            the same data type as specified here. Defaults to float.
    '''
    H = np.zeros((len(seq1) + 1, len(seq2) + 1), dtype = dtype)
    
    for i in range(1, np.shape(H)[0]):
        for j in range(1, np.shape(H)[1]):
            H[i][j] = max(
                H[i-1][j-1] + score_function(seq1[i-1], seq2[j-1]),
                H[i-1][j] - gap_penalty,
                H[i][j-1] - gap_penalty,
                0
            )
    
    return np.max(H)