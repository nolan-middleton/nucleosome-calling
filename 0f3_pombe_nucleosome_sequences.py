#############################################################################
# This file will find sequences pertaining to the nucleosomes identied in   #
# S. pombe by Moyle-Heyrman et al.                                          #
#############################################################################

print("===== 0f3_pombe_nucleosome_sequences.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Defs
def do_sequence_analysis(dyads, genome, directory):
    positions = {}
    seqs = {}
    for chrom in dyads:
        positions[chrom] = np.zeros((len(dyads[chrom]), 146), float)
        for i in range(len(dyads[chrom])):
            positions[chrom][i,:] = np.arange(
                dyads[chrom][i] - 72.5,
                dyads[chrom][i] + 73
            )
        lefts = (positions[chrom] - 0.5).astype(int)
        seqs[chrom] = genome[chrom][np.column_stack((lefts, lefts[:,-1] + 1))]
        seqs[chrom] = np.column_stack((dyads[chrom], seqs[chrom]))
    
    if (not os.path.isdir(directory)):
        os.mkdir(directory)
    
    io.output_tsv(seqs,directory+"/seqs",fmt="%s",chroms=True)
    
    return seqs

#%% Load Data

print("> Loading data...")
genome = {}
fa = io.read_fasta("../Data/ASM294v2.fa")
for i in range(np.shape(fa)[0]):
    if (str(fa[i,0]) != "chrM"):
        genome[str(fa[i,0])] = np.array(list(str(fa[i,1]).upper()))

fa = None # To clear it out of memory...
chroms = list(genome.keys())

U = np.loadtxt(
    "../Data/sd01.txt",
    dtype = str,
    usecols = [0,1]
)
unique = {}
for chrom in np.unique(U[:,0]):
    unique[str(chrom)] = U[U[:,0] == chrom,1].astype(int) - 1

R = np.loadtxt(
    "../Data/sd02.txt",
    dtype = str,
    usecols = [0,1]
)
redundant = {}
for chrom in np.unique(R[:,0]):
    redundant[str(chrom)] = R[R[:,0] == chrom,1].astype(int) - 1

for chrom in chroms:
    redundant[chrom] = redundant[chrom][
        (redundant[chrom] >= 73) & (redundant[chrom] < len(genome[chrom])-73)
    ]
    unique[chrom] = unique[chrom][
        (unique[chrom] >= 73) & (unique[chrom] < len(genome[chrom])-73)
    ]

#%% Get Sequences in All Nucleosomes

D = "pombe_sequence_analysis"
if (not os.path.isdir(D)):
    os.mkdir(D)

unique_seqs = do_sequence_analysis(unique, genome, D + "/unique")
redundant_seqs = do_sequence_analysis(redundant, genome, D + "/redundant")

#%% End

print("===== 0f3_pombe_nucleosome_sequences.py - Exiting Properly =====")