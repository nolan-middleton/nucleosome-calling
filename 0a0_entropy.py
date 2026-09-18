#############################################################################
# This file will evaluate the nucleosome placement based on the probability #
# distribution and the entropy to judge "fuzziness" and how strongly the    #
# nucleosomes are positioned.                                               #
#############################################################################

print("===== 0a0_entropy.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an

#%% Load Data

print("> Loading data...")
dirs = [
    "NucleosomePattern/cl_rt_log2",
    "NucleosomePattern/cl_ice_log2",
    "NucleosomePattern/cl_ice_abs_diff",
    "NucleosomePattern/composite_models/composite_ice"
]
num_dirs = len(dirs)

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

# Load in test nucleosome placements
print(">> Putative dyad calls...")
putative_dyads = [
    io.load_tsv(
        D + "/viterbi_dyad_calls",
        chromosomes = io.load_list(D + "/viterbi_dyad_calls/chromosomes.txt"),
        dtype = int
    )
    for D in dirs
]

print(">> Probabilities...")
probs = [
    io.load_tsv(
        D + "/dynamic_probability_distribution",
        chromosomes = io.load_list(
            D + "/dynamic_probability_distribution/chromosomes.txt"
        )
    )
    for D in dirs
]

print(">> Genome and blocks...")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)
for chrom in blocks:
    blocks[chrom] = blocks[chrom].astype(int)

#%% Main

H = []
for i in range(num_dirs):
    print("> " + dirs[i] + "...")
    H.append({})
    for chrom in chroms:
        theseDyads = putative_dyads[i][chrom]
        theseProbs = probs[i][chrom]
        
        H[-1][chrom] = np.zeros(len(theseDyads))
        for j in range(len(theseDyads)):
            d = theseDyads[j]
            P = theseProbs[(theseProbs[:,0]>=d-73)&(theseProbs[:,0]<=d+73), 1]
            P /= np.sum(P)
            
            H[-1][chrom][j] = -np.sum(P*np.log2(P))

#%% Save

for i in range(num_dirs):
    H[i] = an.column_join([putative_dyads[i], H[i]])
    io.output_tsv(H[i], dirs[i] + "/entropy")
    io.output_wig(H[i], dirs[i] + "/entropy.wig")

#%% Main 2

H = []
for i in range(num_dirs):
    print("> " + dirs[i] + "...")
    H.append({})
    for chrom in chroms:
        theseDyads = putative_dyads[i][chrom]
        theseProbs = probs[i][chrom]
        
        H[-1][chrom] = np.zeros(len(theseDyads))
        for j in range(len(theseDyads)):
            d = theseDyads[j]
            P = theseProbs[(theseProbs[:,0]>=d-73)&(theseProbs[:,0]<=d+73), 1]
            P /= len(P)
            
            H[-1][chrom][j] = -np.sum(P*np.log2(P)) - \
                (1-np.sum(P))*np.log2(1-np.sum(P))

#%% Save 2

for i in range(num_dirs):
    H[i] = an.column_join([putative_dyads[i], H[i]])
    io.output_tsv(H[i], dirs[i] + "/entropy_withNo")
    io.output_wig(H[i], dirs[i] + "/entropy_withNo.wig")

#%% End

print("===== 0a0_entropy.py - Exiting Properly =====")