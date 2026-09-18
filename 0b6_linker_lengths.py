#############################################################################
# This file will find all of the linker lengths of nucleosome placements.   #
#############################################################################

print("===== 0b6_linker_lengths.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

#%% Load Data

print("> Loading data...")
dirs = [
    "NucleosomePattern/cl_rt_log2",
    "NucleosomePattern/cl_ice_log2",
    "NucleosomePattern/cl_ice_abs_diff",
    "NucleosomePattern/composite_models/composite_ice"
]

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

# Load in test nucleosome placements
print(">> Putative dyad calls...")
dyads = [
    io.load_tsv(
        D + "/viterbi_dyad_calls",
        chromosomes = io.load_list(D + "/viterbi_dyad_calls/chromosomes.txt"),
        dtype = int
    )
    for D in dirs
]

dyads += [
    io.load_tsv(
        D + "/dyad_calls",
        chromosomes=io.load_list(D+"/dyad_calls/chromosomes.txt"),
        dtype = int
    )
    for D in dirs
]

dirs_greedy = [d + "/GreedyLinkerLens" for d in dirs]
for d in dirs_greedy:
    if (not os.path.isdir(d)):
        os.mkdir(d)

dirs += dirs_greedy


print(">> Genome and blocks...")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)
for chrom in blocks:
    blocks[chrom] = blocks[chrom].astype(int)

print(">> Nucleosomes...")
all_nucleosomes = io.load_tsv("PreprocessedData/all_nucleosomes", chroms)
weiner = io.load_tsv("PreprocessedData/OtherData/weiner", chroms)

dyads += [
    an.make_type(an.column_split(all_nucleosomes)[2], int),
    an.make_type(an.column_split(weiner)[2], int)
]

if (not os.path.isdir("NucleosomePattern/ReferenceLinkerLengths")):
    os.mkdir("NucleosomePattern/ReferenceLinkerLengths")

if (not os.path.isdir("NucleosomePattern/ReferenceLinkerLengths/all")):
    os.mkdir("NucleosomePattern/ReferenceLinkerLengths/all")

if (not os.path.isdir("NucleosomePattern/ReferenceLinkerLengths/weiner")):
    os.mkdir("NucleosomePattern/ReferenceLinkerLengths/weiner")

dirs += [
    "NucleosomePattern/ReferenceLinkerLengths/all",
    "NucleosomePattern/ReferenceLinkerLengths/weiner"
]

#%% Main

for i in range(len(dirs)):
    print("> " + dirs[i] + "...")
    linker_lens = {}
    for chrom in dyads[i]:
        S = np.sort(dyads[i][chrom])
        linker_lens[chrom] = S[1:] - S[:-1] - 146 - 1
    io.output_tsv(linker_lens, dirs[i]+"/linker_lens", fmt="%d")
    np.savetxt(
        dirs[i] + "/linker_lens.tsv",
        an.dict_concat(linker_lens),
        delimiter = "\t",
        fmt = "%d"
    )

#%% End

print("===== 0b6_linker_lengths.py - Exiting Properly =====")