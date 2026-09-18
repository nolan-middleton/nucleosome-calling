#############################################################################
# This file will find overlapping nucleosomes between the greedy, Viterbi,  #
# Weiner, and Brogaard maps.                                                #
#############################################################################

print("===== 0d0_overlaps.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Defs
def triple_venn(A, B, C):
    a = 0
    b = 0
    c = 0
    ab = 0
    ac = 0
    bc = 0
    abc = 0
    
    for chrom in A:
        for loc in A[chrom]:
            in_b = True in (abs(B[chrom] - loc) <= 73)
            in_c = True in (abs(C[chrom] - loc) <= 73)
            
            if (in_b and in_c):
                abc += 1
            elif (in_b):
                ab += 1
            elif (in_c):
                ac += 1
            else:
                a += 1
    
    for chrom in B:
        for loc in B[chrom]:
            in_a = True in (abs(A[chrom] - loc) <= 73)
            in_c = True in (abs(C[chrom] - loc) <= 73)
            
            if (not in_a) and (in_c):
                bc += 1
            elif (not in_a):
                b += 1
    
    for chrom in C:
        c += len(C[chrom])
    
    c -= ac + bc + abc
    
    return (a, b, c, ab, ac, bc, abc)

#%% Load Data

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = an.make_type(io.load_tsv("PreprocessedData/blocks", chroms), int)

directory = "NucleosomePattern/composite_models/composite_ice"

greedy = io.load_tsv(
    directory + "/dyad_calls",
    chroms,
    dtype = int
)

greedy_107 = io.load_tsv(
    directory + "/dyad_calls_107bp",
    chroms,
    dtype = int
)

viterbi = io.load_tsv(
    directory + "/viterbi_dyad_calls",
    chroms,
    dtype = int
)

viterbi_overlap = io.load_tsv(
    directory + "/viterbi_40b_overlap_dyad_calls",
    chroms,
    dtype = int
)

viterbi_overlap2 = io.load_tsv(
    directory + "/viterbi_20b_overlap_dyad_calls",
    chroms,
    dtype = int
)

weinerChroms=io.load_list("PreprocessedData/OtherData/weiner/chromosomes.txt")
weiner = an.column_split(
    io.load_tsv("PreprocessedData/OtherData/weiner", weinerChroms)
)[2]
weiner = an.make_type(weiner, int)
order = an.dict_argsort(weiner)
weiner = an.dict_index(weiner, order)

brogaardChroms=io.load_list("PreprocessedData/all_nucleosomes/chromosomes.txt")
brogaard = an.make_type(
    an.column_split(
        io.load_tsv("PreprocessedData/all_nucleosomes",brogaardChroms)
    )[2],
    int
)

# Remove the first 7000bp of chrXVI, brogaard is missing that entire region
greedy["chrXVI"] = greedy["chrXVI"][greedy["chrXVI"] > 7000]
greedy_107["chrXVI"] = greedy_107["chrXVI"][greedy_107["chrXVI"] > 7000]
viterbi["chrXVI"] = viterbi["chrXVI"][viterbi["chrXVI"] > 7000]
viterbi_overlap["chrXVI"] = viterbi_overlap["chrXVI"][
    viterbi_overlap["chrXVI"] > 7000
]
viterbi_overlap2["chrXVI"] = viterbi_overlap2["chrXVI"][
    viterbi_overlap2["chrXVI"] > 7000
]
weiner["chrXVI"] = weiner["chrXVI"][weiner["chrXVI"] > 7000]
brogaard["chrXVI"] = brogaard["chrXVI"][brogaard["chrXVI"] > 7000]

#%% Main

greedy_weiner_brogaard = triple_venn(greedy, weiner, brogaard)
viterbi_weiner_brogaard = triple_venn(viterbi, weiner, brogaard)
viterbi20bp_weiner_brogaard = triple_venn(viterbi_overlap2, weiner, brogaard)

with open("VennDiagrams.txt", "w") as file:
    file.write("GreedyWeinerBrogaard: ")
    file.write("\t".join([str(a) for a in greedy_weiner_brogaard]) + "\n")
    file.write("ViterbiWeinerBrogaard: ")
    file.write("\t".join([str(a) for a in viterbi_weiner_brogaard]) + "\n")
    file.write("Viterbi20bpWeinerBrogaard: ")
    file.write("\t".join([str(a) for a in viterbi20bp_weiner_brogaard]) + "\n")

#%% End

print("===== 0d0_overlaps.py - Exiting Properly =====")