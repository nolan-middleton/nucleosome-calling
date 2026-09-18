#############################################################################
# This file will evaluate the nucleosome placement based on the center to   #
# center distances from the predicted nucleosomes to the real nucleosomes.  #
#############################################################################

print("===== 0c6_center_to_center_rotational.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

#%% Load Data

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = an.make_type(io.load_tsv("PreprocessedData/blocks", chroms), int)

directory = "NucleosomePattern/composite_models/composite_ice"

viterbi = io.load_tsv(
    directory + "/viterbi_dyad_calls",
    chroms,
    dtype = str
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

redund = np.loadtxt(
    "../Data/Brogaard_RedundantDyads.tsv",
    dtype = str,
    delimiter = "\t"
)
redundant = {}
for chrom in np.unique(redund[:,0]):
    redundant[str(chrom)] = redund[redund[:,0] == chrom, 1].astype(int)

rotational_scores = np.loadtxt(
    directory + "/rotational_scores.tsv",
    delimiter = "\t",
    dtype = str
)

TSS = io.load_tsv(
    "PreprocessedData/TSS",
    io.load_list("PreprocessedData/TSS/chromosomes.txt"),
    dtype = str
)

TSS_sites = an.make_type(an.column_split(TSS)[0], int)
TSS_strands = an.column_split(TSS)[4]

#%% Main

print("> Getting distances...")
rotational_scores = rotational_scores[
    np.argsort(rotational_scores[:,2].astype(float)),
    :
]
N = np.shape(rotational_scores)[0] // 10
P1 = rotational_scores[:N,:]
P9 = rotational_scores[-N:,:]

redund_dists = np.zeros((1001,3), dtype = int)
weiner_dists = np.zeros((1001,3), dtype = int)

redund_dists[:,0] = np.arange(-500, 501)
weiner_dists[:,0] = np.arange(-500, 501)

for i in range(N):
    chrom = str(P1[i,0])
    dyad = int(P1[i,1])
    
    dists = redundant[chrom] - dyad
    redund_dists[dists[abs(dists) <= 500] + 500,1] += 1
    
    dists = weiner[chrom] - dyad
    weiner_dists[dists[abs(dists) <= 500] + 500,1] += 1
    
    chrom = str(P9[i,0])
    dyad = int(P9[i,1])
    
    dists = redundant[chrom] - dyad
    redund_dists[dists[abs(dists) <= 500] + 500,2] += 1
    
    dists = weiner[chrom] - dyad
    weiner_dists[dists[abs(dists) <= 500] + 500,2] += 1

print("> Getting TSS...")
TSS = np.zeros((1151,3), dtype = int)
TSS[:,0] = np.arange(-500,651)
for chrom in TSS_sites:
    for i in range(len(TSS_sites[chrom])):
        site = TSS_sites[chrom][i]
        strand = TSS_strands[chrom][i]
        
        theseP1 = P1[P1[:,0] == chrom,1].astype(int)
        dists = theseP1 - site
        if (strand == "-"):
            dists = -dists
        
        TSS[dists[(dists >= -500) & (dists <= 650)] + 500, 1] += 1
        
        theseP9 = P9[P9[:,0] == chrom,1].astype(int)
        dists = theseP9 - site
        if (strand == "-"):
            dists = -dists
        
        TSS[dists[(dists >= -500) & (dists <= 650)] + 500, 2] += 1

#%% Outputs

print("> Outputting...")
if (not os.path.isdir(directory + "/rotational_center_to_center")):
    os.mkdir(directory + "/rotational_center_to_center")

np.savetxt(
    directory + "/rotational_center_to_center/redundant.tsv",
    redund_dists,
    delimiter = "\t",
    fmt = "%d"
)
np.savetxt(
    directory + "/rotational_center_to_center/weiner.tsv",
    weiner_dists,
    delimiter = "\t",
    fmt = "%d"
)
np.savetxt(
    directory + "/rotational_center_to_center/TSS.tsv",
    TSS,
    delimiter = "\t",
    fmt = "%d"
)

#%% End

print("===== 0c6_center_to_center_rotational.py - Exiting Properly =====")