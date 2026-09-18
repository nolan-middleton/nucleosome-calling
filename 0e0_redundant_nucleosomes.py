#############################################################################
# This file will find overlapping nucleosomes between the greedy, Viterbi,  #
# Weiner, and Brogaard maps.                                                #
#############################################################################

print("===== 0e0_redundant_nucleosomes.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Defs
def getDist(dyads, reference):
    freqs = np.zeros(201, dtype = int)
    for chrom in reference:
        for i in range(len(reference[chrom])):
            theseDists = dyads[chrom] - reference[chrom][i]
            theseDists = theseDists[(theseDists >= -100) & (theseDists <= 100)]
            if (len(theseDists) > 0):
                freqs[theseDists + 100] += 1
    return freqs

def getOffset(dyads, reference):
    freqs = np.zeros(201, dtype = int)
    for chrom in dyads:
        for i in range(len(dyads[chrom])):
            theseDists = dyads[chrom][i] - reference[chrom]
            theseDists = theseDists[(theseDists >= -100) & (theseDists <= 100)]
            if (len(theseDists) > 0):
                index = np.argmin(abs(theseDists))
                freqs[theseDists[index] + 100] += 1
    return freqs

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

viterbi = io.load_tsv(
    directory + "/viterbi_dyad_calls",
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

redund = np.loadtxt(
    "../Data/Brogaard_RedundantDyads.tsv",
    dtype = str,
    delimiter = "\t"
)
redundant = {}
for chrom in np.unique(redund[:,0]):
    redundant[str(chrom)] = redund[redund[:,0] == chrom, 1].astype(int)

rotational_positions = np.loadtxt(
    "../Data/Nucleosome_rotational_categories.tsv",
    dtype = str,
    skiprows = 1,
    delimiter = "\t"
)
anti_rotational_positions = rotational_positions[
    rotational_positions[:,1] == "Minor_In",
    0
].astype(int)
rotational_positions = rotational_positions[
    rotational_positions[:,1] == "Minor_Out",
    0
].astype(int)

#%% Positions

notUnique = {}
for chrom in redundant:
    notUnique[chrom] = redundant[chrom][
        ~np.isin(redundant[chrom], brogaard[chrom])
    ]

redundRef = np.zeros((201, 6), dtype = int)
redundRef[:,0] = np.arange(-100,101)
redundRef[:,1] = getDist(greedy, redundant)
redundRef[:,2] = getDist(viterbi, redundant)
redundRef[:,3] = getDist(weiner, redundant)
redundRef[:,4] = getDist(brogaard, notUnique)
redundRef[:,5] = getDist(notUnique, brogaard)

directory = "NucleosomePattern/composite_models/composite_ice/distance_evals/"
cols = [
    "distance",
    "greedyVsRedund",
    "viterbiVsRedund",
    "weinerVsRedund",
    "brogaardVsRedund",
    "redundVsBrogaard"
]
np.savetxt(
    directory + "/redundant_freqTable.tsv",
    redundRef,
    delimiter = "\t",
    comments = "",
    header = "\t".join(cols),
    fmt = "%d"
)

#%%% Periodic Proportions

freqs = np.zeros(201, dtype = int)
for chrom in redundant:
    for i in range(len(redundant[chrom])):
        theseDists = redundant[chrom][i] - redundant[chrom]
        theseDists = theseDists[
            (theseDists >= -100) & (theseDists <= 100) & (theseDists != 0)
        ]
        if (len(theseDists) > 0):
            freqs[theseDists + 100] += 1

redundRef = np.column_stack((redundRef, freqs))
cols += ["redundVsRedund"]
redundRef[100,:] = 0

periodic_proportions = []
anti_periodic_proportions = []

for i in range(1, np.shape(redundRef)[1]):
    periodic_proportions.append(
        nuc.get_periodic_proportion(redundRef[:,[0,i]], rotational_positions)
    )
    anti_periodic_proportions.append(
        nuc.get_periodic_proportion(
            redundRef[:,[0,i]],
            anti_rotational_positions
        )
    )

with open(directory + "/redundant_periodic.tsv", "w") as file:
    for i in range(1, np.shape(redundRef)[1]):
        file.write(cols[i]+"\t")
        file.write("\t".join([str(s) for s in periodic_proportions[i-1]]))
        file.write("\n")

with open(directory + "/redundant_anti_periodic.tsv", "w") as file:
    for i in range(1, np.shape(redundRef)[1]):
        file.write(cols[i]+"\t")
        file.write("\t".join([str(s) for s in anti_periodic_proportions[i-1]]))
        file.write("\n")

#%% Offsets

notUnique = {}
for chrom in redundant:
    notUnique[chrom] = redundant[chrom][
        ~np.isin(redundant[chrom], brogaard[chrom])
    ]

redundRef = np.zeros((201, 6), dtype = int)
redundRef[:,0] = np.arange(-100,101)
redundRef[:,1] = getOffset(greedy, redundant)
redundRef[:,2] = getOffset(viterbi, redundant)
redundRef[:,3] = getOffset(weiner, redundant)
redundRef[:,4] = getOffset(brogaard, notUnique)
redundRef[:,5] = getOffset(notUnique, brogaard)

directory = "NucleosomePattern/composite_models/composite_ice/distance_evals/"
cols = [
    "distance",
    "greedyVsRedund",
    "viterbiVsRedund",
    "weinerVsRedund",
    "brogaardVsRedund",
    "redundVsBrogaard"
]
np.savetxt(
    directory + "/redundant_offsetFreqTable.tsv",
    redundRef,
    delimiter = "\t",
    comments = "",
    header = "\t".join(cols),
    fmt = "%d"
)

#%% Center to Centers

c2c = np.zeros((101,2), dtype = int)
c2c[:,0] = np.arange(101)

for chrom in redundant:
    theseDyads = np.sort(redundant[chrom])
    for i in range(len(theseDyads) - 1):
        for j in range(i + 1, len(theseDyads)):
            d = theseDyads[j] - theseDyads[i]
            if (d <= 100):
                c2c[d,1] += 1
            else:
                break

directory = "NucleosomePattern/composite_models/composite_ice/distance_evals/"
np.savetxt(
    directory + "/redundant_center_to_centers.tsv",
    c2c,
    delimiter = "\t",
    fmt = "%d"
)

#%% End

print("===== 0e0_redundant_nucleosomes.py - Exiting Properly =====")