#############################################################################
# This file will classify the "blocks" of rotationally-consistent           #
# nucleosomes.                                                              #
#############################################################################

print("===== 0b4_rotational_categories_classify.py - Starting =====")

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
num_dirs = len(dirs)

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

dirs += [
    "NucleosomePattern/ReferenceRotationalCategories/all",
    "NucleosomePattern/ReferenceRotationalCategories/weiner"
]

dirs += ["RandomPattern/RotationalCategories/"+str(i) for i in range(10000)]

print(">> TSS and PAS")
TSS = io.load_tsv("PreprocessedData/TSS", chroms, dtype = str)
PAS = io.load_tsv("PreprocessedData/PAS", chroms, dtype = str)

genes = {}
for chrom in chroms:
    genes[chrom] = TSS[chrom][:,1]
    
    theseTSS = TSS[chrom][:,0].astype(int)
    thesePAS = -np.ones(len(theseTSS), dtype = int)
    for i in range(len(genes[chrom])):
        thisPAS = PAS[chrom][PAS[chrom][:,1] == genes[chrom][i],0]
        if (len(thisPAS) > 0):
            thesePAS[i] = thisPAS[0]
    
    bad = thesePAS == -1
    
    theseCoords = np.column_stack((theseTSS, thesePAS))
    
    genes[chrom] = np.column_stack(
        (genes[chrom], np.min(theseCoords,axis=1), np.max(theseCoords,axis=1))
    )[~bad]

#%% Main Loop

for D in dirs:
    print("> " + D + "...")
    data = io.load_json(D + "/rotational_categories.json")
    
    classes = {}
    for chrom in data:
        thisList = []
        for entry in data[chrom]:
            found = False
            for item in entry:
                search = (item > genes[chrom][:,1].astype(int)) & \
                    (item < genes[chrom][:,2].astype(int))
                if (True in search):
                    found = True
                    break
            thisList.append([len(entry), found])
        classes[chrom] = thisList
    
    io.output_json(D + "/rotational_classes.json", classes)
    
    data = io.load_json(D + "/rotational_categories_lenient.json")
    
    classes = {}
    for chrom in data:
        thisList = []
        for entry in data[chrom]:
            found = False
            for item in entry:
                search = (item > genes[chrom][:,1].astype(int)) & \
                    (item < genes[chrom][:,2].astype(int))
                if (True in search):
                    found = True
                    break
            thisList.append([len(entry), found])
        classes[chrom] = thisList
    
    io.output_json(D + "/rotational_classes_lenient.json", classes)

#%% End

print("===== 0b4_rotational_categories_classify.py - Exiting Properly =====")