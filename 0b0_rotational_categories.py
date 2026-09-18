#############################################################################
# This file will identify "blocks" of nucleosomes with similar rotational   #
# phasing.                                                                  #
#############################################################################

print("===== 0b0_rotational_categories.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Generator

gen = np.random.default_rng(786743472378946232374238943890478239074238974)

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

refDirs = [
    "NucleosomePattern/ReferenceRotationalCategories/all",
    "NucleosomePattern/ReferenceRotationalCategories/weiner"
]

if (not os.path.isdir("NucleosomePattern/ReferenceRotationalCategories")):
    os.mkdir("NucleosomePattern/ReferenceRotationalCategories")

for path in refDirs:
    if (not os.path.isdir(path)):
        os.mkdir(path)

dirs += refDirs

normal = len(dirs)

#%% Random Dyads

print("> Getting random dyads...")

if (not os.path.isdir("RandomPattern/RotationalCategories")):
    os.mkdir("RandomPattern/RotationalCategories")

random_dyads = []
for i in range(10000):
    print(">> " + str(i) + "...")
    if (not os.path.isdir("RandomPattern/RotationalCategories/" + str(i))):
        os.mkdir("RandomPattern/RotationalCategories/" + str(i))
    
    random_dyads.append({})
    for chrom in chroms:
        n = len(dyads[2][chrom])
        available_positions = np.arange(len(genome[chrom]))
        
        bad = blocks[chrom][blocks[chrom][:,2] == 1,:]
        
        for j in range(np.shape(bad)[0]):
            available_positions = available_positions[
                (available_positions<bad[j,0]) | (available_positions>bad[j,1])
            ]
        
        random_dyads[-1][chrom] = np.zeros(n, int) - 1
        m = 0
        isDone = False
        while (not isDone):
            theseLocs = available_positions[
                gen.integers(0, len(available_positions), n - m)
            ]
            
            j = 0
            while (j < len(theseLocs)):
                theseLocs = np.concatenate(
                    (
                        theseLocs[:j+1],
                        theseLocs[j+1:][abs(theseLocs[j+1:]-theseLocs[j])>146]
                    )
                )
                j += 1
            
            random_dyads[-1][chrom][m:m+len(theseLocs)] = theseLocs
            m += len(theseLocs)
            
            badLocs = np.concatenate(
                tuple([np.arange(L - 146, L + 146 + 1) for L in theseLocs])
            )
            
            available_positions = available_positions[
                ~np.isin(available_positions, badLocs)
            ]
            
            if (m >= n) or (len(available_positions) <= 0):
                isDone = True
                break
        
        random_dyads[-1][chrom] = random_dyads[-1][chrom][
            random_dyads[-1][chrom] >= 0
        ]
    
    io.output_tsv(
        random_dyads[-1],
        "RandomPattern/RotationalCategories/" + str(i) + "/dyads",
        fmt = "%d",
        chroms = True
    )

dyads += random_dyads

dirs += ["RandomPattern/RotationalCategories/" + str(i) for i in range(10000)]

#%% Main

for i in range(len(dirs)):
    print("> " + dirs[i] + "...")
    theseCats = {}
    for chrom in dyads[i]:
        theseCats[chrom] = []
        thisList = []
        for dyad in dyads[i][chrom].tolist():
            if (thisList == []):
                thisList.append(dyad)
            else:
                d = dyad - thisList[-1]
                offset = d % 10.1
                if (min(offset, 10.1 - offset) <= 1.5):
                    thisList.append(dyad)
                else:
                    theseCats[chrom].append(thisList)
                    thisList = [dyad]
        theseCats[chrom].append(thisList)
    
    print(">> JSON...")
    io.output_json(dirs[i] + "/rotational_categories.json", theseCats)
    
    if (i <= normal):
        maxLen = 0
        for chrom in theseCats:
            maxLen = max(maxLen, len(theseCats[chrom]))
        
        wigCats = [{} for n in range(maxLen)]
        for chrom in theseCats:
            for n in range(len(theseCats[chrom])):
                wigCats[n][chrom] = theseCats[chrom][n]
    
        print("> Wig track...")
        wigCat = {}
        for n in range(len(wigCats)):
            for chrom in wigCats[n]:
                if (chrom not in wigCat):
                    wigCat[chrom] = np.zeros((0,2), int)
                
                wigCat[chrom] = np.vstack(
                    (
                        wigCat[chrom],
                        np.column_stack(
                            (
                                wigCats[n][chrom],
                                np.zeros(len(wigCats[n][chrom]),int) + (n % 2) + 1
                            )
                        )
                    )
                )
        
        io.output_wig(wigCat, dirs[i] + "/rotational_categories.wig")     
        
        print(">> Two wigs...")
        wigCat1 = np.zeros((0,4), dtype = "<U32")
        wigCat2 = np.zeros((0,4), dtype = "<U32")
        for chrom in wigCat:
            one = wigCat[chrom][:,1] == 1
            wigCat1 = np.vstack(
                (
                    wigCat1,
                    np.column_stack(
                        (
                            np.repeat(chrom, np.sum(one)),
                            wigCat[chrom][one,0] - 73,
                            wigCat[chrom][one,0] + 73,
                            wigCat[chrom][one,0]
                        )
                    )
                )
            )
            wigCat2 = np.vstack(
                (
                    wigCat2,
                    np.column_stack(
                        (
                            np.repeat(chrom, np.sum(~one)),
                            wigCat[chrom][~one,0] - 73,
                            wigCat[chrom][~one,0] + 73,
                            wigCat[chrom][~one,0]
                        )
                    )
                )
            )
            
            io.output_wig_interval_from_table(
                wigCat1,
                dirs[i] + "/rotational_categories1.wig",
                highlight = True
            )
            io.output_wig_interval_from_table(
                wigCat2,
                dirs[i] + "/rotational_categories2.wig",
                highlight = True
            )

#%% End

print("===== 0b0_rotational_categories.py - Exiting Properly =====")