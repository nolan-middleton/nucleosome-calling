#############################################################################
# This file will generate a set of random positions across the genome as    #
# "non-nucleosome" positions for use as a baseline model.                   #
#############################################################################

print("===== 020_generate_random_positions.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import os as os

# Initialize random seed
gen = np.random.default_rng(8923740982379047329041793)

# Define how many random datasets you want
N = 10

# Width of nucleosome window
w = 146

# Make initial directory
if (not os.path.isdir("RandomPattern")):
    os.mkdir("RandomPattern")
if (not os.path.isdir("RandomPattern/data")):
    os.mkdir("RandomPattern/data")

with open("RandomPattern/datasets.txt", "w") as file:
    file.write(str(N))

#%% Load Data

# Read in genome and blocks
print("> Loading data...")
chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM") # We don't want chrM
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)

badBlocks = {}
for chrom in blocks:
    badBlocks[chrom] = blocks[chrom][blocks[chrom][:,2] == 1,0:2]

# Get nucleosome positions - we want to know how many points to grab
nucChroms = io.load_list("PreprocessedData/nucleosomes/chromosomes.txt")
nucleosomes = io.load_tsv("PreprocessedData/nucleosomes", nucChroms)

#%% Main Loop

# Find the number of positions to generate and width of nucleosome window
m = 0
for chrom in nucleosomes:
    m += np.shape(nucleosomes[chrom])[0]

for i in range(N):
    print("> " + str(i + 1) + "...")
    
    # Find random nucleosome positions
    random_positions = {}
    random_chroms = np.array(chroms)[gen.integers(0,len(chroms),m)]
    for chrom in np.unique(random_chroms):
        random_positions[str(chrom)] = np.zeros((0,3))
    
    for j in range(m):
        thisChrom = str(random_chroms[j])
        isBad = True
        while (isBad):
            thisStart = gen.integers(0,len(genome[thisChrom]) - w)
            isBad = False
            
            if (np.shape(badBlocks[thisChrom])[0] > 0):
                starts = badBlocks[thisChrom][:,0]
                stops = badBlocks[thisChrom][:,1]
                overlap = np.clip(
                    np.min(
                        np.column_stack(
                            (
                                thisStart + w - starts + 1,
                                stops - thisStart + 1
                            )
                        ),
                        axis = 1
                    ),
                    0,
                    w + 1
                )
                
                if (np.sum(overlap) > 0):
                    isBad = True
            if (np.shape(random_positions[thisChrom])[0] > 0):
                starts = random_positions[thisChrom][:,0]
                stops = random_positions[thisChrom][:,1]
                overlap = np.clip(
                    np.min(
                        np.column_stack(
                            (
                                thisStart + w - starts + 1,
                                stops - thisStart + 1
                            )
                        ),
                        axis = 1
                    ),
                    0,
                    w + 1
                )
                
                if (np.sum(overlap) > 0):
                    isBad = True
        
        random_positions[thisChrom] = np.vstack(
            (
                random_positions[thisChrom],
                [thisStart, thisStart + w, thisStart + w/2]
            )
        )
    
    io.output_tsv(random_positions,"RandomPattern/data/"+str(i),chroms=True)

#%% End

print("===== 020_generate_random_positions.py - Exiting Properly =====")