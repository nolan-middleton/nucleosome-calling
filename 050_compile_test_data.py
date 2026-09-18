#############################################################################
# This file will compile a dataset of 1000 data points: 500 nucleosomes and #
# 500 random to be used in model evaluation.                                #
#############################################################################

print("===== 050_compile_test_data.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc

# Initialize random seed
gen = np.random.default_rng(529463839389817297342823637354273427934)

# How many of each type (nucleosome or random) position do you want?
M = 500

#%% Load Data

# Read in genome and blocks
print("> Loading data...")
chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM") # We don't want chrM
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)

nucChroms = io.load_list("PreprocessedData/nucleosomes/chromosomes.txt")
nucleosomes = io.load_tsv("PreprocessedData/nucleosomes", nucChroms)

# Read in random datasets
N = nuc.get_N()

random_data = {}
for i in range(N):
    chroms = io.load_list("RandomPattern/data/" + str(i) + "/chromosomes.txt")
    for chrom in chroms:
        if (chrom not in random_data):
            random_data[chrom] = np.zeros((0,3))
    this_data = io.load_tsv("RandomPattern/data/" + str(i), chroms)
    
    for chrom in chroms:
        random_data[chrom]=np.vstack((random_data[chrom],this_data[chrom]))

randChroms = list(random_data.keys())

#%% Compile Random Datasets

print("> Compiling random dataset...")
output_table = np.array([["answer", "chrom", "start", "stop", "dyad"]])

print(">> Nucleosome positions...")
chroms = np.zeros(M, dtype = "<U32")
indices = np.zeros(M, dtype = int)
starts = np.zeros(M, dtype = int)
stops = np.zeros(M, dtype = int)
dyads = np.zeros(M, dtype = int)
for m in range(M):
    chroms[m] = nucChroms[gen.integers(0, len(nucChroms))]
    indices[m] = gen.integers(0, np.shape(nucleosomes[chroms[m]])[0])
    starts[m] = nucleosomes[chroms[m]][indices[m],0]
    stops[m] = nucleosomes[chroms[m]][indices[m],1]
    dyads[m] = nucleosomes[chroms[m]][indices[m],2]
    
output_table = np.vstack(
    (
        output_table,
        np.column_stack((np.repeat(True, M), chroms, starts, stops, dyads))
    )
)

print(">> Random positions...")
chroms = np.zeros(M, dtype = "<U32")
indices = np.zeros(M, dtype = int)
starts = np.zeros(M, dtype = int)
stops = np.zeros(M, dtype = int)
dyads = np.zeros(M, dtype = int)
for m in range(M):
    chroms[m] = randChroms[gen.integers(0, len(nucChroms))]
    indices[m] = gen.integers(0, np.shape(random_data[chroms[m]])[0])
    starts[m] = random_data[chroms[m]][indices[m],0]
    stops[m] = random_data[chroms[m]][indices[m],1]
    dyads[m] = random_data[chroms[m]][indices[m],2]
    
output_table = np.vstack(
    (
        output_table,
        np.column_stack((np.repeat(False, M), chroms, starts, stops, dyads))
    )
)

#%% Outputs

print("> Outputting...")
np.savetxt("evaluationPositions.tsv", output_table, fmt="%s", delimiter="\t")

#%% End

print("===== 050_compile_test_data.py - Exiting Properly =====")