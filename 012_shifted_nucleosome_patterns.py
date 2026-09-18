#############################################################################
# This file will find the nucleosome damage pattern with a shift of 10.1bp  #
#############################################################################

print("===== 012_shifted_nucleosome_patterns.py - Starting =====")

#%% Setup

# Imports
import os as os
import DataIO as io
import Nucleosomes as nuc

#%% Load Data

# Get data
print("> Loading data...")
dataset_names = [
    "unfloored_cl",
    "cl_rt_log2",
    "cl_ice_log2",
    "cl_rt_diff",
    "cl_ice_diff",
    "cl_rt_abs_diff",
    "cl_ice_abs_diff"
]
data = io.load_dataset(
    "PreprocessedData",
    dataset_names
)
genome, blocks = data[0:2]
minus_data = [data[2 + 2*i] for i in range(len(dataset_names))]
plus_data = [data[3 + 2*i] for i in range(len(dataset_names))]

nucChroms = io.load_list("PreprocessedData/all_nucleosomes/chromosomes.txt")
nucleosomes = io.load_tsv("PreprocessedData/all_nucleosomes", nucChroms)

# Make output directory
if (not os.path.isdir("NucleosomePattern_shifted")):
    os.mkdir("NucleosomePattern_shifted")

#%% Main Loop

for chrom in nucleosomes:
    nucleosomes[chrom] += 10

this_dir = "NucleosomePattern_shifted/"
if (not os.path.isdir(this_dir)):
    os.mkdir(this_dir)

results = []
for I in range(len(dataset_names)):
    dataset_name = dataset_names[I]
    print(">> " + dataset_name + "...")
    if (not os.path.isdir(this_dir + "/" + dataset_name)):
        os.mkdir(this_dir + "/" + dataset_name)
    
    results.append(
        nuc.make_nucleosome_pattern(
            nucleosomes,
            minus_data[I],
            plus_data[I],
            this_dir + "/" + dataset_name
        )
    )

#%% End

print("===== 012_shifted_nucleosome_patterns.py - Exiting Properly =====")