#############################################################################
# This file will find the damage pattern of the random sites.               #
#############################################################################

print("===== 021_find_random_pattern.py - Starting =====")

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

# Read in how many position sets we'll handle
N = nuc.get_N()

#%% Main Loop

# Loop through every dataset
for I in range(len(dataset_names)):
    dataset_name = dataset_names[I]
    print("> " + dataset_name + "...")
    if (not os.path.isdir("RandomPattern/" + dataset_name)):
        os.mkdir("RandomPattern/" + dataset_name)
    # Loop through every random position set and compile the random models
    for i in range(N):
        print(">> " + str(i + 1) + "/" + str(N) + "...")
        data_dir = "RandomPattern/data/" + str(i)
        
        randChroms = io.load_list(data_dir + "/chromosomes.txt")
        random_spans = io.load_tsv(data_dir, randChroms)
        
        # Make output directory
        out_dir = "RandomPattern/" + dataset_name + "/" + str(i)
        if (not os.path.isdir(out_dir)):
            os.mkdir(out_dir)
        
        results = nuc.make_nucleosome_pattern(
            random_spans,
            minus_data[I],
            plus_data[I],
            out_dir
        )

#%% End

print("===== 021_find_random_pattern.py - Exiting Properly =====")