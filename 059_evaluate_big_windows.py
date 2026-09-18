#############################################################################
# This file will evaluate the scores given to all of the models across a    #
# 160bp window, finding peaks and deciding if the algorithm correctly       #
# placed the nucleosome.                                                    #
#############################################################################

print("===== 059_evaluate_big_windows.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Defs
def find_overlap(offset_array):
    overlaps = np.zeros(np.shape(offset_array))
    
    overlaps[offset_array < 0] = 147 + offset_array[offset_array < 0]
    overlaps[offset_array >= 0] = 147 - offset_array[offset_array >= 0]
    
    return np.clip(overlaps, 0, 147)

def find_rotational_offset(offset_array):
    diff = np.abs(offset_array).astype(float)
    while (np.max(diff) > 10.3):
        diff[diff > 10.3] -= 10.3
    
    alt_diff = diff - 10.3
    diff[diff > np.abs(alt_diff)] = alt_diff[diff > np.abs(alt_diff)]
    
    return np.abs(diff)

# Variables
datasets = [
    "unfloored_cl",
    "cl_rt_log2",
    "cl_ice_log2",
    "cl_rt_diff",
    "cl_ice_diff",
    "cl_rt_abs_diff",
    "cl_ice_abs_diff"
]

gen = np.random.default_rng(872689294382963979343294329423)

# The bootstrapping samples
N = 50000

#%% Load Data

test_points = np.loadtxt(
    "evaluationPositions.tsv",
    delimiter = "\t",
    dtype = str,
    skiprows = 1
)
M = np.shape(test_points)[0]

nucleosome_rows = test_points[:,0] == "True"

W = 80

#%% Bootstrapping

print("> Bootstrapping...")
bootstrapped_nucleosome_overlaps = np.zeros(N)
bootstrapped_nucleosome_rotations = np.zeros(N)
bootstrapped_random_overlaps = np.zeros(N)
bootstrapped_random_rotations = np.zeros(N)
for i in range(N):
    if (i % 10000 == 0):
        print(">>> " + str(i) + "/" + str(N) + "...")
    thesePositions = gen.integers(-W, W, M, endpoint = True)
    
    overlaps = find_overlap(thesePositions)
    bootstrapped_nucleosome_overlaps[i] = np.mean(overlaps[nucleosome_rows])
    bootstrapped_random_overlaps[i] = np.mean(overlaps[~nucleosome_rows])
    
    rots = find_rotational_offset(thesePositions)
    bootstrapped_nucleosome_rotations[i] = np.mean(rots[nucleosome_rows])
    bootstrapped_random_rotations[i] = np.mean(rots[~nucleosome_rows])

np.savetxt(
    "NucleosomePattern/bigWindowNucleosomeBootstrap.tsv",
    np.column_stack(
        (bootstrapped_nucleosome_overlaps, bootstrapped_nucleosome_rotations)
    ),
    delimiter = "\t"
)
np.savetxt(
    "NucleosomePattern/bigWindowRandomBootstrap.tsv",
    np.column_stack(
        (bootstrapped_random_overlaps, bootstrapped_random_rotations)
    ),
    delimiter = "\t"
)

#%% Main Loop

print("> Individual models...")
for dataset in datasets:
    print(">> " + dataset + "...")
    
    evaluation_data = np.loadtxt(
        "NucleosomePattern/" + dataset + "/models/bigWindowEvalBayes.tsv",
        delimiter = "\t",
        skiprows = 1
    )
    
    nucleosome_positions = np.argmax(evaluation_data, axis = 1) - W
    overlaps = find_overlap(nucleosome_positions)
    rotational_positions = find_rotational_offset(nucleosome_positions)
    
    mean_nucleosome_overlap = np.mean(overlaps[nucleosome_rows])
    mean_nucleosome_rot = np.mean(rotational_positions[nucleosome_rows])
    
    nucleosome_overlap_phat = np.sum(
        bootstrapped_nucleosome_overlaps >= mean_nucleosome_overlap
    ) / N
    nucleosome_rotation_phat = np.sum(
        np.abs(bootstrapped_nucleosome_rotations) <= abs(mean_nucleosome_rot)
    ) / N
    
    mean_random_overlap = np.mean(overlaps[~nucleosome_rows])
    mean_random_rotation = np.mean(rotational_positions[~nucleosome_rows])
    
    random_overlap_phat = np.sum(
        bootstrapped_random_overlaps >= mean_random_overlap
    ) / N
    random_rotation_phat = np.sum(
        np.abs(bootstrapped_random_rotations) <= abs(mean_random_rotation)
    ) / N
    
    np.savetxt(
        "NucleosomePattern/" + dataset + "/models/bigWindowEvalResults.tsv",
        np.column_stack((nucleosome_positions,overlaps,rotational_positions)),
        delimiter = "\t",
        header = "offset\toverlap\trotational_offset",
        comments = ""
    )
    
    io.output_json(
        "NucleosomePattern/" + dataset + "/models/bigWindowEvalStats.json",
        {
            "nucleosome": {
                "mean_overlap": mean_nucleosome_overlap,
                "mean_rotation": mean_nucleosome_rot,
                "overlap_phat": nucleosome_overlap_phat,
                "rotation_phat": nucleosome_rotation_phat
            },
            "random": {
                "mean_overlap": mean_random_overlap,
                "mean_rotation": mean_random_rotation,
                "overlap_phat": random_overlap_phat,
                "rotation_phat": random_rotation_phat
            }
        }
    )

#%% End

print("===== 059_evaluate_big_windows.py - Exiting Properly =====")