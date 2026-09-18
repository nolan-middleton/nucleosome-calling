#############################################################################
# This file will optimize the weight factors for combining different scores #
# for different models.                                                     #
#############################################################################

print("===== 102_combined_shuffled_scores.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

#%% Load Data

print("> Loading data..")
chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

scores = [
    io.load_tsv(
        "NucleosomePattern/" + model + "/shuffled_genome_scores",
        chroms
    )
    for model in ["cl_ice_log2", "cl_ice_abs_diff"]
]

w = io.load_json(
    "NucleosomePattern/composite_models/composite_ice/info.json"
)["optimal_weights"]

#%% Main

print("> Combining...")
combined = {}
for chrom in scores[0]:
    combined[chrom] = np.column_stack(
        (
            scores[0][chrom][:,0],
            w*scores[0][chrom][:,1] + (1-w)*scores[1][chrom][:,1]
        )
    )

io.output_tsv(
    combined,
    "NucleosomePattern/composite_models/composite_ice/shuffled_scores"
)

#%% End

print("===== 077_combined_score_optimization.py - Exiting Properly =====")