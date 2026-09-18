#############################################################################
# This file will place the nucleosomes using the scores for the cl_ice_log2 #
# model.                                                                    #
#############################################################################

print("===== 072_cl_ice_log2_nucleosome_placement.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Variables
dataset = "cl_ice_log2"

#%% Load Data

print("> Loading data...")
chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

scores = io.load_tsv("NucleosomePattern/"+dataset+"/genome_scores", chroms)

#%% Greedy Placement (Threshold = 1)

print("> Placements (threshold = 1)...")
placements = nuc.greedy_placement(scores, threshold = 1)
print(">> Outputting...")
nuc.output_placements(placements, "NucleosomePattern/"+dataset+"/dyad_calls")

#%% Greedy Placement (No Threshold)

print("> Placements (no threshold)...")
total_placements = nuc.greedy_placement(scores, threshold = -np.inf)
print(">> Outputting...")
nuc.output_placements(
    total_placements,
    "NucleosomePattern/" + dataset + "/total_dyad_calls"
)

#%% Dynamic Probability Distribution

print("> Dynamic probability distribution...")
dynamic_dist = nuc.get_dynamic_probability_distribtion(scores)
print(">> Outputting...")
io.output_tsv(
    dynamic_dist,
    "NucleosomePattern/"+dataset+"/dynamic_probability_distribution",
    chroms = True
)
io.output_wig(
    dynamic_dist,
    "NucleosomePattern/"+dataset+"/dynamic_probability_distribution.wig"
)

#%% Viterbi Placement

print("> Placements (Viterbi)...")
viterbi_placements = nuc.viterbi_placement(scores)
print(">> Outputting...")
nuc.output_placements(
    viterbi_placements,
    "NucleosomePattern/" + dataset + "/viterbi_dyad_calls"
)

#%% End

print("===== 072_cl_ice_log2_nucleosome_placement.py - Exiting Properly =====")