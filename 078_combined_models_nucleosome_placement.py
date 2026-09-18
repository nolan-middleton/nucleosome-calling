#############################################################################
# This file will calculate the combined scores and place the nucleosomes in #
# the combined models.                                                  #
#############################################################################

print("===== 078_combined_models_nucleosome_placement.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

#%% Load Data

# Load data
print("> Loading data...")
chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

combo_dir = "NucleosomePattern/composite_models"
combos = {}
for combo in io.load_list(combo_dir + "/combinations.txt"):
    combos[combo] = io.load_json(combo_dir + "/" + combo + "/info.json")

models = []
for combo in combos:
    models += combos[combo]["models"]
models = np.unique(models)

scores = {}
for model in models:
    scores[model] = io.load_tsv(
        "NucleosomePattern/" + model + "/genome_scores",
        chroms
    )

#%% Combine Scores

print("> Combining scores...")
combined_scores = {}
for combo in combos:
    print(">> " + combo + "...")
    combined_scores[combo] = {}
    
    a = combos[combo]["optimal_weights"]
    
    for chrom in chroms:
        print(">>> " + chrom + "...")
        locs, indices1, indices2 = np.intersect1d(
            scores[combos[combo]["models"][0]][chrom][:,0],
            scores[combos[combo]["models"][1]][chrom][:,0],
            return_indices = True
        )
        values = (
            a*scores[combos[combo]["models"][0]][chrom][indices1,1]
        ) + (
            (1 - a)*scores[combos[combo]["models"][1]][chrom][indices1,1]
        )
        combined_scores[combo][chrom] = np.column_stack((locs, values))

#%% Outputs

print("> Outputting...")
for combo in combined_scores:
    print(">> " + combo + "...")
    thisDir = combo_dir + "/" + combo
    theseScores = combined_scores[combo]
    
    print(">>> Outputting scores...")
    io.output_tsv(theseScores, thisDir + "/genome_scores", chroms = True)
    io.output_wig(theseScores, thisDir + "/genome_scores.wig")
    
    print(">>> Placements (threshold = 1)...")
    placements = nuc.greedy_placement(
        theseScores,
        threshold = 1
    )
    nuc.output_placements(placements, thisDir + "/dyad_calls")
    
    print(">>> Placements (threshold = 1)...")
    placements = nuc.greedy_placement(
        theseScores,
        threshold = 1,
        exclusion_radius = 107
    )
    nuc.output_placements(placements, thisDir + "/dyad_calls_107bp")
    
    print(">>> Placements (no threshold)...")
    placements = nuc.greedy_placement(theseScores, threshold = -np.inf)
    nuc.output_placements(placements, thisDir + "/total_dyad_calls")
    
    # print(">>> Dynamic probability distribution...")
    # dynamic_dist = nuc.get_dynamic_probability_distribtion(theseScores)
    # io.output_tsv(
    #     dynamic_dist,
    #     thisDir + "/dynamic_probability_distribution",
    #     chroms = True
    # )
    # io.output_wig(dynamic_dist,thisDir+"/dynamic_probability_distribution.wig")
    
    print("> Placements (Viterbi)...")
    placements = nuc.viterbi_placement(theseScores)
    nuc.output_placements(placements, thisDir + "/viterbi_dyad_calls")
    
    print("> Placements (Validation)...")
    validation = nuc.viterbi_placement_with_overlap(theseScores, overlap = 0)
    
    print("> Placements (Viterbi with 40bp overlap)...")
    placements = nuc.viterbi_placement_with_overlap(theseScores)
    nuc.output_placements(placements,thisDir+"/viterbi_40b_overlap_dyad_calls")
    
    print("> Placements (Viterbi with 20bp overlap)...")
    placements = nuc.viterbi_placement_with_overlap(theseScores, overlap = 20)
    nuc.output_placements(placements,thisDir+"/viterbi_20b_overlap_dyad_calls")

#%% End

S="===== 078_combined_models_nucleosome_placement.py - Exiting Properly ====="
print(S)