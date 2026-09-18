#############################################################################
# This file will optimize the weight factors for combining different scores #
# for different models.                                                     #
#############################################################################

print("===== 077_combined_score_optimization.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os
import scipy.optimize as opt

# Defs
def objective(a, scores1, scores2):
    combined_scores = a*scores1 + (1 - a)*scores2
    locs = np.argmax(combined_scores, axis = 1)
    
    return -np.sum(locs == R)
    

# Variables
combos = {
    "composite_ice": ["cl_ice_log2", "cl_ice_abs_diff"]
}

N = 1000 # The number of nucleosomes to randomly sample
R = 80 # The number of bases to expand on either side of the true dyad
granularity = 1000 # The number of points for the graph

rng = np.random.default_rng(12789361789031374128903712890731289730)

#%% Load Data

print("> Loading data..")
chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

models = []
for combo in combos:
    models += combos[combo]
models = np.unique(models)

scores = {}
for model in models:
    scores[model] = io.load_tsv(
        "NucleosomePattern/" + model + "/genome_scores",
        chroms
    )

nucleosomes = io.load_tsv("PreprocessedData/nucleosomes", chroms, dtype = int)

#%% Select Nucleosomes

print("> Selecting " + str(N) + " random nucleosome positions...")
selected_nucs = {}
chromIndices = rng.integers(0, len(chroms), N, dtype = int)
for index in np.unique(chromIndices):
    thisChrom = chroms[index]
    selected_nucs[thisChrom] = nucleosomes[thisChrom][
        rng.integers(
            0,
            np.shape(nucleosomes[thisChrom])[0],
            np.sum(chromIndices == index)
        ),
        2
    ]

selected_positions = {}
for chrom in selected_nucs:
    selected_positions[chrom] = np.repeat(
        selected_nucs[chrom][:,np.newaxis],
        2*R + 1,
        axis = 1
    ) + np.repeat(
        np.arange(-R, R + 1, dtype = int)[np.newaxis, :],
        len(selected_nucs[chrom]),
        axis = 0
    )

#%% Get Window Scores

print("> Getting window scores...")
window_scores = {}
for model in scores:
    print(">> " + model + "...")
    window_scores[model] = {}
    for chrom in selected_positions:
        print("=> " + chrom + "...")
        M = np.shape(selected_positions[chrom])[0]
        thisArray = np.zeros((M, 2*R+1))
        
        for i in range(M):
            theseIndices = np.isin(
                scores[model][chrom][:,0],
                selected_positions[chrom][i,:]
            )
            thisArray[i,:] = scores[model][chrom][theseIndices,1]
        window_scores[model][chrom] = thisArray

scores1 = {}
scores2 = {}
for combo in combos:
    scores1[combo] = np.vstack(
        tuple(
            [
                window_scores[combos[combo][0]][chrom]
                for chrom in window_scores[combos[combo][0]]
            ]
        )
    )
    scores2[combo] = np.vstack(
        tuple(
            [
                window_scores[combos[combo][1]][chrom]
                for chrom in window_scores[combos[combo][1]]
            ]
        )
    )

#%% Optimization

# Minimizing the parameter
print("> Minimizing")
minimizers = {}
for combo in combos:
    print(">> " + combo + "...")
    minimization = opt.minimize_scalar(
        objective,
        bounds = (0,1),
        args = (scores1[combo], scores2[combo])
    )
    
    if (not minimization.success):
        print("/!\\ Minimization did not succeed!")
    
    minimizers[combo] = minimization.x

#%% Graph Objective Function

print("> Objective function graph...")
x = np.linspace(0, 1, granularity)
graph = {}
for combo in combos:
    print(">> " + combo + "...")
    y = np.array([-objective(a, scores1[combo], scores2[combo]) for a in x])
    graph[combo] = np.column_stack((x, y))

#%% Outputs

print("> Outputting...")
out_dir = "NucleosomePattern/composite_models"
if (not os.path.isdir(out_dir)):
    os.mkdir(out_dir)

io.output_list(list(combos.keys()), out_dir + "/combinations.txt")
io.output_tsv(selected_nucs, out_dir + "/selected_nucleosomes", fmt = "%d")

for combo in combos:
    thisDir = out_dir + "/" + combo
    if (not os.path.isdir(thisDir)):
        os.mkdir(thisDir)
    
    io.output_json(
        thisDir + "/info.json",
        {
            "models": combos[combo],
            "optimal_weights": minimizers[combo]
        }
    )
    np.savetxt(thisDir + "/graph.tsv", graph[combo], delimiter = "\t")

#%% End

print("===== 077_combined_score_optimization.py - Exiting Properly =====")