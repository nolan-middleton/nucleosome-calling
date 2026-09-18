#############################################################################
# This file will score the ENTIRE GENOME for the unfloored_cl model.        #
#############################################################################

print("===== 060_unfloored_cl_score_genome.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Variables
dataset = "unfloored_cl"
nucleosome_std_relation = lambda MLEs, means : MLEs["k"]*means + MLEs["l"]
model = "linear"
nucleosome_model = nuc.nucleosome_linear_model

random_std_relation = lambda MLEs, means : MLEs["k"]*means + MLEs["l"]

#%% Load Data

print("> Loading data...")
genome, blocks, data_minus, data_plus = io.load_dataset(
    "PreprocessedData",
    [dataset]
)
chroms = []
for chrom in genome:
    chroms.append(chrom)
chroms.remove("chrM") # No nucleosome in mitochondrial genome

random_MLE = io.load_json("RandomPattern/" + dataset + "/0/MLE_values.json")
random_out_dir = "RandomPattern/" + dataset + "/genome_ints"
if (not os.path.isdir(random_out_dir)):
    os.mkdir(random_out_dir)

pattern = np.loadtxt(
    "NucleosomePattern/" + dataset + "/damagePattern.tsv",
    delimiter = "\t"
)
offsets = pattern[0,:]

nucleosome_MLE = io.load_json(
    "NucleosomePattern/" + dataset + "/models/" + model + "/MLE_values.json"
)

#%% Score Genome

print("> Scoring genome...")
bayes, nucInts, randInts = nuc.score_genome(
    genome,
    blocks,
    chroms,
    data_minus,
    data_plus,
    offsets,
    nucleosome_model,
    nucleosome_std_relation,
    nucleosome_MLE,
    random_std_relation,
    random_MLE
)

io.output_tsv(bayes, "NucleosomePattern/" + dataset + "/genome_scores")
io.output_tsv(nucInts, "NucleosomePattern/"+dataset+"/genome_nucleosome_ints")
io.output_tsv(randInts, "NucleosomePattern/"+dataset+"/genome_random_ints")

io.output_wig(bayes, "NucleosomePattern/" + dataset + "/scores.wig")

#%% End

print("===== 060_unfloored_cl_score_genome.py - Exiting Properly =====")