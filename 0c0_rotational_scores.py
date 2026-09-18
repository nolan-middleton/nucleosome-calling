#############################################################################
# This file will score the nucleosomes rotational positioning by using the  #
# nucleosome score.                                                         #
#############################################################################

print("===== 0c0_rotational_scores.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

#%% Load Data

print("> Loading data...")
dirs = [
    "NucleosomePattern/cl_rt_log2",
    "NucleosomePattern/cl_ice_log2",
    "NucleosomePattern/cl_ice_abs_diff",
    "NucleosomePattern/composite_models/composite_ice"
]
num_dirs = len(dirs)

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

# Load in test nucleosome placements and scores
print(">> Putative dyad calls...")
dyads = [
    io.load_tsv(
        D + "/viterbi_dyad_calls",
        chromosomes = io.load_list(D + "/viterbi_dyad_calls/chromosomes.txt"),
        dtype = int
    )
    for D in dirs
]

print(">> Scores...")
scores = [io.load_tsv(D+"/genome_scores", chromosomes = chroms) for D in dirs]

# Load in genome data
print(">> Genome and blocks...")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)
for chrom in blocks:
    blocks[chrom] = blocks[chrom].astype(int)

print(">> Rotational positions...")
rotational_positions = np.loadtxt(
    "../Data/Nucleosome_rotational_categories.tsv",
    dtype = str,
    skiprows = 1,
    delimiter = "\t"
)
anti_rotational_positions = rotational_positions[
    rotational_positions[:,1] == "Minor_In",
    0
].astype(int)
rotational_positions = rotational_positions[
    rotational_positions[:,1] == "Minor_Out",
    0
].astype(int)

out_indices = rotational_positions + 73
in_indices = anti_rotational_positions + 73

#%% Main

RS = []
for i in range(len(dirs)):
    print("> " + dirs[i] + "...")
    
    RS.append({})
    
    for chrom in dyads[i]:
        print(">> " + chrom + "...")
        RS[-1][chrom] = np.zeros((len(dyads[i][chrom]), 2))
        RS[-1][chrom][:,0] = dyads[i][chrom]
        
        for j in range(np.shape(RS[-1][chrom])[0]):
            thisStart = RS[-1][chrom][j,0] - 73
            thisStop = RS[-1][chrom][j,0] + 73
            
            theseScores = scores[i][chrom][
                (scores[i][chrom][:,0] >= thisStart) & \
                    (scores[i][chrom][:,0] <= thisStop),
                1
            ]
            
            if (len(theseScores) != 147):
                RS[-1][chrom][j,1] = np.nan
            else:
                RS[-1][chrom][j,1] = np.sum(theseScores[out_indices]) - \
                    np.sum(theseScores[in_indices])
        RS[-1][chrom] = RS[-1][chrom][~np.isnan(RS[-1][chrom][:,1]), :]
    
    io.output_tsv(RS[-1], dirs[i] + "/rotational_scores")
    io.output_wig(RS[-1], dirs[i] + "/rotational_scores.wig")
    
    master_table = np.vstack(
        tuple(
            [
                np.column_stack(
                    (
                        np.repeat(chrom, np.shape(RS[-1][chrom])[0]),
                        RS[-1][chrom][:,0].astype(int),
                        RS[-1][chrom][:,1]
                    )
                )
                for chrom in RS[-1]
            ]
        )
    )
    
    np.savetxt(
        dirs[i] + "/rotational_scores.tsv",
        master_table,
        delimiter = "\t",
        fmt = "%s"
    )

#%% End

print("===== 0c0_rotational_scores.py - Exiting Properly =====")