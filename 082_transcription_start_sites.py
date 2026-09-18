#############################################################################
# This file will evaluate the nucleosome placement based on the scores and  #
# nucleosome placements around the transcription start sites.               #
#############################################################################

print("===== 082_transcription_start_sites.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Defs
def get_all_TSS_dyads(R, dyads, TSS_sites, TSS_strands):
    D = an.get_distance_matrix(dyads, TSS_sites, signed = True)
    
    output = {}
    for chrom in D:
        output[chrom] = np.zeros(
            (len(TSS_sites[chrom]), R[1] - R[0] + 1),
            dtype = bool
        )
        for j in range(np.shape(D[chrom])[1]):
            if (TSS_strands[chrom][j] == "-"):
                D[chrom][:,j] *= -1
            theseDists = np.sort(
                D[chrom][:,j][
                    np.logical_and(D[chrom][:,j]>=R[0], D[chrom][:,j]<=R[1])
                ]
            )
            output[chrom][j,theseDists] = True
    
    return output

def get_TSS_frequencies(R, dyads, TSS_sites, TSS_strands):
    frequencies = np.arange(R[0],R[1]+1)
    frequencies = np.column_stack(
        (frequencies, np.zeros(len(frequencies), dtype = int))
    )
    
    D = an.get_distance_matrix(dyads, TSS_sites, signed = True)
    
    for chrom in D:
        for j in range(np.shape(D[chrom])[1]):
            if (TSS_strands[chrom][j] == "-"):
                D[chrom][:,j] *= -1
            theseDists = np.sort(
                D[chrom][:,j][
                    np.logical_and(D[chrom][:,j]>=R[0], D[chrom][:,j]<=R[1])
                ]
            )
            frequencies[np.isin(frequencies[:,0], theseDists),1] += 1
    return frequencies

def get_all_TSS_scores(R, scores, TSS_sites, TSS_strands):
    output = {}
    for chrom in TSS_sites:
        output[chrom] = np.zeros(
            (len(TSS_sites[chrom]), R[1] - R[0] + 1),
            dtype = bool
        ) - np.nan
        for i in range(len(TSS_sites[chrom])):
            if (TSS_strands[chrom][i] == "+"):
                S = scores[chrom][
                    np.logical_and(
                        scores[chrom][:,0] >= TSS_sites[chrom][i] + R[0],
                        scores[chrom][:,0] <= TSS_sites[chrom][i] + R[1]
                    ),
                    1
                ]
            else:
                S = scores[chrom][
                    np.logical_and(
                        scores[chrom][:,0] >= TSS_sites[chrom][i] - R[1],
                        scores[chrom][:,0] <= TSS_sites[chrom][i] - R[0]
                    ),
                    1
                ][::-1]
            
            if (len(S) == R[1] - R[0] + 1):
                output[chrom][i,:] = S
            
    return output

def get_TSS_scores(R, scores, TSS_sites, TSS_strands):
    total = np.arange(R[0],R[1]+1)
    N = len(total)
    total = np.column_stack((total, np.zeros(N)))
    n = 0
    
    for chrom in TSS_sites:
        for i in range(len(TSS_sites[chrom])):
            if (TSS_strands[chrom][i] == "+"):
                S = scores[chrom][
                    np.logical_and(
                        scores[chrom][:,0] >= TSS_sites[chrom][i] + R[0],
                        scores[chrom][:,0] <= TSS_sites[chrom][i] + R[1]
                    ),
                    1
                ]
            else:
                S = scores[chrom][
                    np.logical_and(
                        scores[chrom][:,0] >= TSS_sites[chrom][i] - R[1],
                        scores[chrom][:,0] <= TSS_sites[chrom][i] - R[0]
                    ),
                    1
                ][::-1]
            
            if (len(S) == N):
                n += 1
                total[:,1] += S
    
    total[:,1] /= n
    
    return total

def get_TSS_median_scores(R, scores, TSS_sites, TSS_strands):
    locs = np.arange(R[0],R[1]+1)
    N = len(locs)
    
    dTable = np.zeros((N, 0))
    
    for chrom in TSS_sites:
        for i in range(len(TSS_sites[chrom])):
            if (TSS_strands[chrom][i] == "+"):
                S = scores[chrom][
                    np.logical_and(
                        scores[chrom][:,0] >= TSS_sites[chrom][i] + R[0],
                        scores[chrom][:,0] <= TSS_sites[chrom][i] + R[1]
                    ),
                    1
                ]
            else:
                S = scores[chrom][
                    np.logical_and(
                        scores[chrom][:,0] >= TSS_sites[chrom][i] - R[1],
                        scores[chrom][:,0] <= TSS_sites[chrom][i] - R[0]
                    ),
                    1
                ][::-1]
            
            if (len(S) == N):
                dTable = np.column_stack((dTable, S))
    
    med = np.column_stack((locs, np.median(dTable, axis = 1)))
    
    return med

def get_TSS_sum(R, counts, TSS_sites, TSS_strands):
    total = np.arange(R[0],R[1]+1)
    N = len(total)
    total = np.column_stack((total, np.zeros(N)))
    
    for chrom in TSS_sites:
        for i in range(len(TSS_sites[chrom])):
            if (TSS_strands[chrom][i] == "+"):
                S = counts[chrom][
                    np.logical_and(
                        counts[chrom][:,0] >= TSS_sites[chrom][i] + R[0],
                        counts[chrom][:,0] <= TSS_sites[chrom][i] + R[1]
                    ),
                    1
                ]
            else:
                S = counts[chrom][
                    np.logical_and(
                        counts[chrom][:,0] >= TSS_sites[chrom][i] - R[1],
                        counts[chrom][:,0] <= TSS_sites[chrom][i] - R[0]
                    ),
                    1
                ][::-1]
            
            if (len(S) == N):
                total[:,1] += S
    
    return total

# Variables
R = (-500, 650) # The range of values around the TSS to graph

#%% Load Data

print("> Loading data...")

print(">> TSSs...")
TSS = io.load_tsv(
    "PreprocessedData/TSS",
    io.load_list("PreprocessedData/TSS/chromosomes.txt"),
    dtype = str
)

TSS_sites = an.make_type(an.column_split(TSS)[0], int)
TSS_strands = an.column_split(TSS)[4]

print(">> Genome and blocks...")
chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)

print(">> Scores and placements...")
dirs = [
    # "NucleosomePattern/cl_rt_log2",
    # "NucleosomePattern/cl_ice_log2",
    # "NucleosomePattern/cl_ice_abs_diff",
    "NucleosomePattern/composite_models/composite_ice"
]
num_dirs = len(dirs)

scores=[io.load_tsv(dirs[i]+"/genome_scores",chroms) for i in range(num_dirs)]

greedy = [io.load_tsv(D + "/dyad_calls", chroms, dtype = int) for D in dirs]
viterbi=[io.load_tsv(D+"/viterbi_dyad_calls",chroms,dtype=int) for D in dirs]
viterbi_overlap=[
    io.load_tsv(D+"/viterbi_40b_overlap_dyad_calls",chroms,dtype=int)
    for D in dirs
]
viterbi_overlap2=[
    io.load_tsv(D+"/viterbi_20b_overlap_dyad_calls",chroms,dtype=int)
    for D in dirs
]
# probs=[io.load_tsv(D+"/dynamic_probability_distribution",chroms) for D in dirs]

print(">> Weiner et al. data...")
weinerChroms=io.load_list("PreprocessedData/OtherData/weiner/chromosomes.txt")
weiner_dyads = an.column_split(
    io.load_tsv("PreprocessedData/OtherData/weiner", weinerChroms)
)[2]
weiner_dyads = an.make_type(weiner_dyads, int)
order = an.dict_argsort(weiner_dyads)
weiner_dyads = an.dict_index(weiner_dyads, order)

print(">> Weiner et al. reads...")
readTrials = ["1_1", "1_2", "1_3", "2_1", "2_2", "2_3", "3_1", "3_2"]
weiner_reads = [
    io.load_tsv(
        "PreprocessedData/OtherData/WeinerReads_" + trial,
        io.load_list(
            "PreprocessedData/OtherData/WeinerReads_"+trial+"/chromosomes.txt"
        )
    )
    for trial in readTrials
]

print(">> Brogaard et al. data...")
strong_dyads = an.make_type(an.column_split(
    io.load_tsv(
        "PreprocessedData/nucleosomes",
        io.load_list("PreprocessedData/nucleosomes/chromosomes.txt")
    )
)[2], int)

decent_dyads = an.make_type(an.column_split(
    io.load_tsv(
        "PreprocessedData/decent_nucleosomes",
        io.load_list("PreprocessedData/decent_nucleosomes/chromosomes.txt")
    )
)[2], int)

all_dyads = an.make_type(an.column_split(
    io.load_tsv(
        "PreprocessedData/all_nucleosomes",
        io.load_list("PreprocessedData/all_nucleosomes/chromosomes.txt")
    )
)[2], int)

#%% Main Loop

for i in range(num_dirs):
    print("> " + dirs[i] + "...")
    out_dir = dirs[i] + "/TSS_evals"
    if (not os.path.isdir(out_dir)):
        os.mkdir(out_dir)
    
    print(">> Greedy...")
    frequencies = get_TSS_frequencies(R, greedy[i], TSS_sites, TSS_strands)
    np.savetxt(
        out_dir + "/greedy_frequencies.tsv",
        frequencies,
        fmt = "%d",
        delimiter = "\t"
    )
    individual = get_all_TSS_dyads(R, greedy[i], TSS_sites, TSS_strands)
    io.output_tsv(
        individual,
        out_dir + "/greedy_individual",
        fmt = "%d",
        chroms = True
    )
    
    print(">> Viterbi...")
    frequencies = get_TSS_frequencies(R, viterbi[i], TSS_sites, TSS_strands)
    np.savetxt(
        out_dir + "/viterbi_frequencies.tsv",
        frequencies,
        fmt = "%d",
        delimiter = "\t"
    )
    individual = get_all_TSS_dyads(R, viterbi[i], TSS_sites, TSS_strands)
    io.output_tsv(
        individual,
        out_dir + "/viterbi_individual",
        fmt = "%d",
        chroms = True
    )
    
    print(">> Viterbi with 40bp Overlap...")
    frequencies = get_TSS_frequencies(R, viterbi_overlap[i], TSS_sites, TSS_strands)
    np.savetxt(
        out_dir + "/viterbi_40b_overlap_frequencies.tsv",
        frequencies,
        fmt = "%d",
        delimiter = "\t"
    )
    individual = get_all_TSS_dyads(R, viterbi_overlap[i], TSS_sites, TSS_strands)
    io.output_tsv(
        individual,
        out_dir + "/viterbi_40b_overlap_individual",
        fmt = "%d",
        chroms = True
    )
    
    print(">> Viterbi with 20bp Overlap...")
    frequencies = get_TSS_frequencies(R, viterbi_overlap2[i], TSS_sites, TSS_strands)
    np.savetxt(
        out_dir + "/viterbi_20b_overlap_frequencies.tsv",
        frequencies,
        fmt = "%d",
        delimiter = "\t"
    )
    individual = get_all_TSS_dyads(R, viterbi_overlap2[i], TSS_sites, TSS_strands)
    io.output_tsv(
        individual,
        out_dir + "/viterbi_20b_overlap_individual",
        fmt = "%d",
        chroms = True
    )
    
    print(">> Scores...")
    S = get_TSS_scores(R, scores[i], TSS_sites, TSS_strands)
    np.savetxt(out_dir + "/mean_score.tsv", S, delimiter = "\t")
    all_S = get_all_TSS_scores(R, scores[i], TSS_sites, TSS_strands)
    io.output_tsv(all_S, out_dir + "/individual_scores", chroms = True)
    median_S = get_TSS_median_scores(R, scores[i], TSS_sites, TSS_strands)
    np.savetxt(out_dir + "/median_score.tsv", median_S, delimiter = "\t")
    
    # print(">> Probabilities...")
    # P = get_TSS_scores(R, probs[i], TSS_sites, TSS_strands)
    # np.savetxt(out_dir + "/mean_prob.tsv", P, delimiter = "\t")
    # all_P = get_all_TSS_scores(R, probs[i], TSS_sites, TSS_strands)
    # io.output_tsv(all_P, out_dir + "/individual_probs", chroms = True)

#%% Reference Datasets

print("> Reference datasets...")
out_dir = "NucleosomePattern/TSS_reference_datasets"
if (not os.path.isdir(out_dir)):
    os.mkdir(out_dir)

print(">> Brogaard et al. dyad frequencies...")
frequencies = get_TSS_frequencies(R, strong_dyads, TSS_sites, TSS_strands)
np.savetxt(
    out_dir + "/strong_frequencies.tsv",
    frequencies,
    fmt = "%d",
    delimiter = "\t"
)
individual = get_all_TSS_dyads(R, strong_dyads, TSS_sites, TSS_strands)
io.output_tsv(individual, out_dir + "/strong_individual", fmt = "%d")

frequencies = get_TSS_frequencies(R, decent_dyads, TSS_sites, TSS_strands)
np.savetxt(
    out_dir + "/decent_frequencies.tsv",
    frequencies,
    fmt = "%d",
    delimiter = "\t"
)
individual = get_all_TSS_dyads(R, decent_dyads, TSS_sites, TSS_strands)
io.output_tsv(individual, out_dir + "/decent_individual", fmt = "%d")

frequencies = get_TSS_frequencies(R, all_dyads, TSS_sites, TSS_strands)
np.savetxt(
    out_dir + "/all_frequencies.tsv",
    frequencies,
    fmt = "%d",
    delimiter = "\t"
)
individual = get_all_TSS_dyads(R, all_dyads, TSS_sites, TSS_strands)
io.output_tsv(individual, out_dir + "/all_individual", fmt = "%d")

print(">> Weiner et al. dyad frequencies...")
frequencies = get_TSS_frequencies(R, weiner_dyads, TSS_sites, TSS_strands)
np.savetxt(
    out_dir + "/weiner_frequencies.tsv",
    frequencies,
    fmt = "%d",
    delimiter = "\t"
)
individual = get_all_TSS_dyads(R, weiner_dyads, TSS_sites, TSS_strands)
io.output_tsv(individual, out_dir + "/weiner_individual", fmt = "%d")

print(">> Weiner et al. read coverage...")
coverage = [
    get_TSS_sum(R, readData, TSS_sites, TSS_strands)
    for readData in weiner_reads
]

total_coverage = np.zeros(np.shape(coverage[0]))
total_coverage[:,0] = np.arange(R[0], R[1]+1)
total_coverage[:,1]=np.sum(np.vstack(tuple([C[:,1] for C in coverage])),axis=0)
np.savetxt(
    out_dir + "/read_coverage.tsv",
    total_coverage,
    delimiter = "\t"
)


#%% End

print("===== 082_transcription_start_sites.py - Exiting Properly =====")