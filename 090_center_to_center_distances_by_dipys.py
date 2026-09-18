#############################################################################
# This file will evaluate the nucleosome placement based on the center to   #
# center distances from the predicted nucleosomes to the real nucleosomes.  #
#############################################################################

print("===== 090_center_to_center_distances_by_dipys.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

def get_D(putative_dyads, reference_dyads):
    this_D = an.get_distance_matrix(
        putative_dyads,
        reference_dyads,
        signed = True
    )
    selected = an.dict_argmin(
        an.dict_abs(this_D),
        axis = 1
    )
    
    return_dict = {}
    for chrom in selected:
        return_dict[chrom] = this_D[chrom][
            np.arange(0,np.shape(this_D[chrom])[0]),
            selected[chrom]
        ]
    
    return return_dict

# Variables
distance_threshold = 65

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

# Load in test nucleosome placements
print(">> Putative dyad calls...")
putative_dyads = [
    io.load_tsv(
        D + "/dyad_calls",
        chromosomes = io.load_list(D + "/dyad_calls/chromosomes.txt"),
        dtype = int
    )
    for D in dirs
] + [
    io.load_tsv(
        D + "/viterbi_dyad_calls",
        chromosomes = io.load_list(D + "/dyad_calls/chromosomes.txt"),
        dtype = int
    )
    for D in dirs
]
putative_dyad_names = ["greedy" for D in dirs] + ["viterbi" for D in dirs]
putative_dyad_dirs = [D + "/distance_evals" for D in dirs]
putative_dyad_dirs += [D + "/distance_evals" for D in dirs]

print(">> Reference dyads...")
ref_dirs = [
    "PreprocessedData/nucleosomes",
    "PreprocessedData/decent_nucleosomes",
    "PreprocessedData/all_nucleosomes"
]
reference_dyads = [
    an.make_type(
        an.column_split(
            io.load_tsv(D, io.load_list(D + "/chromosomes.txt"))
        )[2],
        int
    )
    for D in ref_dirs
]
reference_dyad_names = ["strong", "decent", "all"]

print(">> Weiner et al. data...")
weinerChroms=io.load_list("PreprocessedData/OtherData/weiner/chromosomes.txt")
weiner_dyads = an.column_split(
    io.load_tsv("PreprocessedData/OtherData/weiner", weinerChroms)
)[2]
weiner_dyads = an.make_type(weiner_dyads, int)
order = an.dict_argsort(weiner_dyads)
weiner_dyads = an.dict_index(weiner_dyads, order)

putative_dyads += [weiner_dyads]
putative_dyad_names += ["Weiner"]
putative_dyad_dirs += ["NucleosomePattern/weiner_brogaard_distance_evals"]
reference_dyads += [weiner_dyads]
reference_dyad_names += ["Weiner"]

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

#%% Getting Nucleosome Sequences

reference_seqs = []
reference_dinucs = []
total_dinucs = []
Q1_dyads = []
Q3_dyads = []
for nucs in reference_dyads:
    reference_seqs.append({})
    reference_dinucs.append({})
    for chrom in nucs:
        starts = nucs[chrom] - 73
        ends = nucs[chrom] + 73 + 1
        
        valid = (starts >= 0) & (ends < len(genome[chrom]))
        nucs[chrom] = nucs[chrom][valid]
        
        reference_seqs[-1][chrom] = np.zeros(
            (np.shape(nucs[chrom])[0], 147),
            "<U1"
        )
        for j in range(147):
            reference_seqs[-1][chrom][:,j] = genome[chrom][nucs[chrom]-73+j]
        
        reference_dinucs[-1][chrom] = np.zeros(np.shape(nucs[chrom])[0], int)
        for j in range(146):
            reference_dinucs[-1][chrom] += np.isin(
                reference_seqs[-1][chrom][:,j]+reference_seqs[-1][chrom][:,j+1],
                ["TT", "TC", "CT", "CC", "AA", "GA", "AG", "GG"]
            )
    total_dinucs.append(
        np.concatenate(
            tuple(
                [reference_dinucs[-1][chrom] for chrom in reference_dinucs[-1]]
            )
        )
    )
    
    Q1 = np.quantile(total_dinucs[-1], 0.25)
    Q3 = np.quantile(total_dinucs[-1], 0.75)
    
    Q1_dyads.append({})
    Q3_dyads.append({})
    for chrom in nucs:
        Q1_dyads[-1][chrom] = nucs[chrom][reference_dinucs[-1][chrom] <= Q1]
        Q3_dyads[-1][chrom] = nucs[chrom][reference_dinucs[-1][chrom] >= Q3]

for i in range(len(reference_dinucs)):
    io.output_tsv(
        reference_dinucs[i],
        "NucleosomePattern/" + reference_dyad_names[i] + "_dinucs",
        fmt = "%d"
    )
    
    io.output_json(
        "NucleosomePattern/" + reference_dyad_names[i] + "_dinucs_stats.json",
        {
            "mean": np.mean(total_dinucs[i]),
            "std": np.std(total_dinucs[i]),
            "Q1": np.quantile(total_dinucs[i], 0.25),
            "Q3": np.quantile(total_dinucs[i], 0.75),
            "n": len(total_dinucs[i]),
            "med": np.quantile(total_dinucs[i], 0.5),
            "data": total_dinucs[i].tolist()
        }
    )

#%% Main

print("> Getting distances...")
Q1_distances = [[get_D(P, R) for P in putative_dyads] for R in Q1_dyads]
Q3_distances = [[get_D(P, R) for P in putative_dyads] for R in Q3_dyads]

print("> Getting frequency tables...")
Q1_freq_tables = [
    [nuc.get_distance_freq_table(P, R) for P in putative_dyads]
    for R in Q1_dyads
]
Q3_freq_tables = [
    [nuc.get_distance_freq_table(P, R) for P in putative_dyads]
    for R in Q3_dyads
]

D_Q1_freq_tables = [
    [nuc.get_distance_freq_table_from_distances(D) for D in L]
    for L in Q1_distances
]
D_Q3_freq_tables = [
    [nuc.get_distance_freq_table_from_distances(D) for D in L]
    for L in Q3_distances
]

print("> Getting periodic proportions...")
Q1_periodic_proportions = [
    [nuc.get_periodic_proportion(D, rotational_positions) for D in L]
    for L in D_Q1_freq_tables
]
Q3_periodic_proportions = [
    [nuc.get_periodic_proportion(D, rotational_positions) for D in L]
    for L in D_Q3_freq_tables
]
Q1_anti_periodic_proportions = [
    [nuc.get_periodic_proportion(D, anti_rotational_positions) for D in L]
    for L in D_Q1_freq_tables
]
Q3_anti_periodic_proportions = [
    [nuc.get_periodic_proportion(D, anti_rotational_positions) for D in L]
    for L in D_Q3_freq_tables
]

E = np.sum(
    np.abs(rotational_positions) <= distance_threshold
) / (2*distance_threshold + 1)
anti_E = np.sum(
    np.abs(anti_rotational_positions) <= distance_threshold
) / (2*distance_threshold + 1)

#%% Outputs

print("> Outputting...")

print(">> Initial setup...")
for i in range(len(putative_dyads)):
    out_dir = putative_dyad_dirs[i]
    if (not os.path.isdir(out_dir)):
        os.mkdir(out_dir)
    
    putName = putative_dyad_names[i]
    fname = out_dir + "/Q1_periodic_proportions.txt"
    with open(fname, "w") as file:
        file.write("# p=" + str(E) + "\n")
    fname = out_dir + "/Q3_periodic_proportions.txt"
    with open(fname, "w") as file:
        file.write("# p=" + str(E) + "\n")
    fname = out_dir + "/Q1_anti_periodic_proportions.txt"
    with open(fname, "w") as file:
        file.write("# p=" + str(anti_E) + "\n")
    fname = out_dir + "/Q3_anti_periodic_proportions.txt"
    with open(fname, "w") as file:
        file.write("# p=" + str(anti_E) + "\n")

print(">> Periodic proportions...")
for i in range(len(putative_dyads)):
    out_dir = putative_dyad_dirs[i]
    putName = putative_dyad_names[i]
    fname = out_dir + "/Q1_periodic_proportions.txt"
    with open(fname,"a") as file:
        for I in range(len(Q1_periodic_proportions)):
            refName = reference_dyad_names[I]
            if (refName != "Weiner") or (putName != "Weiner"):
                file.write(refName+"\t"+putName+"\t")
                file.write(str(Q1_periodic_proportions[I][i][0])+"\t")
                file.write(str(Q1_periodic_proportions[I][i][1])+"\n")
    fname = out_dir + "/Q3_periodic_proportions.txt"
    with open(fname,"a") as file:
        for I in range(len(Q3_periodic_proportions)):
            refName = reference_dyad_names[I]
            if (refName != "Weiner") or (putName != "Weiner"):
                file.write(refName+"\t"+putName+"\t")
                file.write(str(Q3_periodic_proportions[I][i][0])+"\t")
                file.write(str(Q3_periodic_proportions[I][i][1])+"\n")
    fname = out_dir + "/Q1_anti_periodic_proportions.txt"
    with open(fname,"a") as file:
        for I in range(len(Q1_anti_periodic_proportions)):
            refName = reference_dyad_names[I]
            if (refName != "Weiner") or (putName != "Weiner"):
                file.write(refName+"\t"+putName+"\t")
                file.write(str(Q1_anti_periodic_proportions[I][i][0])+"\t")
                file.write(str(Q1_anti_periodic_proportions[I][i][1])+"\n")
    fname = out_dir + "/Q3_anti_periodic_proportions.txt"
    with open(fname,"a") as file:
        for I in range(len(Q3_anti_periodic_proportions)):
            refName = reference_dyad_names[I]
            if (refName != "Weiner") or (putName != "Weiner"):
                file.write(refName+"\t"+putName+"\t")
                file.write(str(Q3_anti_periodic_proportions[I][i][0])+"\t")
                file.write(str(Q3_anti_periodic_proportions[I][i][1])+"\n")

for I in range(len(Q1_distances)):
    refName = reference_dyad_names[I]
    print(">> " + refName + "...")
    for i in range(len(putative_dyads)):
        putName = putative_dyad_names[i]
        print(">>> " + putName + "...")
        
        out_dir = putative_dyad_dirs[i]
        
        if (refName != "Weiner") or (putName != "Weiner"):
            io.output_tsv(
                Q1_distances[I][i],
                out_dir+"/"+refName+"_"+putName+"_Q1_distances",
                fmt = "%d",
                chroms = True
            )
            np.savetxt(
                out_dir+"/"+refName+"_"+putName+"_Q1_freq_table.tsv",
                Q1_freq_tables[I][i],
                fmt = "%d",
                delimiter = "\t"
            )
            np.savetxt(
                out_dir+"/"+refName+"_"+putName+"_Q1_offset_freq_table.tsv",
                D_Q1_freq_tables[I][i],
                fmt = "%d",
                delimiter = "\t"
            )
            io.output_tsv(
                Q3_distances[I][i],
                out_dir+"/"+refName+"_"+putName+"_Q3_distances",
                fmt = "%d",
                chroms = True
            )
            np.savetxt(
                out_dir+"/"+refName+"_"+putName+"_Q3_freq_table.tsv",
                Q3_freq_tables[I][i],
                fmt = "%d",
                delimiter = "\t"
            )
            np.savetxt(
                out_dir+"/"+refName+"_"+putName+"_Q3_offset_freq_table.tsv",
                D_Q3_freq_tables[I][i],
                fmt = "%d",
                delimiter = "\t"
            )

#%% End

print(
    "===== 090_center_to_center_distances_by_dipys.py - Exiting Properly ====="
)