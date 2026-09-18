#############################################################################
# This file will evaluate the nucleosome placement based on the center to   #
# center distances from the predicted nucleosomes to the real nucleosomes.  #
#############################################################################

print("===== 0f1_center_to_center_a_minus_3.py - Starting =====")

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

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = an.make_type(io.load_tsv("PreprocessedData/blocks", chroms), int)

directory = "NucleosomePattern/composite_models/composite_ice"

viterbi = io.load_tsv(
    directory + "/all_sequence_analysis/seqs",
    chroms,
    dtype = str
)

weinerChroms=io.load_list("PreprocessedData/OtherData/weiner/chromosomes.txt")
weiner = an.column_split(
    io.load_tsv("PreprocessedData/OtherData/weiner", weinerChroms)
)[2]
weiner = an.make_type(weiner, int)
order = an.dict_argsort(weiner)
weiner = an.dict_index(weiner, order)

brogaardChroms=io.load_list("PreprocessedData/all_nucleosomes/chromosomes.txt")
brogaard = an.make_type(
    an.column_split(
        io.load_tsv("PreprocessedData/all_nucleosomes",brogaardChroms)
    )[2],
    int
)

redund = np.loadtxt(
    "../Data/Brogaard_RedundantDyads.tsv",
    dtype = str,
    delimiter = "\t"
)
redundant = {}
for chrom in np.unique(redund[:,0]):
    redundant[str(chrom)] = redund[redund[:,0] == chrom, 1].astype(int)

rotational_scores = np.loadtxt(
    directory + "/rotational_scores.tsv",
    delimiter = "\t",
    dtype = str
)

#%% Getting Nucleosome Sequences

AT3 = {}
not_AT3 = {}
for chrom in viterbi:
    selected = (viterbi[chrom][:,71] == "A") | (viterbi[chrom][:,77] == "T")
    AT3[chrom] = viterbi[chrom][selected, 0].astype(int)
    not_AT3[chrom] = viterbi[chrom][~selected, 0].astype(int)

#%% Main

print("> Getting distances...")
brogaard_AT3_dists = nuc.get_distance_freq_table_from_distances(
    get_D(brogaard, AT3)
)
weiner_AT3_dists = nuc.get_distance_freq_table_from_distances(
    get_D(weiner, AT3)
)
redund_AT3_dists = nuc.get_distance_freq_table_from_distances(
    get_D(redundant, AT3)
)

brogaard_notAT3_dists = nuc.get_distance_freq_table_from_distances(
    get_D(brogaard, not_AT3)
)
weiner_notAT3_dists = nuc.get_distance_freq_table_from_distances(
    get_D(weiner, not_AT3)
)
redund_notAT3_dists = nuc.get_distance_freq_table_from_distances(
    get_D(redundant, not_AT3)
)

#%% Outputs

print("> Outputting...")
if (not os.path.isdir(directory + "/reverse_AT3_distance_evals")):
    os.mkdir(directory + "/reverse_AT3_distance_evals")

np.savetxt(
    directory + "/reverse_AT3_distance_evals/brogaard_AT3.tsv",
    brogaard_AT3_dists,
    delimiter = "\t",
    fmt = "%d"
)
np.savetxt(
    directory + "/reverse_AT3_distance_evals/weiner_AT3.tsv",
    weiner_AT3_dists,
    delimiter = "\t",
    fmt = "%d"
)
np.savetxt(
    directory + "/reverse_AT3_distance_evals/redund_AT3.tsv",
    redund_AT3_dists,
    delimiter = "\t",
    fmt = "%d"
)

np.savetxt(
    directory + "/reverse_AT3_distance_evals/brogaard_notAT3.tsv",
    brogaard_notAT3_dists,
    delimiter = "\t",
    fmt = "%d"
)
np.savetxt(
    directory + "/reverse_AT3_distance_evals/weiner_notAT3.tsv",
    weiner_notAT3_dists,
    delimiter = "\t",
    fmt = "%d"
)
np.savetxt(
    directory + "/reverse_AT3_distance_evals/redund_notAT3.tsv",
    redund_notAT3_dists,
    delimiter = "\t",
    fmt = "%d"
)

#%% End

print("===== 0f1_center_to_center_a_minus_3.py - Exiting Properly =====")