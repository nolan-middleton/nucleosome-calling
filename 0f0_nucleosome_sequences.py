#############################################################################
# This file will find sequences pertaining to the nucleosomes we placed.    #
#############################################################################

print("===== 0f0_nucleosome_sequences.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Defs
def do_sequence_analysis(dyads, genome, dinucs, directory, W = 73):
    positions = {}
    nuc_dinucs = {}
    seqs = {}
    for chrom in dyads:
        positions[chrom] = np.zeros((len(dyads[chrom]), 2*W), float)
        for i in range(len(dyads[chrom])):
            positions[chrom][i,:] = np.arange(
                dyads[chrom][i] - W + 0.5,
                dyads[chrom][i] + W
            )
        lefts = (positions[chrom] - 0.5).astype(int)
        
        valid = (lefts >= 0) & (lefts < len(dinucs[chrom]) - 1)
        lefts[~valid] = 0
        
        nuc_dinucs[chrom] = dinucs[chrom][lefts]
        nuc_dinucs[chrom][~valid] = ""
        seqs[chrom] = genome[chrom][np.column_stack((lefts, lefts[:,-1] + 1))]
        seqs[chrom] = np.column_stack((dyads[chrom], seqs[chrom]))
        seqs[chrom][
            ~np.column_stack(
                (
                    np.ones(len(dyads[chrom]), dtype = bool),
                    valid,
                    valid[:,-1]
                )
            )
        ] = ""

    dinuc_table = an.dict_stack(nuc_dinucs)

    WW_freqs = np.sum(dinuc_table == "AA", axis = 0) + \
        np.sum(dinuc_table == "AT", axis = 0) + \
        np.sum(dinuc_table == "TA", axis = 0) + \
        np.sum(dinuc_table == "TT", axis = 0)
    
    SS_freqs = np.sum(dinuc_table == "GG", axis = 0) + \
        np.sum(dinuc_table == "GC", axis = 0) + \
        np.sum(dinuc_table == "CG", axis = 0) + \
        np.sum(dinuc_table == "CC", axis = 0)
    
    if (not os.path.isdir(directory)):
        os.mkdir(directory)
    
    io.output_tsv(seqs,directory+"/seqs",fmt="%s",chroms=True)
    np.savetxt(directory + "/WW.tsv", WW_freqs, delimiter = "\t")
    np.savetxt(directory + "/SS.tsv", SS_freqs, delimiter = "\t")
    
    return (WW_freqs, SS_freqs, seqs)

#%% Load Data

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = an.make_type(io.load_tsv("PreprocessedData/blocks", chroms), int)

directory = "NucleosomePattern/composite_models/composite_ice"

greedy = io.load_tsv(
    directory + "/dyad_calls",
    chroms,
    dtype = int
)

viterbi = io.load_tsv(
    directory + "/viterbi_dyad_calls",
    chroms,
    dtype = int
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

#%% Find dinucleotides

dinucs = {}
for chrom in genome:
    dinucs[chrom] = genome[chrom][:-1] + genome[chrom][1:]

#%% Get Sequences in All Nucleosomes

D = directory + "/all_sequence_analysis"

N = 0
for chrom in viterbi:
    N += len(viterbi[chrom])

all_WW, all_SS, all_seqs = do_sequence_analysis(viterbi, genome, dinucs, D)

with open(D + "/n.txt", "w") as file:
    file.write(str(N) + "\n")

#%% Get Sequences in Weakly Rotationally-Positioned Nucleosomes

D = directory + "/Q1_sequence_analysis"

Q = int(np.shape(rotational_scores)[0] / 4)

Q1=rotational_scores[np.argsort(rotational_scores[:,2].astype(float)),:][:Q]
Q1_dyads = {}
for chrom in np.unique(Q1[:,0]):
    Q1_dyads[str(chrom)] = np.sort(Q1[Q1[:,0] == chrom,1].astype(int))

Q1_WW, Q1_SS, Q1_seqs = do_sequence_analysis(Q1_dyads, genome, dinucs, D)

with open(D + "/n.txt", "w") as file:
    file.write(str(Q) + "\n")

#%% Get Sequences in Strongly Rotationally-Positioned Nucleosomes

D = directory + "/Q3_sequence_analysis"

Q = int(np.shape(rotational_scores)[0] / 4)

Q3=rotational_scores[np.argsort(rotational_scores[:,2].astype(float)),:][-Q:]
Q3_dyads = {}
for chrom in np.unique(Q3[:,0]):
    Q3_dyads[str(chrom)] = np.sort(Q3[Q3[:,0] == chrom,1].astype(int))

Q3_WW, Q3_SS, Q3_seqs = do_sequence_analysis(Q3_dyads, genome, dinucs, D)

with open(D + "/n.txt", "w") as file:
    file.write(str(Q) + "\n")

#%% Get Sequences in Very Weakly Rotationally-Positioned Nucleosomes

D = directory + "/P1_sequence_analysis"

P = int(np.shape(rotational_scores)[0] / 10)

P1=rotational_scores[np.argsort(rotational_scores[:,2].astype(float)),:][:P]
P1_dyads = {}
for chrom in np.unique(P1[:,0]):
    P1_dyads[str(chrom)] = np.sort(P1[P1[:,0] == chrom,1].astype(int))

P1_WW, P1_SS, P1_seqs = do_sequence_analysis(P1_dyads, genome, dinucs, D)

with open(D + "/n.txt", "w") as file:
    file.write(str(P) + "\n")

#%% Get Sequences in Very Strongly Rotationally-Positioned Nucleosomes

D = directory + "/P3_sequence_analysis"

P = int(np.shape(rotational_scores)[0] / 10)

P3=rotational_scores[np.argsort(rotational_scores[:,2].astype(float)),:][-P:]
P3_dyads = {}
for chrom in np.unique(P3[:,0]):
    P3_dyads[str(chrom)] = np.sort(P3[P3[:,0] == chrom,1].astype(int))

P3_WW, P3_SS, P3_seqs = do_sequence_analysis(P3_dyads, genome, dinucs, D)

with open(D + "/n.txt", "w") as file:
    file.write(str(P) + "\n")

#%% Get Sequences for Brogaard

D = "NucleosomePattern/ReferenceSequenceAnalysis"
if (not os.path.isdir(D)):
    os.mkdir(D)

D += "/brogaard_all"

for chrom in brogaard:
    brogaard[chrom] = brogaard[chrom][
        (brogaard[chrom] >= 73) & (brogaard[chrom] < len(genome[chrom]) - 73)
    ]

WW, SS, seqs = do_sequence_analysis(brogaard, genome, dinucs, D)

#%% Get Sequences for Weiner

D = "NucleosomePattern/ReferenceSequenceAnalysis"
if (not os.path.isdir(D)):
    os.mkdir(D)

D += "/weiner"

for chrom in weiner:
    weiner[chrom] = weiner[chrom][
        (weiner[chrom] >= 73) & (weiner[chrom] < len(genome[chrom]) - 73)
    ]

WW_weiner,SS_weiner,seqs_weiner = do_sequence_analysis(weiner,genome,dinucs,D)

#%% Get Sequences for Redundant

D = "NucleosomePattern/ReferenceSequenceAnalysis"
if (not os.path.isdir(D)):
    os.mkdir(D)

D += "/redund"

for chrom in weiner:
    redundant[chrom] = redundant[chrom][
        (redundant[chrom] >= 73) & (redundant[chrom] < len(genome[chrom])-73)
    ]

WW_redund,SS_redund,seqs_redund=do_sequence_analysis(redundant,genome,dinucs,D)

#%% Get Wide Sequences in Very Weakly Rotationally-Positioned Nucleosomes

D = directory + "/P1_sequence_analysis_wide"

P = int(np.shape(rotational_scores)[0] / 10)

P1=rotational_scores[np.argsort(rotational_scores[:,2].astype(float)),:][:P]
P1_dyads = {}
for chrom in np.unique(P1[:,0]):
    P1_dyads[str(chrom)] = np.sort(P1[P1[:,0] == chrom,1].astype(int))

P1_WW, P1_SS, P1_seqs = do_sequence_analysis(P1_dyads, genome, dinucs, D, 500)

with open(D + "/n.txt", "w") as file:
    file.write(str(P) + "\n")

#%% Get Wide Sequences in Very Strongly Rotationally-Positioned Nucleosomes

D = directory + "/P3_sequence_analysis_wide"

P = int(np.shape(rotational_scores)[0] / 10)

P3=rotational_scores[np.argsort(rotational_scores[:,2].astype(float)),:][-P:]
P3_dyads = {}
for chrom in np.unique(P3[:,0]):
    P3_dyads[str(chrom)] = np.sort(P3[P3[:,0] == chrom,1].astype(int))

P3_WW, P3_SS, P3_seqs = do_sequence_analysis(P3_dyads, genome, dinucs, D, 500)

with open(D + "/n.txt", "w") as file:
    file.write(str(P) + "\n")

#%% Get Wide Sequences in Weakly Rotationally-Positioned Nucleosomes

D = directory + "/Q1_sequence_analysis_wide"

Q = int(np.shape(rotational_scores)[0] / 4)

Q1=rotational_scores[np.argsort(rotational_scores[:,2].astype(float)),:][:Q]
Q1_dyads = {}
for chrom in np.unique(Q1[:,0]):
    Q1_dyads[str(chrom)] = np.sort(Q1[Q1[:,0] == chrom,1].astype(int))

Q1_WW, Q1_SS, Q1_seqs = do_sequence_analysis(Q1_dyads, genome, dinucs, D, 500)

with open(D + "/n.txt", "w") as file:
    file.write(str(Q) + "\n")

#%% Get Wide Sequences in Strongly Rotationally-Positioned Nucleosomes

D = directory + "/Q3_sequence_analysis_wide"

Q = int(np.shape(rotational_scores)[0] / 4)

Q3=rotational_scores[np.argsort(rotational_scores[:,2].astype(float)),:][-Q:]
Q3_dyads = {}
for chrom in np.unique(Q3[:,0]):
    Q3_dyads[str(chrom)] = np.sort(Q3[Q3[:,0] == chrom,1].astype(int))

Q3_WW, Q3_SS, Q3_seqs = do_sequence_analysis(Q3_dyads, genome, dinucs, D, 500)

with open(D + "/n.txt", "w") as file:
    file.write(str(Q) + "\n")

#%% Get Wide Sequences for All

D = directory + "/all_wide_sequence_analysis"

for chrom in viterbi:
    viterbi[chrom] = viterbi[chrom][
        (viterbi[chrom] >= 500) & (viterbi[chrom] < len(genome[chrom]) - 500)
    ]

N = 0
for chrom in viterbi:
    N += len(viterbi[chrom])

all_WW, all_SS, all_seqs = do_sequence_analysis(viterbi,genome,dinucs,D,500)

with open(D + "/n.txt", "w") as file:
    file.write(str(N) + "\n")

#%% Get Wide Sequences for Brogaard

D = "NucleosomePattern/WideReferenceSequenceAnalysis"
if (not os.path.isdir(D)):
    os.mkdir(D)

D += "/brogaard_all"

for chrom in brogaard:
    brogaard[chrom] = brogaard[chrom][
        (brogaard[chrom] >= 500) & (brogaard[chrom] < len(genome[chrom])-500)
    ]

WW, SS, seqs = do_sequence_analysis(brogaard, genome, dinucs, D, 500)

#%% Get Sequences for Weiner

D = "NucleosomePattern/WideReferenceSequenceAnalysis"
if (not os.path.isdir(D)):
    os.mkdir(D)

D += "/weiner"

for chrom in weiner:
    weiner[chrom] = weiner[chrom][
        (weiner[chrom] >= 500) & (weiner[chrom] < len(genome[chrom]) - 500)
    ]

WW_wein,SS_wein,seqs_wein = do_sequence_analysis(weiner,genome,dinucs,D,500)

#%% Get Sequences for Redundant

D = "NucleosomePattern/WideReferenceSequenceAnalysis"
if (not os.path.isdir(D)):
    os.mkdir(D)

D += "/redund"

for chrom in weiner:
    redundant[chrom] = redundant[chrom][
        (redundant[chrom] >= 500) & (redundant[chrom] < len(genome[chrom])-500)
    ]

WW_red,SS_red,seqs_red = do_sequence_analysis(redundant,genome,dinucs,D,500)

#%% End

print("===== 0f0_nucleosome_sequences.py - Exiting Properly =====")