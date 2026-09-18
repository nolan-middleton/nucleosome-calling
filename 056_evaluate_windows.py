#############################################################################
# This file will evaluate the scores given to all of the models across a    #
# 20bp window, finding peaks and deciding if the algorithm correctly placed #
# the nucleosome.                                                           #
#############################################################################

print("===== 054_evaluate_windows.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

def find_peaks(Y):
    peaks = []
    if (Y[0] > Y[1]):
        peaks.append(0)
    
    for i in range(1, len(Y) - 1):
        if (Y[i-1] < Y[i]) and (Y[i] > Y[i+1]):
            peaks.append(i)
    
    if (Y[-1] > Y[-2]):
        peaks.append(len(Y) - 1)
    
    return peaks

def get_class(peaks):
    return np.min(np.abs(peaks))

datasets = [
    "unfloored_cl",
    "cl_rt_log2",
    "cl_ice_log2",
    "cl_rt_diff",
    "cl_ice_diff",
    "cl_rt_abs_diff",
    "cl_ice_abs_diff"
]

models = ["flat", "linear", "quadratic"]

#%% Load Data

print("> Loading data...")

test_positions = np.loadtxt(
    "evaluationPositions.tsv",
    delimiter = "\t",
    dtype = str,
    skiprows = 1
)
M = np.shape(test_positions)[0]

window_scores = {}
for dataset in datasets:
    window_scores[dataset] = {}
    for model in models:
        this_dir = "NucleosomePattern/" + dataset + "/models/" + model
        window_scores[dataset][model] = np.column_stack(
            tuple(
                [
                    np.loadtxt(
                        this_dir+"/windowEvalScores_"+str(offset)+".tsv",
                        delimiter = "\t",
                        skiprows = 1,
                        usecols = 0
                    )
                    for offset in range(-10, 11)
                ]
            )
        )

#%% Find Peaks

peaks = {}
for dataset in window_scores:
    peaks[dataset] = {}
    for model in window_scores[dataset]:
        peaks[dataset][model] = [
            np.array(find_peaks(window_scores[dataset][model][i,:])) - 10
            for i in range(M)
        ]

#%% Classes

classes = {}
for dataset in peaks:
    classes[dataset] = {}
    for model in peaks[dataset]:
        classes[dataset][model] = np.array(
            [get_class(peaklist) for peaklist in peaks[dataset][model]]
        )

#%% Outputs

for dataset in classes:
    np.savetxt(
        "NucleosomePattern/" + dataset + "/models/windowClasses.tsv",
        np.column_stack(
            tuple([classes[dataset][model] for model in classes[dataset]])
        ),
        delimiter = "\t",
        fmt = "%d",
        header = "".join([str(model)+"\t" for model in classes[dataset]])[:-1],
        comments = ""
    )

#%% End

print("===== 054_evaluate_windows.py - Exiting Properly =====")