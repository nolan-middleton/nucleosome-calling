#############################################################################
# This file will take all of the model scores generated and evaluate their  #
# performance by calculating different model evaluation statistics.         #
#############################################################################

print("===== 052_evaluate_model_performance.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an

# Test thresholds
thresholds = [0, 0.5, 1, 2]

# Datasets to load
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

# The number of background models used
N = nuc.get_N()

# The test points we evaluated
test_points = np.loadtxt(
    "evaluationPositions.tsv",
    dtype = str,
    delimiter = "\t"
)
test_points_cols = test_points[0,:]
test_points = test_points[1:,:]

answers = test_points[:,test_points_cols == "answer"][:,0] == "True"

#%% Load Data

print("> Loading data...")
model_scores = {}
for dataset in datasets:
    models_dir = "NucleosomePattern/" + dataset + "/models"
    
    model_scores[dataset] = {}
    for model in models:
        model_scores[dataset][model] = np.loadtxt(
            models_dir + "/" + model + "/evaluationScores.tsv",
            dtype = float,
            delimiter = "\t",
            skiprows = 1
        )

#%% Model Statistics

print("> Model statistics...")
stats_table = np.array(
    [
        [
            "dataset",
            "model",
            "background",
            "threshold",
            "statistic",
            "value"
        ]
    ],
    dtype = "<U32"
)

stats = [
    "TP",
    "TN",
    "FP",
    "FN",
    "accuracy",
    "sensitivity",
    "specificity",
    "precision",
    "NPV",
    "phi"
]

for dataset in datasets:
    # Load data
    print(">> " + dataset + "...")
    
    scores = model_scores[dataset]
    
    # Calculate statistics for each threshold
    model_stats = {}
    for model in scores:
        model_stats[model] = {}
        for threshold in thresholds:
            TP = np.sum(
                scores[model][answers] > threshold, axis = 0
            ).astype("int64")
            TN = np.sum(
                scores[model][~answers] <= threshold,
                axis = 0
            ).astype("int64")
            FP = np.sum(
                scores[model][~answers] > threshold,
                axis = 0
            ).astype("int64")
            FN = np.sum(
                scores[model][answers] <= threshold,
                axis = 0
            ).astype("int64")
            model_stats[model][str(threshold)] = {
                "TP": TP,
                "TN": TN,
                "FP": FP,
                "FN": FN,
                "accuracy": (TP + TN)/(TP + TN + FP + FN),
                "sensitivity": TP / (TP + FN),
                "specificity": TN / (TN + FP),
                "precision": TP / (TP + FP),
                "NPV": TN / (TN + FN),
                "phi": (TP*TN-FP*FN)/np.sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
            }
    
    stats_table = np.vstack(
        (
            stats_table,
            np.column_stack(
                (
                    np.repeat(
                        dataset,
                        len(models) * N * len(thresholds) * len(stats)
                    ),
                    np.repeat(
                        models,
                        N * len(thresholds) * len(stats)
                    ),
                    np.tile(
                        np.repeat(range(N), len(thresholds) * len(stats)),
                        len(models)
                    ),
                    np.tile(
                        np.repeat(thresholds, len(stats)),
                        len(models) * N
                    ),
                    np.tile(
                        stats,
                        len(models) * N * len(thresholds)
                    ),
                    [
                        model_stats[model][str(threshold)][stat][n]
                        for model in models
                        for n in range(N)
                        for threshold in thresholds
                        for stat in stats
                    ]
                )
            )
        )
    )

print(">> Outputting...")
np.savetxt(
    "evaluationTable.tsv",
    stats_table,
    fmt = "%s",
    delimiter = "\t"
)

#%% ROC and PR Curves

print("> ROC and PR curves...")
curves_table = np.array(
    [
        [
            "dataset",
            "model",
            "background",
            "precision",
            "sensitivity",
            "FPR"
        ]
    ]
)

for dataset in datasets:
    print(">> " + dataset + "...")
    for model in models:
        print(">>> " + model + "...")
        for n in range(N):
            print("=> " + str(n) + "...")
            theseScores = model_scores[dataset][model][:,n]
            thesePoints = np.unique(theseScores)
            L = len(thesePoints)
            
            precisions = np.zeros(L)
            sensitivities = np.zeros(L)
            FPRs = np.zeros(L)
            for i in range(L):
                thresh = thesePoints[i]
                TP = np.sum(theseScores[answers] > thresh)
                TN = np.sum(theseScores[~answers] <= thresh)
                FP = np.sum(theseScores[~answers] > thresh)
                FN = np.sum(theseScores[answers] <= thresh)
                
                if (TP == 0):
                    precisions[i] = 0
                    sensitivities[i] = 0
                else:
                    precisions[i] = TP / (TP + FP)
                    sensitivities[i] = TP / (TP + FN)
                if (FP == 0):
                    FPRs[i] = 0
                else:
                    FPRs[i] = FP / (TN + FP)
            
            curves_table = np.vstack(
                (
                    curves_table,
                    np.column_stack(
                        (
                            np.repeat(dataset, L),
                            np.repeat(model, L),
                            np.repeat(n, L),
                            precisions,
                            sensitivities,
                            FPRs
                        )
                    )
                )
            )

print(">> Outputting...")
np.savetxt(
    "ROC_PR_curveTable.tsv",
    curves_table,
    fmt = "%s",
    delimiter = "\t"
)

#%% End

print("===== 052_evaluate_model_performance.py - Exiting Properly =====")