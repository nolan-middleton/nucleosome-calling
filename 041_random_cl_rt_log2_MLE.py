#############################################################################
# This file will find the maximum likelihood estimates for the nucleosome   #
# random models assuming a constant baseline likelihood.                    #
#############################################################################

print("===== 041_random_cl_rt_log2_MLE.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import Analysis as an
import Nucleosomes as nuc
import DataIO as io
import os as os

# Defs
def flat_mean(A):
    return np.zeros(146) + A[0]

def std_relation(A, means):
    '''
    A = [k]
    
    sigma = k
    '''
    return np.zeros(146) + A[0]

dataset_name = "cl_rt_log2"

N = nuc.get_N()

#%% Main Loop

for i in range(N):
    print("> " + str(i + 1) + "/" + str(N) + "...")
    out_dir = "RandomPattern/" + dataset_name + "/" + str(i)
    chroms = io.load_list("RandomPattern/data/"+str(i)+"/chromosomes.txt")
    
    print(">> Loading data...")
    dataset = nuc.load_nucleosome_pattern_dataset(out_dir, chroms)
    pattern = dataset[0]
    minus_vals, plus_vals, selected_minus, selected_plus = dataset[1:]
    
    pooled_minus_vals = np.zeros((0, np.shape(pattern)[1]))
    pooled_plus_vals = np.zeros((0, np.shape(pattern)[1]))
    pooled_minus_dipys = np.zeros((0, np.shape(pattern)[1]), dtype = bool)
    pooled_plus_dipys = np.zeros((0, np.shape(pattern)[1]), dtype = bool)
    
    for chrom in chroms:
        pooled_minus_vals = np.vstack(
            (pooled_minus_vals, minus_vals[chrom])
        )
        pooled_plus_vals = np.vstack(
            (pooled_plus_vals, plus_vals[chrom])
        )
        pooled_minus_dipys = np.vstack(
            (pooled_minus_dipys, selected_minus[chrom].astype(bool))
        )
        pooled_plus_dipys = np.vstack(
            (pooled_plus_dipys, selected_plus[chrom].astype(bool))
        )
    
    print(">> Overall MLE...")
    overall_results = nuc.pooled_MLE(
        [np.mean(pattern[1,:]), np.mean(np.sqrt(pattern[2,:]))],
        1,
        pooled_minus_vals,
        pooled_plus_vals,
        pooled_minus_dipys,
        pooled_plus_dipys,
        flat_mean,
        [],
        std_relation
    )
    
    overall_MLE = {
        "a": overall_results[0],
        "k": overall_results[1]
    }
    
    print(">> Hyperparameters...")
    hyper_MLE = nuc.hyperparameter_MLE(
        overall_MLE["a"],
        minus_vals,
        plus_vals,
        selected_minus,
        selected_plus,
        flat_mean,
        [],
        std_relation,
        [overall_MLE["k"]]
    )
    
    print(">> Outputting...")
    nuc.output_MLEs(out_dir+"/MLE_values.json", overall_MLE, hyper_MLE)

#%% End

print("===== 041_random_cl_rt_log2_MLE.py - Exiting Properly =====")