#############################################################################
# This file will find the maximum likelihood estimates for the nucleosome   #
# cl_rt_diff model under a variety of different underlying likelihoods.     #
#############################################################################

print("===== 033_nucleosome_cl_rt_diff_MLE.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import Analysis as an
import Nucleosomes as nuc
import DataIO as io
import os as os

# Defs
def F(x):
    # The x*sigmoid(x)
    return x*np.exp(x)/(1+np.exp(x))
def std_relation(A, means):
    '''
    A = [k, l, m, n]
    
    sigma = k*F(n*(means-m)) + l, where F(x)=sigmoid(x)*x
    '''
    return A[0]*F(A[3]*(means - A[2])) + A[1]

out_dir = "NucleosomePattern/cl_rt_diff/models"
if (not os.path.isdir(out_dir)):
    os.mkdir(out_dir)

#%% Load Data

print("> Loading data...")
nucChroms = io.load_list("PreprocessedData/nucleosomes/chromosomes.txt")

data = nuc.load_nucleosome_pattern_dataset(
    "NucleosomePattern/cl_rt_diff",
    nucChroms
)

pattern = data[0]
minus_vals, plus_vals, minus_selected, plus_selected = data[1:]

offsets = pattern[0,:]

X = np.arange(-len(offsets)/2, len(offsets)/2) + 0.5

pooled_minus_vals = np.zeros((0, np.shape(pattern)[1]))
pooled_plus_vals = np.zeros((0, np.shape(pattern)[1]))
pooled_minus_selected = np.zeros((0, np.shape(pattern)[1]), dtype = bool)
pooled_plus_selected = np.zeros((0, np.shape(pattern)[1]), dtype = bool)

for chrom in nucChroms:
    pooled_minus_vals = np.vstack((pooled_minus_vals, minus_vals[chrom]))
    pooled_plus_vals = np.vstack((pooled_plus_vals, plus_vals[chrom]))
    
    pooled_minus_selected = np.vstack(
        (pooled_minus_selected, minus_selected[chrom])
    )
    pooled_plus_selected = np.vstack(
        (pooled_plus_selected, plus_selected[chrom])
    )

#%% Flat MLEs

print("> Flat MLEs...")
flat_results = nuc.pooled_MLE(
    [5, 0.5, 60.67, 10.07, 1.6],
    1,
    pooled_minus_vals,
    pooled_plus_vals,
    pooled_minus_selected,
    pooled_plus_selected,
    nuc.nucleosome_flat_mean,
    [offsets, X],
    std_relation,
    [],
    [{"type": "ineq", "fun": lambda A : A[1]}],
    {"maxiter": 5000}
)

flat_MLE = {
    "a": flat_results[0],
    "k": flat_results[1],
    "l": flat_results[2],
    "m": flat_results[3],
    "n": flat_results[4]
}

#%% Linear MLEs

print("> Linear MLEs...")
linear_results = nuc.pooled_MLE(
    [5, 0, 0.5, 60.67, 10.07, 1.6],
    2,
    pooled_minus_vals,
    pooled_plus_vals,
    pooled_minus_selected,
    pooled_plus_selected,
    nuc.nucleosome_linear_mean,
    [offsets, X],
    std_relation,
    [],
    [{"type": "ineq", "fun": lambda A : A[2]}],
    {"maxiter": 5000}
)

linear_MLE = {
    "a": linear_results[0],
    "b": linear_results[1],
    "k": linear_results[2],
    "l": linear_results[3],
    "m": linear_results[4],
    "n": linear_results[5]
}

#%% Quadratic MLEs

print("> Quadratic MLEs...")

# The optimal p parameter here looks really small, so we'll have to re-scale it
def rescaled_nucleosome_quadratic_mean(A):
    return nuc.nucleosome_quadratic_mean([A[0]] + [A[1] / 1e4] + A[2:])

quadratic_results = nuc.pooled_MLE(
    [6.6, -8.5, 20.5, 0.39, 60.67, 10.07, 2],
    3,
    pooled_minus_vals,
    pooled_plus_vals,
    pooled_minus_selected,
    pooled_plus_selected,
    rescaled_nucleosome_quadratic_mean,
    [offsets, X],
    std_relation,
    [],
    [{"type": "ineq", "fun": lambda A : A[3]}],
    {"maxiter": 5000}
)

quadratic_MLE = {
    "a": quadratic_results[0],
    "p": quadratic_results[1] / 1e4,
    "h": quadratic_results[2],
    "k": quadratic_results[3],
    "l": quadratic_results[4],
    "m": quadratic_results[5],
    "n": quadratic_results[6]
}

#%% Hyperparameters

print("> Flat hyperparameters...")
flat_hyper = nuc.hyperparameter_MLE(
    flat_MLE["a"],
    minus_vals,
    plus_vals,
    minus_selected,
    plus_selected,
    nuc.nucleosome_flat_mean,
    [offsets, X],
    std_relation,
    [flat_MLE["k"], flat_MLE["l"], flat_MLE["m"], flat_MLE["n"]]
)

print("> Linear hyperparameters...")
linear_hyper = nuc.hyperparameter_MLE(
    linear_MLE["a"],
    minus_vals,
    plus_vals,
    minus_selected,
    plus_selected,
    nuc.nucleosome_linear_mean,
    [linear_MLE["b"], offsets, X],
    std_relation,
    [linear_MLE["k"], linear_MLE["l"], linear_MLE["m"], flat_MLE["n"]]
)

print("> Quadratic hyperparameters...")
quadratic_hyper = nuc.hyperparameter_MLE(
    quadratic_MLE["a"],
    minus_vals,
    plus_vals,
    minus_selected,
    plus_selected,
    nuc.nucleosome_quadratic_mean,
    [quadratic_MLE["p"], quadratic_MLE["h"], offsets, X],
    std_relation,
    [quadratic_MLE["k"], quadratic_MLE["l"], linear_MLE["m"], flat_MLE["n"]]
)

#%% Outputs

print("> Outputting...")
# Flat results
flat_out_dir = out_dir + "/flat"
if (not os.path.isdir(flat_out_dir)):
    os.mkdir(flat_out_dir)

nuc.output_MLEs(flat_out_dir + "/MLE_values.json", flat_MLE, flat_hyper)

# Linear results
linear_out_dir = out_dir + "/linear"
if (not os.path.isdir(linear_out_dir)):
    os.mkdir(linear_out_dir)

nuc.output_MLEs(linear_out_dir + "/MLE_values.json", linear_MLE, linear_hyper)

# Quadratic results
quadratic_out_dir = out_dir + "/quadratic"
if (not os.path.isdir(quadratic_out_dir)):
    os.mkdir(quadratic_out_dir)

nuc.output_MLEs(
    quadratic_out_dir + "/MLE_values.json",
    quadratic_MLE,
    quadratic_hyper
)

#%% End

print("===== 033_nucleosome_cl_rt_diff_MLE.py - Exiting Properly =====")