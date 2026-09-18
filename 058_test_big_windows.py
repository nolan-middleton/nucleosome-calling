#############################################################################
# This file will evaluate all of the models using a 80bp window around all  #
# the 1000 random datapoints.                                               #
#############################################################################

print("===== 058_test_big_windows.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an

# Defs
def F(x):
    # The x*sigmoid(x)
    return x*np.exp(x)/(1+np.exp(x))

# Variables
datasets = [
    "unfloored_cl",
    "cl_rt_log2",
    "cl_ice_log2",
    "cl_rt_diff",
    "cl_ice_diff",
    "cl_rt_abs_diff",
    "cl_ice_abs_diff"
]

nucleosome_std_relations = {
    "unfloored_cl": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_rt_log2": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_ice_log2": lambda MLEs,means:MLEs["k"]*(means-MLEs["m"])**2+MLEs["l"],
    "cl_rt_diff":
        lambda MLEs,means:MLEs["k"]*F(MLEs["n"]*(means-MLEs["m"]))+MLEs["l"],
    "cl_ice_diff":
        lambda MLEs,means:MLEs["k"]*F(MLEs["n"]*(means-MLEs["m"]))+MLEs["l"],
    "cl_rt_abs_diff": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_ice_abs_diff": lambda MLEs, means : MLEs["k"]*means + MLEs["l"]
}
nucleosome_models = {
    "unfloored_cl": ["linear", nuc.nucleosome_linear_model],
    "cl_rt_log2": ["linear", nuc.nucleosome_linear_model],
    "cl_ice_log2": ["linear", nuc.nucleosome_linear_model],
    "cl_rt_diff": ["linear", nuc.nucleosome_linear_model],
    "cl_ice_diff": ["quadratic", nuc.nucleosome_quadratic_model],
    "cl_rt_abs_diff": ["quadratic", nuc.nucleosome_quadratic_model],
    "cl_ice_abs_diff": ["quadratic", nuc.nucleosome_quadratic_model]
}

random_std_relations = {
    "unfloored_cl": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_rt_log2": lambda MLEs, means : MLEs["k"] + np.zeros(len(means)),
    "cl_ice_log2": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_rt_diff": lambda MLEs, means : MLEs["k"] + np.zeros(len(means)),
    "cl_ice_diff": lambda MLEs, means : MLEs["k"] + np.zeros(len(means)),
    "cl_rt_abs_diff": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_ice_abs_diff": lambda MLEs, means : MLEs["k"]*means + MLEs["l"]
}

# The window size, on EITHER SIDE of the actual dyad
W = 80

#%% Load Data

print("> Loading datasets...")
loaded_data = io.load_dataset("PreprocessedData", datasets)

genome, blocks = loaded_data[0:2]

all_data_minus = {}
all_data_plus = {}
for i in range(len(datasets)):
    dataset = datasets[i]
    all_data_minus[dataset],all_data_plus[dataset]=loaded_data[(2+2*i):(4+2*i)]

X = np.arange(-73, 73) + 0.5

test_positions = np.loadtxt(
    "evaluationPositions.tsv",
    delimiter = "\t",
    dtype = str,
    skiprows = 1
)
M = np.shape(test_positions)[0]

#%% Main Loop

save_file_header = "".join(
    np.char.add(
        np.char.add(
            np.char.add(
                "marginal_",
                (np.array(range(2*W + 1)) - W).astype(str)
            ),
            "\t"
        ),
        np.char.add(
            np.char.add(
                "log_factor_",
                (np.array(range(2*W + 1)) - W).astype(str)
            ),
            "\t"
        )
    )
)[:-1]

for dataset in datasets:
    print("> " + dataset + "...")
    data_minus = all_data_minus[dataset]
    data_plus = all_data_plus[dataset]
    
    nucleosome_std_relation = nucleosome_std_relations[dataset]
    model_name = nucleosome_models[dataset][0]
    nucleosome_model = nucleosome_models[dataset][1]
    random_std_relation = random_std_relations[dataset]
    
    # Load data
    print(">> Loading data...")
    pattern = np.loadtxt(
        "NucleosomePattern/" + dataset + "/damagePattern.tsv",
        delimiter = "\t"
    )
    offsets = pattern[0,:]
    
    nucleosome_MLE = io.load_json(
        "NucleosomePattern/"+dataset+"/models/"+model_name+"/MLE_values.json"
    )
    
    random_MLE = io.load_json("RandomPattern/"+dataset+"/0/MLE_values.json")
    
    print(">> Forming data array...")
    values_minus = np.zeros((0, len(offsets) + 2*W))
    values_plus = np.zeros((0, len(offsets) + 2*W))
    dipys_minus = np.zeros((0, len(offsets) + 2*W), dtype = bool)
    dipys_plus = np.zeros((0, len(offsets) + 2*W), dtype = bool)
    
    for i in range(M):
        chrom = test_positions[i,1]
        start = int(test_positions[i,2]) - W
        stop = int(test_positions[i,3]) + W
        
        positions = np.arange(start, stop) + 0.5
        
        values_minus = np.vstack((values_minus, np.zeros(len(offsets)+2*W)))
        dipys_minus = np.vstack(
            (dipys_minus, np.zeros(len(offsets) + 2*W, dtype = bool))
        )
        dummy, positionIndices, dataIndices = np.intersect1d(
            positions,
            data_minus[chrom][:,0],
            return_indices = True
        )
        values_minus[i,:][positionIndices] = data_minus[chrom][dataIndices,1]
        dipys_minus[i,:][positionIndices] = True
        
        values_plus = np.vstack((values_plus, np.zeros(len(offsets)+2*W)))
        dipys_plus = np.vstack(
            (dipys_plus, np.zeros(len(offsets) + 2*W, dtype = bool))
        )
        dummy, positionIndices, dataIndices = np.intersect1d(
            positions,
            data_plus[chrom][:,0],
            return_indices = True
        )
        values_plus[i,:][positionIndices] = data_plus[chrom][dataIndices,1]
        dipys_plus[i,:][positionIndices] = True
    
    test_data = np.column_stack(
        (values_minus, dipys_minus, values_plus, dipys_plus)
    )
    
    print(">>> Outputting...")
    colnames = sum(
        [
            [name + str(i) + "\t" for i in range(len(offsets))]
            for name in [
                "minus_values_",
                "minus_dipys_",
                "plus_values_",
                "plus_dipys_"
            ]
        ],
        start = []
    )
    np.savetxt(
        "NucleosomePattern/" + dataset + "/bigWindowEvalData.tsv",
        test_data,
        delimiter = "\t",
        header = "".join(colnames)[:-1],
        comments = ""
    )
    
    print(">> Nucleosome integration...")
    nucleosome_integrals = np.zeros((M,2*(2*W + 1)))
    nucleosome_bounds = nuc.nucleosome_bounds(
        nucleosome_MLE,
        nucleosome_model,
        offsets,
        X,
        nucleosome_std_relation
    )
    nuc_P_int = nuc.prior_int(nucleosome_bounds, nucleosome_MLE)
    
    for w in range(2*W + 1):
        print(">>> Offset: " + str(w - W) + "...")
        for i in range(M):
            theseIntResults = an.scaled_integral(
                nuc.nucleosome_integrand_objective,
                nuc.nucleosome_integrand,
                [
                    values_minus[i,:][w:(w+len(offsets))],
                    values_plus[i,:][w:(w+len(offsets))],
                    dipys_minus[i,:][w:(w+len(offsets))],
                    dipys_plus[i,:][w:(w+len(offsets))],
                    offsets,
                    X,
                    nucleosome_MLE,
                    nucleosome_model,
                    nucleosome_std_relation
                ],
                nucleosome_bounds
            )
            nucleosome_integrals[i,2*w] = theseIntResults[0] / nuc_P_int[0]
            nucleosome_integrals[i,1+2*w] = theseIntResults[1] - nuc_P_int[1]
    
    print(">>> Outputting...")
    np.savetxt(
        "NucleosomePattern/" + dataset + "/models/bigWindowEvalInts.tsv",
        nucleosome_integrals,
        delimiter = "\t",
        header = save_file_header,
        comments = ""
    )
    io.output_list(
        nuc_P_int,
        "NucleosomePattern/" + dataset + "/models/bigWindowEvalPInt.txt"
    )
    
    print(">> Random integration...")
    random_integrals = np.zeros((M,2*(2*W + 1)))
    random_bounds = nuc.random_bounds(random_MLE, random_std_relation)
    
    rand_P_int = nuc.prior_int(random_bounds, random_MLE)

    for w in range(2*W + 1):
        print("=> Offset: " + str(w - W) + "...")
        for i in range(M):
            if (i % 250 == 0):
                print("==> " + str(i) + "/" + str(M) + "...")
            thisMarginal, thisLogFactor = an.scaled_integral(
                nuc.random_integrand_objective,
                nuc.random_integrand,
                [
                    values_minus[i,:][w:(w+len(offsets))],
                    values_plus[i,:][w:(w+len(offsets))],
                    dipys_minus[i,:][w:(w+len(offsets))],
                    dipys_plus[i,:][w:(w+len(offsets))],
                    random_MLE,
                    random_std_relation
                ],
                random_bounds
            )
            
            random_integrals[i,2*w] = thisMarginal / rand_P_int[0]
            random_integrals[i,1+2*w] = thisLogFactor - rand_P_int[1]
    
    print(">>> Outputting...")
    np.savetxt(
        "RandomPattern/" + dataset + "/bigWindowEvalInts.tsv",
        random_integrals,
        delimiter = "\t",
        header = save_file_header,
        comments = ""
    )
    io.output_list(
        rand_P_int,
        "RandomPattern/" + dataset + "/bigWindowEvalPInt.txt"
    )
    
    print(">> Bayes factors...")
    bayes = np.zeros((M, 0))
    for w in range(2*W + 1):
        bayes = np.column_stack(
            (
                bayes,
                np.log(
                    nucleosome_integrals[:,2*w]
                ) + nucleosome_integrals[:,1+2*w] - np.log(
                    random_integrals[:,2*w]
                ) - random_integrals[:,1+2*w]
            )
        )
    np.savetxt(
        "NucleosomePattern/" + dataset + "/models/bigWindowEvalBayes.tsv",
        bayes,
        delimiter = "\t",
        header = "".join([str(w) + "\t" for w in range(-W, W+1)])[:-1],
        comments = ""
    )

#%% End

print("===== 058_test_big_windows.py - Exiting Properly =====")