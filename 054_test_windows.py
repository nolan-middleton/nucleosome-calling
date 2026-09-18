#############################################################################
# This file will evaluate all of the models using a 20bp window around all  #
# the 1000 random datapoints.                                               #
#############################################################################

print("===== 054_test_windows.py - Starting =====")

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

random_std_relations = {
    "unfloored_cl": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_rt_log2": lambda MLEs, means : MLEs["k"] + np.zeros(len(means)),
    "cl_ice_log2": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_rt_diff": lambda MLEs, means : MLEs["k"] + np.zeros(len(means)),
    "cl_ice_diff": lambda MLEs, means : MLEs["k"] + np.zeros(len(means)),
    "cl_rt_abs_diff": lambda MLEs, means : MLEs["k"]*means + MLEs["l"],
    "cl_ice_abs_diff": lambda MLEs, means : MLEs["k"]*means + MLEs["l"]
}

datasets = [
    "unfloored_cl",
    "cl_rt_log2",
    "cl_ice_log2",
    "cl_rt_diff",
    "cl_ice_diff",
    "cl_rt_abs_diff",
    "cl_ice_abs_diff"
]

# The window size, on EITHER SIDE of the actual dyad
W = 10

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

N = nuc.get_N()

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
    random_std_relation = random_std_relations[dataset]
    
    # Load data
    print(">> Loading data...")
    pattern = np.loadtxt(
        "NucleosomePattern/" + dataset + "/damagePattern.tsv",
        delimiter = "\t"
    )
    offsets = pattern[0,:]
    
    flat_MLE = io.load_json(
        "NucleosomePattern/" + dataset + "/models/flat/MLE_values.json"
    )
    linear_MLE = io.load_json(
        "NucleosomePattern/" + dataset + "/models/linear/MLE_values.json"
    )
    quadratic_MLE = io.load_json(
        "NucleosomePattern/" + dataset + "/models/quadratic/MLE_values.json"
    )
    
    random_MLEs = [
        io.load_json("RandomPattern/"+dataset+"/"+str(i)+"/MLE_values.json")
        for i in range(N)
    ]
    
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
            [name + str(i) + "\t" for i in range(len(offsets) + 2*W)]
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
        "NucleosomePattern/" + dataset + "/windowEvalData.tsv",
        test_data,
        delimiter = "\t",
        header = "".join(colnames)[:-1],
        comments = ""
    )
    
    print(">> Flat integration...")
    flat_integrals = np.zeros((M,2*(2*W + 1)))
    flat_bounds = nuc.nucleosome_bounds(
        flat_MLE,
        nuc.nucleosome_flat_model,
        offsets,
        X,
        nucleosome_std_relation
    )
    
    flat_P_int = nuc.prior_int(flat_bounds, flat_MLE)
    
    for w in range(2*W + 1):
        print(">>> Offset: " + str(w - W) + "...")
        for i in range(M):
            if (i % 250 == 0):
                print("=> " + str(i) + "/" + str(M) + "...")
            flat_integrals[i,2*w],flat_integrals[i,1+2*w]=an.scaled_integral(
                nuc.nucleosome_integrand_objective,
                nuc.nucleosome_integrand,
                [
                    values_minus[i,:][w:(w+len(offsets))],
                    values_plus[i,:][w:(w+len(offsets))],
                    dipys_minus[i,:][w:(w+len(offsets))],
                    dipys_plus[i,:][w:(w+len(offsets))],
                    offsets,
                    X,
                    flat_MLE,
                    nuc.nucleosome_flat_model,
                    nucleosome_std_relation
                ],
                flat_bounds
            )
            
            flat_integrals[i,2*w] /= flat_P_int[0]
            flat_integrals[i,1+2*w] -= flat_P_int[1]
    
    print(">>> Outputting...")
    np.savetxt(
        "NucleosomePattern/"+dataset+"/models/flat/windowEvalInts.tsv",
        flat_integrals,
        delimiter = "\t",
        header = save_file_header,
        comments = ""
    )
    io.output_list(
        flat_P_int,
        "NucleosomePattern/" + dataset + "/models/flat/windowEvalPInt.txt"
    )
    
    print(">> Linear integration...")
    linear_integrals = np.zeros((M,2*(2*W + 1)))
    linear_bounds = nuc.nucleosome_bounds(
        linear_MLE,
        nuc.nucleosome_linear_model,
        offsets,
        X,
        nucleosome_std_relation
    )
    
    linear_P_int = nuc.prior_int(linear_bounds, linear_MLE)
    
    for w in range(2*W + 1):
        print(">>> Offset: " + str(w - W) + "...")
        for i in range(M):
            if (i % 250 == 0):
                print("=> " + str(i) + "/" + str(M) + "...")
            intResults = an.scaled_integral(
                nuc.nucleosome_integrand_objective,
                nuc.nucleosome_integrand,
                [
                    values_minus[i,:][w:(w+len(offsets))],
                    values_plus[i,:][w:(w+len(offsets))],
                    dipys_minus[i,:][w:(w+len(offsets))],
                    dipys_plus[i,:][w:(w+len(offsets))],
                    offsets,
                    X,
                    linear_MLE,
                    nuc.nucleosome_linear_model,
                    nucleosome_std_relation
                ],
                linear_bounds
            )
            linear_integrals[i,2*w], linear_integrals[i,1+2*w] = intResults
            
            linear_integrals[i,2*w] /= linear_P_int[0]
            linear_integrals[i,1+2*w] -= linear_P_int[1]
    
    print(">>> Outputting...")
    np.savetxt(
        "NucleosomePattern/"+dataset+"/models/linear/windowEvalInts.tsv",
        linear_integrals,
        delimiter = "\t",
        header = save_file_header,
        comments = ""
    )
    io.output_list(
        linear_P_int,
        "NucleosomePattern/" + dataset + "/models/linear/windowEvalPInt.txt"
    )
    
    print(">> Quadratic integration...")
    quadratic_integrals = np.zeros((M,2*(2*W + 1)))
    quadratic_bounds = nuc.nucleosome_bounds(
        quadratic_MLE,
        nuc.nucleosome_quadratic_model,
        offsets,
        X,
        nucleosome_std_relation
    )
    
    quadratic_P_int = nuc.prior_int(quadratic_bounds, quadratic_MLE)
    
    for w in range(2*W + 1):
        print(">>> Offset: " + str(w - W) + "...")
        for i in range(M):
            if (i % 250 == 0):
                print("=> " + str(i) + "/" + str(M) + "...")
            intResults = an.scaled_integral(
                nuc.nucleosome_integrand_objective,
                nuc.nucleosome_integrand,
                [
                    values_minus[i,:][w:(w+len(offsets))],
                    values_plus[i,:][w:(w+len(offsets))],
                    dipys_minus[i,:][w:(w+len(offsets))],
                    dipys_plus[i,:][w:(w+len(offsets))],
                    offsets,
                    X,
                    quadratic_MLE,
                    nuc.nucleosome_quadratic_model,
                    nucleosome_std_relation
                ],
                quadratic_bounds
            )
            quadratic_integrals[i,2*w],quadratic_integrals[i,1+2*w]=intResults
            
            quadratic_integrals[i,2*w] /= quadratic_P_int[0]
            quadratic_integrals[i,1+2*w] -= quadratic_P_int[1]
    
    print(">>> Outputting...")
    np.savetxt(
        "NucleosomePattern/"+dataset+"/models/quadratic/windowEvalInts.tsv",
        quadratic_integrals,
        delimiter = "\t",
        header = save_file_header,
        comments = ""
    )
    io.output_list(
        quadratic_P_int,
        "NucleosomePattern/" + dataset + "/models/quadratic/windowEvalPInt.txt"
    )
    
    print(">> Random integration...")
    random_integrals = [np.zeros((M,2*(2*W + 1))) for n in range(N)]
    random_bounds = [
        nuc.random_bounds(
            random_MLEs[n],
            random_std_relation,
        )
        for n in range(N)
    ]
    
    random_P_ints = [
        nuc.prior_int(random_bounds[n], random_MLEs[n])
        for n in range(N)
    ]
    
    for n in range(N):
        print(">>> Random dataset " + str(n + 1) + "/" + str(N) + "...")
        for w in range(2*W + 1):
            print("=> Offset: " + str(w - W) + "...")
            for i in range(M):
                thisPInt = random_P_ints[n]
                thisMarginal, thisLogFactor = an.scaled_integral(
                    nuc.random_integrand_objective,
                    nuc.random_integrand,
                    [
                        values_minus[i,:][w:(w+len(offsets))],
                        values_plus[i,:][w:(w+len(offsets))],
                        dipys_minus[i,:][w:(w+len(offsets))],
                        dipys_plus[i,:][w:(w+len(offsets))],
                        random_MLEs[n],
                        random_std_relation
                    ],
                    random_bounds[n]
                )
                
                random_integrals[n][i,2*w] = thisMarginal / thisPInt[0]
                random_integrals[n][i,1+2*w] = thisLogFactor - thisPInt[1]
        
        print(">>> Outputting...")
        np.savetxt(
            "RandomPattern/" + dataset + "/" + str(n) + "/windowEvalInts.tsv",
            random_integrals[n],
            delimiter = "\t",
            header = save_file_header,
            comments = ""
        )
        io.output_list(
            random_P_ints[n],
            "RandomPattern/" + dataset + "/" + str(n) + "/windowEvalPInt.txt"
        )
    
    print(">> Bayes factors...")
    for w in range(2*W + 1):
        fname = "windowEvalScores_" + str(w-W) + ".tsv"
        flat_bayes = np.column_stack(
            tuple(
                [
                    np.log(flat_integrals[:,2*w])+flat_integrals[:,1+2*w]-(
                        np.log(
                            random_integrals[n][:,2*w]
                        ) + random_integrals[n][:,1+2*w]
                    )
                    for n in range(N)
                ]
            )
        )
        np.savetxt(
            "NucleosomePattern/"+dataset+"/models/flat/" + fname,
            flat_bayes,
            delimiter = "\t",
            header = "".join([str(n) + "\t" for n in range(N)])[:-1],
            comments = ""
        )
    
        linear_bayes = np.column_stack(
            tuple(
                [
                    np.log(linear_integrals[:,2*w])+linear_integrals[:,1+2*w]-(
                        np.log(
                            random_integrals[n][:,2*w]
                        ) + random_integrals[n][:,1+2*w]
                    )
                    for n in range(N)
                ]
            )
        )
        np.savetxt(
            "NucleosomePattern/" + dataset + "/models/linear/" + fname,
            linear_bayes,
            delimiter = "\t",
            header = "".join([str(n) + "\t" for n in range(N)])[:-1],
            comments = ""
        )
    
        quadratic_bayes = np.column_stack(
            tuple(
                [
                    np.log(
                        quadratic_integrals[:,2*w]
                    ) + quadratic_integrals[:,1+2*w] - (
                        np.log(
                            random_integrals[n][:,2*w]
                        ) + random_integrals[n][:,1+2*w]
                    )
                    for n in range(N)
                ]
            )
        )
        np.savetxt(
            "NucleosomePattern/" + dataset + "/models/quadratic/" + fname,
            quadratic_bayes,
            delimiter = "\t",
            header = "".join([str(n) + "\t" for n in range(N)])[:-1],
            comments = ""
        )

#%% End

print("===== 054_test_windows.py - Exiting Properly =====")