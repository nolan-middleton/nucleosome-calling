#############################################################################
# This file will find sequences pertaining to the nucleosomes identied in   #
# mice by Voong et al.                                                      #
#############################################################################

print("===== 0f4_mouse_nucleosome_sequences.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

#%% Load Data

print("> Loading data...")
U = np.loadtxt(
    "../Data/unique.map_95pc.txt",
    dtype = str,
    usecols = [0,1]
)
unique = {}
for chrom in np.unique(U[:,0]):
    unique["chr" + str(chrom)] = U[U[:,0] == chrom,1].astype(int) - 1
del unique["chr22"]
unique["chrX"] = unique["chr20"]
unique["chrY"] = unique["chr21"]
del unique["chr20"]
del unique["chr21"]

#%% Get Sequences in All Nucleosomes

D = "mouse_sequence_analysis"
if (not os.path.isdir(D)):
    os.mkdir(D)

N = 0
AT = 0
AT_minus3 = 0
n = 0
for chrom in unique:
    print("> " + chrom + "...")
    
    print(">> Loading in genome...")
    genome = []
    with open("../Data/mm9.fa") as file:
        foundIt = False
        for line in file:
            if (not foundIt):
                if (line == ">" + chrom + "\n"):
                    foundIt = True
            elif (line[0] != ">"):
                genome += list(line[:-1].upper())
            else:
                break
    genome = np.array(genome)
    
    print(">> Getting sequence...")
    positions = np.zeros((len(unique[chrom]), 147), int)
    for i in range(len(unique[chrom])):
        positions[i,:] = np.arange(
            unique[chrom][i] - 73,
            unique[chrom][i] + 74
        )
    seqs = genome[positions]
    seqs = np.column_stack((unique[chrom], seqs))
    
    print(">> Calculations...")
    N += len(genome)
    AT += np.sum((genome == "A") | (genome == "T"))
    AT_minus3 += np.sum((seqs[:,71] == "A") | (seqs[:,77] == "T")) # col0=pos
    n += len(unique[chrom])
    
    print(">> Saving...")
    np.savetxt(D + "/" + chrom + ".tsv", seqs, fmt = "%s", delimiter = "\t")

io.output_list(list(unique.keys()), D + "/chromosomes.txt")
print(
    "> AT Frequency: "+str(AT/N)+", A/T at -3/+3 Frequency: "+str(AT_minus3/n)
)
with open("mouse_AT3_data.txt", "w") as file:
    file.write("AT: " + str(AT) + "\n")
    file.write("N: " + str(N) + "\n")
    file.write("AT_pm3: " + str(AT_minus3) + "\n")
    file.write("n_nucs: " + str(n) + "\n")

#%% End

print("===== 0f4_mouse_nucleosome_sequences.py - Exiting Properly =====")