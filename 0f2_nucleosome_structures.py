#############################################################################
# This file will find which RCSB PDB nucleosome structures have an A at the #
# -3 position.                                                              #
#############################################################################

print("===== 0f2_nucleosome_structures.py - Starting =====")

#%% Setup

import numpy as np
import DataIO as io

#%% Load data

table = np.loadtxt(
    "../Data/PDB_structures.tsv",
    delimiter = "\t",
    dtype = str,
    skiprows = 1
)

# We only need one strand
U, I = np.unique(table[:,0], return_index = True)
T = np.column_stack((U, table[I,:][:,1:]))

#%% Reading the PDB Files

minus3s = {}
plus3s = {}
for i in range(np.shape(T)[0]):
    print("> " + str(i + 1) + "/" + str(np.shape(T)[0]) + "...")
    data=io.read_mmCIF_file("../Data/Nucleosome Structures/"+T[i,0]+".cif")
    
    seqs = data["pdbx_poly_seq_scheme"]
    indices = []
    for idx in range(len(seqs["asym_id"])):
        if (seqs["pdb_strand_id"][idx] == T[i,1]):
            indices.append(idx)
    
    minus3 = str(int(T[i,2]) - 3)
    plus3 = str(int(T[i,2]) + 3)
    
    done = 0
    for idx in indices:
        if (done == 2):
            break
        elif (seqs["pdb_seq_num"][idx] == minus3):
            minus3s[str(T[i,0])] = seqs["pdb_mon_id"][idx]
            done += 1
        elif (seqs["pdb_seq_num"][idx] == plus3):
            plus3s[str(T[i,0])] = seqs["pdb_mon_id"][idx]
            done += 1

#%% Totalling

total = 0
for key in plus3s:
    if (minus3s[key] == "DA") or (plus3s[key] == "DT"):
        total += 1

print("> Total: " + str(total) + " of " + str(len(plus3s)) + "...")

#%% End

print("===== 0f2_nucleosome_structures.py - Exiting Properly =====")