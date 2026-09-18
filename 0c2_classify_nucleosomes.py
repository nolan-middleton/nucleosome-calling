#############################################################################
# This file will classify the nucleosomes as in genes, in a PolII gene, +1, #
# etc.                                                                      #
#############################################################################

print("===== 0c2_classify_nucleosomes.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Nucleosomes as nuc
import Analysis as an
import os as os

# Variables
plus_one_tolerance = (-30, 100)
nuc_tol = 37
low_expr_threshold = 1
high_expr_threshold = 10

#%% Load Data

print("> Loading data...")
dirs = [
    "NucleosomePattern/cl_rt_log2",
    "NucleosomePattern/cl_ice_log2",
    "NucleosomePattern/cl_ice_abs_diff",
    "NucleosomePattern/composite_models/composite_ice"
]
num_dirs = len(dirs)

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

# Load in test nucleosome placements
print(">> Putative dyad calls...")
dyads = [
    io.load_tsv(
        D + "/viterbi_dyad_calls",
        chromosomes = io.load_list(D + "/viterbi_dyad_calls/chromosomes.txt"),
        dtype = int
    )
    for D in dirs
]

print(">> Genome and blocks...")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)
for chrom in blocks:
    blocks[chrom] = blocks[chrom].astype(int)

print(">> TSS and PAS")
TSS = io.load_tsv("PreprocessedData/TSS", chroms, dtype = str)
PAS = io.load_tsv("PreprocessedData/PAS", chroms, dtype = str)

genes = {}
for chrom in chroms:
    genes[chrom] = TSS[chrom][:,[1,4]]
    
    theseTSS = TSS[chrom][:,0].astype(int)
    thesePAS = -np.ones(len(theseTSS), dtype = int)
    for i in range(len(genes[chrom])):
        thisPAS = PAS[chrom][PAS[chrom][:,1] == genes[chrom][i,0],0]
        if (len(thisPAS) > 0):
            thesePAS[i] = thisPAS[0]
    
    bad = thesePAS == -1
    
    theseCoords = np.column_stack((theseTSS, thesePAS))
    
    genes[chrom] = np.column_stack(
        (genes[chrom], np.min(theseCoords,axis=1), np.max(theseCoords,axis=1))
    )[~bad,:]

print(">> SGD other...")
sgdOther = np.loadtxt("../Data/sgdOther.txt", delimiter="\t", dtype=str)[:,1:]

other = {}
for chrom in chroms:
    other[chrom] = sgdOther[sgdOther[:,0] == chrom, 1:]

tRNA = {}
snoRNA = {}
LTR = {}
ARS = {}
telomeric = {}
mRNA = {}
centromeric = {}

ncRNA = {}

for chrom in other:
    tRNA[chrom] = other[chrom][
        np.strings.find(other[chrom][:,5], "tRNA") >= 0,
        0:2
    ].astype(int) - 1 # 0-index!
    snoRNA[chrom] = other[chrom][
        np.strings.find(other[chrom][:,5], "snoRNA") >= 0,
        0:2
    ].astype(int) - 1 # 0-index!
    LTR[chrom] = other[chrom][
        np.strings.find(other[chrom][:,5], "LTR") >= 0,
        0:2
    ].astype(int) - 1 # 0-index!
    ARS[chrom] = other[chrom][
       np.strings.find(other[chrom][:,5], "ARS") >= 0,
       0:2
    ].astype(int) - 1 # 0-index!
    telomeric[chrom] = other[chrom][
       np.strings.find(other[chrom][:,5], "Telomeric") >= 0,
       0:2
    ].astype(int) - 1 # 0-index!
    mRNA[chrom] = other[chrom][
       np.strings.find(other[chrom][:,5], "SGD, mRNA") >= 0,
       0:2
    ].astype(int) - 1 # 0-index!
    centromeric[chrom] = other[chrom][
       np.strings.find(other[chrom][:,5], "centromere") >= 0,
       0:2
    ].astype(int) - 1 # 0-index!
    
    ncRNA[chrom] = other[chrom][
       np.strings.find(other[chrom][:,5], "SGD, noncoding_exon") >= 0,
       0:2
    ].astype(int) - 1 # 0-index!

print(">> Gene expression levels...")
gene_expr = np.loadtxt(
    "PreprocessedData/OtherData/gene_expression.tsv",
    delimiter = "\t",
    dtype = str
)

print(">> Nucleosome modifications...")
nuc_mods = {}
for item in os.scandir("PreprocessedData/OtherData/nucleosome_modifications"):
    nuc_mods[item.name] = io.load_tsv(
        item.path,
        io.load_list(item.path + "/chromosomes.txt")
    )

#%% Main Loop

colnames = [
    "in_gene",
    "+1",
    "tRNA",
    "snoRNA",
    "LTR",
    "ARS",
    "telomeric",
    "misc_mRNA",
    "centromeric",
    "ncRNA",
    "in_low_expr_gene",
    "in_high_expr_gene"
] + list(nuc_mods.keys())

gene_names, gene_strands, gene_starts, gene_stops = an.column_split(genes)
gene_starts = an.make_type(gene_starts, int)
gene_stops = an.make_type(gene_stops, int)

minus_genes = {}
plus_genes = {}
low_expr_starts = {}
low_expr_stops = {}
high_expr_starts = {}
high_expr_stops = {}
for chrom in gene_strands:
    minus_genes[chrom] = np.column_stack(
        (
            gene_stops[chrom][gene_strands[chrom]=="-"]-plus_one_tolerance[1],
            gene_stops[chrom][gene_strands[chrom]=="-"]-plus_one_tolerance[0]
        )
    )
    plus_genes[chrom] = np.column_stack(
        (
            gene_stops[chrom][gene_strands[chrom]=="+"]+plus_one_tolerance[0],
            gene_stops[chrom][gene_strands[chrom]=="+"]+plus_one_tolerance[1]
        )
    )
    
    theseLow = np.isin(
        gene_names[chrom],
        gene_expr[gene_expr[:,1].astype(float) < low_expr_threshold,0]
    )
    low_expr_starts[chrom] = gene_starts[chrom][theseLow]
    low_expr_stops[chrom] = gene_stops[chrom][theseLow]
    
    theseHigh = np.isin(
        gene_names[chrom],
        gene_expr[gene_expr[:,1].astype(float) > high_expr_threshold,0]
    )
    high_expr_starts[chrom] = gene_starts[chrom][theseHigh]
    high_expr_stops[chrom] = gene_stops[chrom][theseHigh]

for i in range(len(dirs)):
    print("> " + dirs[i] + "...")
    D = dyads[i]
    classes = {}
    
    for chrom in D:
        print(">> " + chrom + "...")
        classes[chrom] = np.zeros((len(D[chrom]), len(colnames)), dtype="<U32")
        
        for j in range(len(D[chrom])):
            dyad = D[chrom][j]
            
            # in_gene
            classes[chrom][j,0] = True in (
                (dyad >= gene_starts[chrom]) & (dyad <= gene_stops[chrom])
            )
            # +1
            classes[chrom][j,1] = (
                True in (
                    ((dyad >= minus_genes[chrom][:,0]) & \
                     (dyad <= minus_genes[chrom][:,1]))
                )
            ) or (
                True in (
                    ((dyad >= plus_genes[chrom][:,0]) & \
                     (dyad <= plus_genes[chrom][:,1]))
                )
            )
            # tRNA
            classes[chrom][j,2] = True in (
                (dyad >= tRNA[chrom][:,0]) & (dyad <= tRNA[chrom][:,1])
            )
            # snoRNA
            classes[chrom][j,3] = True in (
                (dyad >= snoRNA[chrom][:,0]) & (dyad <= snoRNA[chrom][:,1])
            )
            # LTR
            classes[chrom][j,4] = True in (
                (dyad >= LTR[chrom][:,0]) & (dyad <= LTR[chrom][:,1])
            )
            # ARS
            classes[chrom][j,5] = True in (
                (dyad >= ARS[chrom][:,0]) & (dyad <= ARS[chrom][:,1])
            )
            # telomeric
            classes[chrom][j,6] = True in (
                (dyad>=telomeric[chrom][:,0]) & (dyad<=telomeric[chrom][:,1])
            )
            # misc_mRNA
            classes[chrom][j,7] = True in (
                (dyad >= mRNA[chrom][:,0]) & (dyad <= mRNA[chrom][:,1])
            )
            # centromeric
            classes[chrom][j,8] = True in (
                (dyad>=centromeric[chrom][:,0])&(dyad<=centromeric[chrom][:,1])
            )
            # ncRNA
            classes[chrom][j,9] = True in (
                (dyad >= ncRNA[chrom][:,0]) & (dyad <= ncRNA[chrom][:,1])
            )
            # in_low_expr_gene
            classes[chrom][j,10] = True in (
                (dyad>=low_expr_starts[chrom]) & (dyad<=low_expr_stops[chrom])
            )
            # in_high_expr_gene
            classes[chrom][j,11] = True in (
                (dyad>=high_expr_starts[chrom])&(dyad<=high_expr_stops[chrom])
            )
            
            # Nucleosome Modifications
            for k in range(12, len(colnames)):
                mod = colnames[k]
                closest = np.argmin(np.abs(dyad - nuc_mods[mod][chrom][:,0]))
                if (np.abs(dyad - nuc_mods[mod][chrom][closest,0]) < nuc_tol):
                    classes[chrom][j,k] = nuc_mods[mod][chrom][closest,1]
                else:
                    classes[chrom][j,k] = np.nan
    
    io.output_tsv(
        classes,
        dirs[i] + "/nucleosome_classification",
        fmt = "%s",
        chroms = True
    )
    
    master_table = np.vstack(
        tuple(
            [
                np.column_stack(
                    (
                        np.repeat(chrom, len(D[chrom])),
                        D[chrom],
                        classes[chrom]
                    )
                )
                for chrom in D
            ]
        )
    )
    
    np.savetxt(
        dirs[i] + "/nucleosome_classification.tsv",
        master_table,
        fmt = "%s",
        header = "chrom\tloc\t" + "\t".join(colnames),
        comments = "",
        delimiter = "\t"
    )

#%% End

print("===== 0c2_classify_nucleosomes.py - Exiting Properly =====")