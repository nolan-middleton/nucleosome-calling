#############################################################################
# This file will handle the preprocessing of data from other datasets for   #
# comparing to our analysis.                                                #
#############################################################################

print("===== 001_preprocess_other_datasets.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Analysis as an
import os as os

# Variables
chromosome_names = [
    "",
    "chrI",
    "chrII",
    "chrIII",
    "chrIV",
    "chrV",
    "chrVI",
    "chrVII",
    "chrVIII",
    "chrIX",
    "chrX",
    "chrXI",
    "chrXII",
    "chrXIII",
    "chrXIV",
    "chrXV",
    "chrXVI"
]

out_dir = "PreprocessedData/OtherData"
if (not os.path.isdir(out_dir)):
    os.mkdir(out_dir)

#%% Load Preprocessed Data

print("> Loading initial data...")
chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)

bad_blocks = {}
for chrom in blocks:
    if (np.sum(blocks[chrom][:,2]) > 0):
        bad_blocks[chrom] = blocks[chrom][blocks[chrom][:,2].astype(bool),0:2]

#%% Weiner Dataset

# The Weiner dataset (MNase)
print("> Weiner et al. nucleosome data...")
data = np.loadtxt(
    "../Data/RandoNucleosomemap_mmc3.txt",
    dtype = str,
    delimiter = "\t",
    skiprows = 1
)
chroms = data[:,1].astype(int)
centers = data[:,2].astype(int) - 1 # 0-index!
scores = data[:,3].astype(float)

weiner_data = {}
for chrom_number in np.unique(chroms):
    chrom = chromosome_names[chrom_number]
    Ctrs = centers[chroms == chrom_number]
    Scrs = scores[chroms == chrom_number]
    weiner_data[chrom] = np.column_stack(
        (Ctrs - 73, Ctrs + 73, Ctrs, Scrs)
    )

overlaps = an.find_interval_overlaps(
    an.make_type(an.column_join(an.column_split(weiner_data)[0:2]), int),
    bad_blocks
)
for chrom in overlaps:
    weiner_data[chrom] = weiner_data[chrom][overlaps[chrom] == 0,:]

io.output_tsv(weiner_data, out_dir + "/weiner", chroms = True)

table = np.zeros((0,4), dtype = "<U32")
for chrom in weiner_data:
    table = np.vstack(
        (
            table,
            np.column_stack(
                (
                    np.repeat(chrom, np.shape(weiner_data[chrom])[0]),
                    weiner_data[chrom][:,0:3].astype(int)
                )
            )
        )
    )
io.output_wig_interval_from_table(
    table,
    out_dir + "/weiner_nucleosomes.wig",
    highlight = True
)

#%% Weiner MNase-ChIP-seq Reads

print("> Weiner et al. MNase-ChIP-seq reads...")
weiner_reads = io.preprocess_supplementary_wigs(
    [
        "../Data/GSM1516608_input1.1_tp1_0.wig",
        "../Data/GSM1516556_input1.2_tp1_0.wig",
        "../Data/GSM1516607_input1.3_tp1_0.wig",
        "../Data/GSM1516577_input2.1_tp1_0.wig",
        "../Data/GSM1516571_input2.2_tp1_0.wig",
        "../Data/GSM1516568_input2.3_tp1_0.wig",
        "../Data/GSM1516559_input3.1_tp1_0.wig",
        "../Data/GSM1516603_input3.2_tp1_0.wig",
        
    ],
    blocks,
    [
        out_dir + "/WeinerReads_1_1",
        out_dir + "/WeinerReads_1_2",
        out_dir + "/WeinerReads_1_3",
        out_dir + "/WeinerReads_2_1",
        out_dir + "/WeinerReads_2_2",
        out_dir + "/WeinerReads_2_3",
        out_dir + "/WeinerReads_3_1",
        out_dir + "/WeinerReads_3_2"
    ]
)

#%% Low-Expression Genes

print("> Gene expression levels...")
gene_expression = np.loadtxt(
    "../Data/GeneExpression_holstege98.tsv",
    delimiter = "\t",
    dtype = str,
    skiprows = 1,
    usecols = [0,3]
)
np.savetxt(
    out_dir + "/gene_expression.tsv",
    gene_expression,
    fmt = "%s",
    delimiter = "\t"
)

#%% Nucleosome Modifications

print("> Nucleosome modifications...")
mods = np.loadtxt(
    "../Data/Nucleosome maps/mmc4_0hr_coverage5to40.tsv",
    delimiter = "\t",
    dtype = str
)[:,1:-2]
mods = mods[:,[-2,-1] + list(range(np.shape(mods)[1] - 2))]

colnames = mods[0,:][2:]
mods = mods[1:,:]
nucleosome_mods = {}
for chrom in np.unique(mods[:,0]):
    nucleosome_mods[str(chrom)] = mods[mods[:,0] == chrom, 1:].astype(float)
    nucleosome_mods[str(chrom)][:,0] -= 1 # 0-index!

if (not os.path.isdir(out_dir + "/nucleosome_modifications")):
    os.mkdir(out_dir + "/nucleosome_modifications")

locs = an.column_split(nucleosome_mods)[0]
for i in range(len(colnames)):
    thisDict = an.column_join([locs, an.column_split(nucleosome_mods)[1 + i]])
    io.output_tsv(
        thisDict,
        out_dir + "/nucleosome_modifications/" + colnames[i],
        chroms = True
    )

#%% End

print("===== 001_preprocess_other_datasets.py - Exiting Properly =====")