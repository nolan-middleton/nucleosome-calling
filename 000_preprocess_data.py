#############################################################################
# This file will handle the preprocessing of data for the nucleosome        #
# calling analysis.                                                         #
#############################################################################

print("===== 000_preprocess_data.py - Starting =====")

#%% Setup

# Imports
import numpy as np
import DataIO as io
import Analysis as an

#%% Preprocessing Function

# Preprocess Data
print("> Data preprocessing...")
minus_data, plus_data, blocks, genome = io.preprocess_data(
    [
        "../Data/WT_0h_allreps_dipy_inbetween_bk_minus.wig",
        "../Data/BY1uvRT_inbetween_norm_minus.wig",
        "../Data/WT_nk_bothreps_dipy_inbetween_bk_norm_minus.wig"
    ],
    [
        "../Data/WT_0h_allreps_dipy_inbetween_bk_plus.wig",
        "../Data/BY1uvRT_inbetween_norm_plus.wig",
        "../Data/WT_nk_bothreps_dipy_inbetween_bk_norm_plus.wig"
    ],
    "../Data/sacCer3.fa",
    "PreprocessedData",
    ["cl", "rt", "ice"],
    2,
    3,
    2000,
    500,
    "floor",
    10,
    "../Data/BY4741_deletions.txt"
)

#%% Unfloored Datasets

# Get unfloored dataset to find differences
print("> Getting unfloored data...")
dinucleotides = io.get_dinucleotides(genome)
plus_dipys = io.get_dipyrimidine_locations(dinucleotides, "+")
minus_dipys = io.get_dipyrimidine_locations(dinucleotides, "-")
unfloored_minus, unfloored_plus = io.make_parallel_dataset(
    [
        "../Data/WT_0h_allreps_dipy_inbetween_bk_minus.wig",
        "../Data/BY1uvRT_inbetween_norm_minus.wig",
        "../Data/WT_nk_bothreps_dipy_inbetween_bk_norm_minus.wig"
    ],
    [
        "../Data/WT_0h_allreps_dipy_inbetween_bk_plus.wig",
        "../Data/BY1uvRT_inbetween_norm_plus.wig",
        "../Data/WT_nk_bothreps_dipy_inbetween_bk_norm_plus.wig"
    ],
    ["unfloored_cl", "unfloored_rt", "unfloored_ice"],
    minus_dipys,
    plus_dipys,
    blocks,
    "PreprocessedData"
)

#%% Comparisons

# Get differences and log2 fold ratios for all three datasets
print("> Comparisons...")
prefixes = ["rt", "ice"]

log2_plus = []
log2_minus = []
diff_minus = []
diff_plus = []
abs_diff_minus = []
abs_diff_plus = []
for i in range(len(prefixes)):
    print(">> Comparing cl vs. " + prefixes[i] + "...")
    
    # Log2 fold ratios
    thisLog2Minus, thisLog2Plus = an.make_comparison(
        minus_data[0],
        plus_data[0],
        minus_data[i + 1],
        plus_data[i + 1],
        "cl_" + prefixes[i] + "_log2",
        "PreprocessedData",
        mode = "log2"
    )
    log2_minus.append(thisLog2Minus)
    log2_plus.append(thisLog2Plus)
    
    # Differences
    thisDiffMinus, thisDiffPlus = an.make_comparison(
        unfloored_minus[0],
        unfloored_plus[0],
        unfloored_minus[i + 1],
        unfloored_plus[i + 1],
        "cl_" + prefixes[i] + "_diff",
        "PreprocessedData",
        mode = "diff"
    )
    diff_minus.append(thisDiffMinus)
    diff_plus.append(thisDiffPlus)
    
    # Absolute differences
    thisAbsDiffMinus, thisAbsDiffPlus = an.make_transformed_dataset(
        thisDiffMinus,
        thisDiffPlus,
        np.abs,
        "cl_" + prefixes[i] + "_abs_diff",
        "PreprocessedData"
    )
    abs_diff_minus.append(thisAbsDiffMinus)
    abs_diff_plus.append(thisAbsDiffPlus)

#%% Strongly-Positioned Nucleosome

# Get nucleosome positions
print("> Nucleosomes...")
nucleosomes = io.read_nucleosome_positions(
    "../Data/widom_sac3_nucscore5_dyads.txt"
)
regions = an.column_join(an.column_split(nucleosomes)[0:2])
badNucs = io.get_bad_intervals(regions, blocks)
for chrom in badNucs:
    nucleosomes[chrom] = nucleosomes[chrom][np.logical_not(badNucs[chrom]),:]

io.output_tsv(nucleosomes, "PreprocessedData/nucleosomes", chroms = True)

table = np.zeros((0,4), dtype = "<U32")
for chrom in nucleosomes:
    table = np.vstack(
        (
            table,
            np.column_stack(
                (
                    np.repeat(chrom, np.shape(nucleosomes[chrom])[0]),
                    nucleosomes[chrom].astype(int)
                )
            )
        )
    )
io.output_wig_interval_from_table(
    table,
    "PreprocessedData/nucleosomes.wig",
    highlight = True
)

#%% All Nucleosome Positions

# All nucleosomes
print("> All nucleosomes...")
table=np.loadtxt("../Data/Brogaard_UniqueDyads.tsv",delimiter="\t",dtype=str)
all_nucleosomes = {}
regions = {}
for chrom in np.unique(table[:,0]):
    all_nucleosomes[str(chrom)] = table[table[:,0] == chrom,1:].astype(float)
    
    all_nucleosomes[str(chrom)] = np.column_stack(
        (
            all_nucleosomes[str(chrom)][:,0] - 73,
            all_nucleosomes[str(chrom)][:,0] + 73,
            all_nucleosomes[str(chrom)][:,0]
        )
    )
    regions[str(chrom)] = np.column_stack(
        (
            all_nucleosomes[str(chrom)][:,0] - 73,
            all_nucleosomes[str(chrom)][:,0] + 73
        )
    )

badNucs = io.get_bad_intervals(regions, blocks)
for chrom in badNucs:
    all_nucleosomes[chrom] = all_nucleosomes[chrom][
        np.logical_not(badNucs[chrom]),
        :
    ]

for chrom in all_nucleosomes:
    all_nucleosomes[chrom] = all_nucleosomes[chrom][
        np.logical_and(
            all_nucleosomes[chrom][:,0] >= 0,
            all_nucleosomes[chrom][:,1] < len(genome[chrom])
        ),
        :
    ]

io.output_tsv(all_nucleosomes,"PreprocessedData/all_nucleosomes",chroms=True)

table = np.zeros((0,4), dtype = "<U32")
for chrom in all_nucleosomes:
    table = np.vstack(
        (
            table,
            np.column_stack(
                (
                    np.repeat(chrom, np.shape(all_nucleosomes[chrom])[0]),
                    all_nucleosomes[chrom].astype(int)
                )
            )
        )
    )
io.output_wig_interval_from_table(
    table,
    "PreprocessedData/all_nucleosomes.wig",
    highlight = True
)

#%% Moderately Well-Placed Nucleosomes

# Decent nucleosomes
print("> Decent nucleosomes...")
decent_nucleosomes = io.read_nucleosome_positions(
    "../Data/widom_saccer3_dyads.txt"
)
regions = an.column_join(an.column_split(decent_nucleosomes)[0:2])
badNucs = io.get_bad_intervals(regions, blocks)
for chrom in badNucs:
    decent_nucleosomes[chrom] = decent_nucleosomes[chrom][
        np.logical_not(badNucs[chrom]),
        :
    ]

poor_nucleosomes = io.read_nucleosome_positions(
    "../Data/widom_sac3_nucscoreless1_dyads.txt"
)
for chrom in poor_nucleosomes:
    if (chrom in decent_nucleosomes):
        decent_nucleosomes[chrom] = decent_nucleosomes[chrom][
            ~np.isin(
                decent_nucleosomes[chrom][:,2],
                poor_nucleosomes[chrom][:,2]
            ),
            :
        ]

io.output_tsv(
    decent_nucleosomes,
    "PreprocessedData/decent_nucleosomes",
    chroms = True
)

table = np.zeros((0,4), dtype = "<U32")
for chrom in decent_nucleosomes:
    table = np.vstack(
        (
            table,
            np.column_stack(
                (
                    np.repeat(chrom, np.shape(decent_nucleosomes[chrom])[0]),
                    decent_nucleosomes[chrom].astype(int)
                )
            )
        )
    )
io.output_wig_interval_from_table(
    table,
    "PreprocessedData/decent_nucleosomes.wig",
    highlight = True
)

#%% Transcription Start Sites
print("> TSS...")
TSS_table = np.loadtxt(
    "../Data/TSS_Iyer_NAR2014_saccer3_1based.txt",
    dtype = str,
    delimiter = "\t",
    skiprows = 1,
    usecols = [0,1,2,3,4]
)

TSS_table[:,1] = (TSS_table[:,1].astype(int) - 1).astype(str) #0-index!

strands = np.array([name.split("-")[0][-1] for name in TSS_table[:,2]])
strands[strands == "C"] = "-"
strands[strands == "W"] = "+"

TSS_table = np.column_stack((TSS_table, strands))
TSS_dict = an.table_to_dict(TSS_table)

TSS_sites = an.make_type(an.column_split(TSS_dict)[0], int)
TSS_regions = {}
for chrom in TSS_sites:
    TSS_regions[chrom] = np.column_stack(
        (TSS_sites[chrom] - 500, TSS_sites[chrom] + 500)
    )

badTSS = io.get_bad_intervals(TSS_regions, blocks)
for chrom in TSS_dict:
    TSS_dict[chrom] = TSS_dict[chrom][np.logical_not(badTSS[chrom]),:]

io.output_tsv(TSS_dict, "PreprocessedData/TSS", fmt = "%s", chroms = True)

#%% Polyadenylation Sites
print("> PAS...")
PAS_table = np.loadtxt(
    "../Data/PAS_Iyer_NAR2014_saccer3_1based.txt",
    dtype = str,
    delimiter = "\t",
    skiprows = 1,
    usecols = [0,1,2]
)

PAS_table[:,1] = (PAS_table[:,1].astype(int) - 1).astype(str) #0-index!

strands = np.array([name.split("-")[0][-1] for name in PAS_table[:,2]])
strands[strands == "C"] = "-"
strands[strands == "W"] = "+"

PAS_table = np.column_stack((PAS_table, strands))
PAS_dict = an.table_to_dict(PAS_table)

PAS_sites = an.make_type(an.column_split(PAS_dict)[0], int)
PAS_regions = {}
for chrom in PAS_sites:
    PAS_regions[chrom] = np.column_stack(
        (PAS_sites[chrom] - 500, PAS_sites[chrom] + 500)
    )

badPAS = io.get_bad_intervals(PAS_regions, blocks)
for chrom in PAS_dict:
    PAS_dict[chrom] = PAS_dict[chrom][np.logical_not(badPAS[chrom]),:]

io.output_tsv(PAS_dict, "PreprocessedData/PAS", fmt = "%s", chroms = True)

#%% End

print("===== 000_preprocess_data.py - Exiting Correctly =====")