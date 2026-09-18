#############################################################################
# This file will analyze the phasing of nucleosomes near the +1 nucleosome. #
#############################################################################

print("===== 0c4_plus_one_rotational.py - Starting =====")

#%% Setup

import numpy as np
import DataIO as io
import Analysis as an
import SequenceHandling as seq

# Variables
R = (50, 150) # The number of base pairs to call something a +1 nucleosome

rng = np.random.default_rng(2781349219413979782)

#%% Load data

print("> Loading data...")
directory = "NucleosomePattern/composite_models/composite_ice"

chroms = io.load_list("PreprocessedData/genome/chromosomes.txt")
chroms.remove("chrM")

# Load in test nucleosome placements
print(">> Dyad calls...")
dyads = io.load_tsv(
    directory + "/viterbi_dyad_calls",
    chromosomes=io.load_list(directory+"/viterbi_dyad_calls/chromosomes.txt"),
    dtype = int
)

scores = io.load_tsv(
    directory + "/genome_scores",
    chromosomes=io.load_list(directory+"/genome_scores/chromosomes.txt")
)

rotational_scores = io.load_tsv(
    directory + "/rotational_scores",
    chromosomes = chroms
)

print(">> Genome and blocks...")
genome = io.load_tsv("PreprocessedData/genome/array", chroms, dtype = str)
blocks = io.load_tsv("PreprocessedData/blocks", chroms)
for chrom in blocks:
    blocks[chrom] = blocks[chrom].astype(int)

print(">> TSS...")
TSS = io.load_tsv("PreprocessedData/TSS", chroms, dtype = str)

#%% Getting +1 Nucleosomes

print("> Plus Ones...")
plus_ones = {}
nucs = {}
dinucs = {}
score_table = {}
plus_one_rots = {}
for chrom in TSS:
    n = len(TSS[chrom])
    plus_ones[chrom] = -np.ones(n, dtype = int)
    nucs[chrom] = np.zeros((n, 1151), dtype = int)
    dinucs[chrom] = np.zeros((n, 1150), dtype = "<U2")
    score_table[chrom] = np.zeros((n, 1151))
    plus_one_rots[chrom] = np.zeros(n)
    
    badBlocks = blocks[chrom][blocks[chrom][:,2] == 1,0:2]
    for i in range(n):
        thisTSS = int(TSS[chrom][i,0])
        thisStrand = TSS[chrom][i,4]
        if (thisStrand == "+"):
            searchRange = (thisTSS - R[0], thisTSS + R[1] + 1)
        else:
            searchRange = (thisTSS - R[1], thisTSS + R[0] + 1)
        
        validNucs = dyads[chrom][
            (dyads[chrom] >= searchRange[0]) & (dyads[chrom] < searchRange[1])
        ]
        
        if (len(validNucs) > 0):
            thisPlusOne = validNucs[np.argmin(abs(validNucs - thisTSS))]
            plus_one_rots[chrom][i] = rotational_scores[chrom][
                rotational_scores[chrom][:,0] == thisPlusOne
            ][0,1]
            
            if (thisStrand == "+"):
                dataRange = (thisPlusOne - 500, thisPlusOne + 651)
            else:
                dataRange = (thisPlusOne - 650, thisPlusOne + 501)
            
            if (dataRange[0]>=73) and (dataRange[1]<=len(genome[chrom])-73):
                isBad = False
                for j in range(np.shape(badBlocks)[0]):
                    if (dataRange[1] > badBlocks[j,0] - 73) and \
                       (dataRange[0] <= badBlocks[j,1] + 73):
                        isBad = True
                        break
                
                if (not isBad):
                    plus_ones[chrom][i] = thisPlusOne
                    
                    if (thisStrand == "+"):
                        indices = dyads[chrom] - thisPlusOne + 500
                    else:
                        indices = thisPlusOne - dyads[chrom] + 500
                    nucs[chrom][
                        i,
                        indices[(indices >= 0) & (indices < 1151)]
                    ] = 1
                    
                    if (thisStrand == "+"):
                        score_table[chrom][i,:] = scores[chrom][
                            (scores[chrom][:,0] >= thisPlusOne - 500) & \
                            (scores[chrom][:,0] < thisPlusOne + 651),
                            1
                        ]
                    else:
                        score_table[chrom][i,:] = scores[chrom][
                            (scores[chrom][:,0] >= thisPlusOne - 650) & \
                            (scores[chrom][:,0] < thisPlusOne + 501),
                            1
                        ][::-1]
                    
                    if (thisStrand == "+"):
                        I = np.arange(thisPlusOne-500+0.5, thisPlusOne+650)
                        dinucs[chrom][i,:] = genome[chrom][
                            (I - 0.5).astype(int)
                        ] + genome[chrom][
                            (I + 0.5).astype(int)
                        ]
                    else:
                        I = np.arange(thisPlusOne-650+0.5, thisPlusOne+500)
                        dinucs[chrom][i,:] = seq.reverse_complements(
                            genome[chrom][(I - 0.5).astype(int)] + \
                            genome[chrom][(I + 0.5).astype(int)]
                        )[::-1]

io.output_tsv(plus_ones, directory + "/plusOnes", "%d")
io.output_tsv(plus_one_rots, directory + "/plusOnes_rotationalScores")
io.output_tsv(score_table, directory + "/plusOnes_alignedScores")
io.output_tsv(nucs, directory + "/plusOnes_alignedNucs", "%d")
io.output_tsv(dinucs, directory + "/plusOnes_alignedDinucs", "%s")

#%% Random Nucleosomes

print("> Randoms...")
total = 0
for chrom in plus_ones:
    total += np.shape(plus_ones[chrom])[0]

chromTotals=np.unique(rng.integers(0,len(chroms),total),return_counts=True)[1]

mock = {}
for i in range(len(chroms)):
    chrom = chroms[i]
    mock[chrom] = []
    for k in range(chromTotals[i]):
        valid = False
        while (not valid):
            valid = True
            thisDyad = dyads[chrom][rng.integers(len(dyads[chrom]))]
            if (thisDyad in plus_ones[chrom]):
                valid = False
            elif (thisDyad in mock[chrom]):
                valid = False
        mock[chrom].append(thisDyad)
    mock[chrom] = np.array(mock[chrom], dtype = int)

mock_nucs = {}
mock_dinucs = {}
mock_score_table = {}
mock_rots = {}
for chrom in mock:
    n = len(mock[chrom])
    mock_nucs[chrom] = np.zeros((n, 1151), dtype = int)
    mock_dinucs[chrom] = np.zeros((n, 1150), dtype = "<U2")
    mock_score_table[chrom] = np.zeros((n, 1151))
    mock_rots[chrom] = np.zeros(n)
    
    badBlocks = blocks[chrom][blocks[chrom][:,2] == 1,0:2]
    for i in range(n):
        thisStrand = ["+", "-"][rng.integers(2)]
        
        thisOne = mock[chrom][i]
        
        if (thisStrand == "+"):
            dataRange = (thisOne - 500, thisOne + 651)
        else:
            dataRange = (thisOne - 650, thisOne + 501)
        
        if (dataRange[0] >= 73) and (dataRange[1] <= len(genome[chrom]) - 73):
            isBad = False
            for j in range(np.shape(badBlocks)[0]):
                if (dataRange[1] > badBlocks[j,0] - 73) and \
                   (dataRange[0] <= badBlocks[j,1] + 73):
                    isBad = True
                    break
            
            if (not isBad):
                if (thisStrand == "+"):
                    indices = dyads[chrom] - thisOne + 500
                else:
                    indices = thisOne - dyads[chrom] + 500
                
                mock_rots[chrom][i] = rotational_scores[chrom][
                    rotational_scores[chrom][:,0] == thisOne
                ][0,1]
                
                mock_nucs[chrom][
                    i,
                    indices[(indices >= 0) & (indices < 1151)]
                ] = 1
                
                if (thisStrand == "+"):
                    mock_score_table[chrom][i,:] = scores[chrom][
                        (scores[chrom][:,0] >= thisOne - 500) & \
                        (scores[chrom][:,0] < thisOne + 651),
                        1
                    ]
                else:
                    mock_score_table[chrom][i,:] = scores[chrom][
                        (scores[chrom][:,0] >= thisOne - 650) & \
                        (scores[chrom][:,0] < thisOne + 501),
                        1
                    ][::-1]
                
                if (thisStrand == "+"):
                    I = np.arange(thisOne-500+0.5, thisOne+650)
                    mock_dinucs[chrom][i,:] = genome[chrom][
                        (I - 0.5).astype(int)
                    ] + genome[chrom][
                        (I + 0.5).astype(int)
                    ]
                else:
                    I = np.arange(thisOne-650+0.5, thisOne+500)
                    mock_dinucs[chrom][i,:] = seq.reverse_complements(
                        genome[chrom][(I - 0.5).astype(int)] + \
                        genome[chrom][(I + 0.5).astype(int)]
                    )[::-1]
            else:
                mock[chrom][i] = -1
        else:
            mock[chrom][i] = -1

io.output_tsv(mock, directory + "/mock_plusOnes", "%d")
io.output_tsv(mock_rots, directory + "/mock_plusOnes_rotationalScores")
io.output_tsv(mock_score_table, directory + "/mock_plusOnes_alignedScores")
io.output_tsv(mock_nucs, directory + "/mock_plusOnes_alignedNucs", "%d")
io.output_tsv(mock_dinucs, directory + "/mock_plusOnes_alignedDinucs", "%s")

#%% End

print("===== 0c4_plus_one_rotational.py - Exiting Properly =====")