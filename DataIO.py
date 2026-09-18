##############################################################################
# This file contains library functions for handling data input and output,   #
# including functions for reading raw input files, for outputting .tsv and   #
# .wig files, for preprocessing data, for normalizing data, and for later    #
# import of preprocessed data.                                               #
##############################################################################

#%% Setup

# Imports
import numpy as np
import os as os
import json as json

#%% Raw Data Inputs

def read_mmCIF_file(filename):
    '''
    This function will read in a .cif file from the RCSB PDB and return a
    dictionary of all the data.
    
    Arguments:
        filename: The name of the file to read.
    '''
    # Load lines
    text = ""
    with open(filename) as file:
        for line in file:
            text += line
    
    data = {}
    thisSection = {}
    thisKey = "<top>"
    inLoop = False
    theseCols = []
    thisIndex = 0
    inMultiline = False
    inQuotes = False
    thisQuoteChar = ""
    thisWord = ""
    key = ""
    for char in text:
        wordBreak = False
        if (char == "#"): # Section break
            data[thisKey] = thisSection
            thisSection = {}
            thisKey = ""
            inLoop = False
            theseCols = []
            thisIndex = 0
            inMultiline = False
            inQuotes = False
            thisQuoteChar == ""
            thisWord = ""
            key = ""
        elif (not inQuotes) and ((char == "'") or (char == "\"")):
            inQuotes = True
            thisQuoteChar = char
            wordBreak = True
        elif (inQuotes and char == thisQuoteChar):
            inQuotes = False
            thisQuoteChar = ""
            wordBreak = True
        elif (not inQuotes and char == ";"):
            inMultiline = not inMultiline
            wordBreak = True
        elif (not inQuotes) and (char.isspace()):
            if (not inMultiline):
                wordBreak = True
        else:
            thisWord += char
        
        if (wordBreak):
            if (thisWord == ""): # To handle repeated whitespaces
                pass
            elif (thisWord == "loop_"): # Loop initialization
                inLoop = True
                theseCols = []
                thisIndex = 0
                inMultiline = False
                thisSection = {}
            elif (inLoop and thisWord[0] == "_"): # Loop data rows
                thisKey, key = thisWord[1:].split(".")
                theseCols.append(key)
                thisSection[key] = []
            elif (inLoop):
                thisSection[theseCols[thisIndex]].append(thisWord)
                thisIndex = (thisIndex + 1) % len(theseCols)
            elif (not inLoop and thisWord[0] == "_"): # key
                thisKey, key = thisWord[1:].split(".")
            elif (not inLoop and key != ""):
                thisSection[key] = thisWord
            elif (not inLoop):
                thisSection = thisWord
            
            thisWord = ""
    return data

def read_raw_wig(filename, overlap = "sum"):
    '''
    This function will reads in a .wig file and spit out its values as a
    dictionary, where each key is a chromosome and each entry is a numpy array
    with three columns: chrom, location, value.
    
    Note this function does NOT do anything to the data, just reads it into a
    dictionary of numpy tables.
    
    Arguments:
        filename: The name of the .wig file to read.
        overlap: A string denoting to handle any overlapping regions in the
            .wig file. Options are: "sum", "max", "min", "median" or "mean.
            Defaults to "sum".
    
    The .wig file must start each chromosome section with "variableStep" or
    "fixedStep" and be tab-delimited.
    '''
    if (overlap not in ["sum", "max", "min", "median", "mean"]): 
        raise Exception("Invalid overlap mode.")
    
    # Read in data
    lines = []
    headers = []
    header_indices = []
    with open(filename) as file:
        for line in file:
            L = line.split("\n")[0]
            lines.append(L)
            if ("variableStep" in line) or ("fixedStep" in line):
                header_indices.append(len(lines) - 1)
                headers.append(L)
    
    # Convert to array format
    dataDict = {}
    for i in range(len(headers)):
        H = headers[i]
        mode = H.split(" ")[0]
        chrom = H.split("chrom=")[1].split(" ")[0]
        span=int(H.split("span=")[1].split(" ")[0]) if ("span" in H) else 1
        start=float(H.split("start=")[1].split(" ")[0]) if("start" in H) else 1
        step=float(H.split("step=")[1].split(" ")[0]) if ("step" in H) else 1
        
        start_i = header_indices[i] + 1
        stop_i = header_indices[i+1] if (i < len(headers) - 1) else len(lines)
        
        data = np.zeros(((stop_i - start_i)*span, 2))
        
        if (mode == "variableStep"):
            thisData = np.zeros((span, 2))
            for j in range(start_i, stop_i):
                D = lines[j].split("\t")
                thisData[:,0] = np.arange(float(D[0]), float(D[0])+span)
                thisData[:,1] = float(D[1])
                data[j-start_i,:] = thisData
        elif (mode == "fixedStep"):
            pos = np.arange(start, start + span)
            thisData = np.zeros(span)
            for j in range(start_i, stop_i):
                thisData[:] = float(lines[j])
                data[(j-start_i):(j-start_i+span),0] = pos
                data[(j-start_i):(j-start_i+span),1] = thisData
                pos += step
        
        if (chrom not in dataDict):
            dataDict[chrom] = []
        
        dataDict[chrom].append(data)
    
    # Handle overlaps
    return_dict = {}
    for chrom in dataDict:
        thisArray = np.vstack(tuple(dataDict[chrom]))
        thisArray = thisArray[np.argsort(thisArray[:,0]),:]
        
        U, counts = np.unique(thisArray[:,0], return_counts = True)
        if (len(U) != np.shape(thisArray)[0]):
            if (overlap == "sum"):
                F = np.sum
            elif (overlap == "max"):
                F = np.max
            elif (overlap == "min"):
                F = np.min
            elif (overlap == "median"):
                F = np.median
            elif (overlap == "mean"):
                F = np.mean
            
            C = counts > 1
            degenerates = np.zeros((np.sum(C), 2))
            degenerates[:,0] = U[np.arange(len(U))[C]]
            for i in range(len(degenerates)):
                degenerates[i,1] = F(
                    thisArray[thisArray[:,0] == degenerates[i,0],1]
                )
                thisArray = thisArray[thisArray[:,0] != degenerates[i,0],:]
            
            thisArray = np.vstack((thisArray, degenerates))
            thisArray = thisArray[np.argsort(thisArray[:,0]),:]
        
        return_dict[chrom] = thisArray
    
    return return_dict
            

def read_raw_genome(filename):
    '''
    This function will read in the genome from a .fasta file and produce the
    genome in string form (not in array form).
    
    Arguments:
        filename: the name of the genome .fasta file to read.
    
    The .fasta headers must contain the chromosome name, nothing else.
    '''
    genome = {}
    
    with open(filename) as file:
        chrom = ""
        data = ""
        for line in file:
            if (line[0] == ">"):
                # New chromosome
                if (chrom != ""):
                    genome[chrom] = data
                    data = ""
                chrom = line[1:].split("\n")[0]
            else:
                data += line.split("\n")[0]
        genome[chrom] = data
    
    return genome

def convert_genome_to_array(str_genome):
    '''
    This function will take in a genome dictionary of strings and spit out a
    genome dictionary in array form.
    
    Arguments:
        str_genome: the genome dictionary of strings.
    '''
    arr_genome = {}
    
    for chrom in str_genome:
        arr_genome[chrom] = np.array(list(str_genome[chrom]))
    
    return arr_genome

def get_dinucleotides(genome):
    '''
    This function will generate dinucleotides for the entire genome.
    
    Arguments:
        genome: The genome dictionary in ARRAY FORM
    '''
    dinucs = {}
    
    for chrom in genome:
        dinucs[chrom] = np.char.add(genome[chrom][:-1], genome[chrom][1:])
    
    return dinucs

def get_dipyrimidine_locations(dinucleotides, strand):
    '''
    This function will get the locations of dipyrimidines in the genome.
    
    Arguments:
        dinucleotides: The dinucleotide dictionary.
        strand: "+" or "-" for plus or minus-stranded data.
    '''
    locs = {}
    if (strand == "+"):
        seqs = ["TT", "TC", "CT", "CC"]
    else:
        seqs = ["AA", "GA", "AG", "GG"]
    for chrom in dinucleotides:
        length = len(dinucleotides[chrom])
        num_seqs = len(seqs)
        truth_table = np.zeros((length, num_seqs)).astype(bool)
        for i in range(num_seqs):
            truth_table[:,i] = dinucleotides[chrom] == seqs[i]
        locs[chrom]=np.array(range(length))[np.sum(truth_table,axis=1)>0]+0.5
    
    return locs

def reformat_wig_data(wigData, dipyrimidine_locs):
    '''
    This function will convert the raw data from a .wig file into a more
    useable format, adding missing values as zeros and zero-indexing the
    locations.
    
    Arguments:
        wigData: The data from the .wig file, should be from by read_raw_wig.
        dipyrimidine_locs: The dictionary of dipyrimidine locations.
    '''
    returnData = {}
    for chrom in dipyrimidine_locs:
        dipys = dipyrimidine_locs[chrom]
        
        if (chrom in wigData):
            locs = wigData[chrom][:,0] - 1 # zero-index locations
            vals = wigData[chrom][:,1]
        else:
            locs = np.array([])
            vals = np.array([])
        
        # Every dipyrimidine will get a row as a location, filling in zeros
        # for missing values
        dummy, oldIndices, newIndices = np.intersect1d(
            locs,
            dipys,
            return_indices = True
        )
        newVals = np.zeros(len(dipys))
        newVals[newIndices] = vals[oldIndices]
        
        returnData[chrom] = np.column_stack((dipys, newVals))
    
    return returnData

def read_nucleosome_positions(filename):
    '''
    This function will read in a .txt file of nucleosome dyad positions and
    return a dictionary of nucleosome intervals and dyad positions.
    
    Arguments:
        filename: The name of the .txt file of nucleosome dyad positions.
    '''
    nucleosomes = {}
    
    with open(filename) as file:
        for line in file:
            chrom, dinuc = line.split(":")
            if (chrom not in nucleosomes):
                nucleosomes[chrom] = np.array([[],[], []]).reshape((0,3))
            # Remember to 0-index by subtracting one...
            dyad = np.array(dinuc.split("-")).astype(float)[0] - 1
            nucleosomes[chrom] = np.vstack(
                (nucleosomes[chrom], [dyad - 73, dyad + 73, dyad])
            )
    
    return nucleosomes

def read_intervals(filename, adjustIndex = True):
    '''
    This function will read in a .txt file of intervals expressed as genomic
    coordinates and return a dictionary of start and stop positions.
    
    Arguments:
        filename: The name of the .txt file of the intervals to read in.
        adjustIndex: If set to True, will subtract one from the given
            coordinates to 0-index them. Defaults to True.
    '''
    intervals = {}
    
    with open(filename) as file:
        for line in file:
            chrom, interval = line.split(":")
            start, stop = interval.split("-")
            start = int(start) - int(adjustIndex)
            stop = int(stop[:-1]) - int(adjustIndex)
            
            if (chrom not in intervals):
                intervals[chrom] = np.array([[],[]]).reshape((0,2))
            intervals[chrom] = np.vstack(
                (intervals[chrom], np.array([start, stop]))
            )
    
    return intervals

def read_fasta(filename):
    '''
    This function will read in a FASTA file of sequences and return a numpy
    array of strings, where column 0 is the sequence name and column 1 is the
    sequence itself. If other data exists in the header, it will show up in
    column 2.
    
    Arguments:
        filename: The name of the file to read in.
    '''
    sequences = []
    names = []
    other_data = []
    
    with open(filename) as file:
        thisSeq = ""
        for line in file:
            if (line[0] == ">"):
                # Header
                thisSplit = line[1:-1].split(" ")
                names.append(thisSplit[0])
                other_data.append("".join(thisSplit[1:]))
                
                if (len(names) > 1):
                    # Don't append the sequence to the first header
                    sequences.append(thisSeq)
                    thisSeq = ""
            else:
                thisSeq += line[:-1]
        sequences.append(thisSeq)
    
    returnTable = np.column_stack((names, sequences, other_data))
    
    if (False in (returnTable[:,2] == "")):
        return returnTable
    else:
        return returnTable[:,:-1]

#%% Processed Data Inputs

def load_list(filename):
    '''
    This function will read in a .txt file of list items as strings. Note that
    each item in this file must be on a separate line.
    
    Arguments:
        filename: The name of the .txt file that contains the list.
    '''
    ls = []
    with open(filename) as file:
        for line in file:
            ls.append(line[:-1])
    
    return ls

def load_tsv(folder, chromosomes, dtype = float, ndmin = 0):
    '''
    This function will read in a folder of .tsv files, reading in the file for
    each chromosome provided.
    
    Arguments:
        folder: The directory of .tsv files to load.
        chromosomes: A list of chromosomes. Must correspond to the file names.
        dtype: The data type to read in. Passed through to np.loadtxt().
            Defaults to float.
        ndmin: The minimum number of dimensions to load in for each array.
            Defaults to 0.
    '''
    data = {}
    
    for chrom in chromosomes:
        data[chrom] = np.loadtxt(
            folder + "/" + chrom + ".tsv",
            dtype = dtype,
            delimiter = "\t",
            ndmin = ndmin
        )
    
    return data

def load_dataset(folder, prefixes):
    '''
    This function will read in a preprocessed dataset produced by the
    preprocess_data.py script.
    
    Arguments:
        folder: The top-level directory that the preprocessed data is in.
        prefixes: A list of prefixes of the dataset you want to load.
    
    Note the data will be returned in a tuple, with the following order:
    genome, blocks, prefix0_minus, prefix0_plus, prefix1_minus, prefix1_plus,
    etc.
    '''
    print("~~~ Loading Dataset ~~~")
    
    
    chroms = []
    with open(folder + "/genome/chromosomes.txt") as chromfile:
        for line in chromfile:
            chroms.append(line[:-1])
    
    print("-- genome...")
    genome = load_tsv(folder + "/genome/array", chroms, dtype = str)
    
    print("-- blocks...")
    blocks = load_tsv(folder + "/blocks", chroms)
    
    for chrom in blocks:
        blocks[chrom] = np.array(blocks[chrom], dtype = int)
    
    print("-- datasets (- strand)...")
    datasets_minus = [
        load_tsv(folder + "/" + prefix + "_minus", chroms)
        for prefix in prefixes
    ]
    
    print("-- datasets (+ strand)...")
    datasets_plus = [
        load_tsv(folder + "/" + prefix + "_plus", chroms)
        for prefix in prefixes
    ]
    
    print("- Compiling...")
    returnList = [genome, blocks]
    for i in range(len(prefixes)):
        returnList.append(datasets_minus[i])
        returnList.append(datasets_plus[i])
    
    print("~~~ Loaded Dataset ~~~")
    return tuple(returnList)

def load_json(filename):
    '''
    This function will load a JSON file.
    
    Arguments:
        filename: The name of the JSON file to load.
    '''
    with open(filename) as file:
        data = json.load(file)
    
    return data

#%% Data Outputs

def output_list(ls, destination_file):
    '''
    This function will take in a list of data and spit out a .txt file with
    each item on a separate line.
    
    Arguments:
        ls: The list of items to output.
        destination_file: The .txt file to write to.
    
    Note that the items are passed through to the str() function first.
    '''
    with open(destination_file, "w") as file:
        for item in ls:
            file.write(str(item) + "\n")

def output_tsv(data, out_dir, fmt="%.18e", chroms=False, mkdir=True):
    '''
    This function will convert a data dictionary to an output .tsv file.
    
    Arguments:
        data: A dictionary of data.
        out_dir: The directory to dump the output into.
        fmt: The format string to pass into np.savetxt. Defaults to '%18e'.
        chroms: If True, will also create a 'chromosomes.txt' file
            in the directory with all the data dictionary keys. Defaults to
            False.
        mkdir: If True, will create the outdir directory if it doesn't already
            exist. Defaults to True.
    '''
    if (mkdir):
        if (not os.path.isdir(out_dir)):
            os.mkdir(out_dir)
    
    for chrom in data:
        np.savetxt(
            out_dir + "/" + chrom + ".tsv",
            data[chrom],
            delimiter = "\t",
            fmt = fmt
        )
    
    if (chroms):
        output_list(list(data.keys()), out_dir + "/chromosomes.txt")

def output_wig(data, filename, adjustIndex = True):
    '''
    This function will convert a data dictionary to an output .wig file for
    viewing in IGV. Will convert the locations to ints.
    
    Arguments:
        data: A dictionary of data, must contain a numpy array with 2 columns,
          position (zero-indexed) and value.
        filename: The name of the file to write the output to.
        adjustIndex: True to adjust the index from zero-index to one-index.
            Defaults to True.
    '''
    with open(filename, "w") as file:
        for chrom in data:
            file.write("variableStep chrom=" + chrom + " span=1\n")
            size = np.shape(data[chrom])[0]
            
            string = "".join(
                np.char.add(
                    np.char.add(
                        (
                            data[chrom][:,0] + int(adjustIndex)
                        ).astype(int).astype(str),
                        np.repeat("\t",size)
                    ),
                    np.char.add(
                        data[chrom][:,1].astype(str),
                        np.repeat("\n", size)
                    )
                )
            )
            
            file.write(string)

def output_wig_raw(data, filename):
    '''
    This function will convert a data dictionary to an output .wig file
    without making ANY adjustments to the data (will leave locations as
    floating-point numbers).
    
    Arguments:
        data: A dictionary of data, must contain a numpy array with 2 columns:
          position and value.
        filename: The name of the file to write the output to.
    
    Note that no adjustments are made to the position (it is left in whateever
    index you left it in).
    '''
    with open(filename, "w") as file:
        for chrom in data:
            file.write("variableStep chrom=" + chrom + " span=1\n")
            n = np.shape(data[chrom])[0]
            
            string = "".join(
                np.char.add(
                    np.char.add(
                        data[chrom][:,0].astype(str),
                        np.repeat("\t", n)
                    ),
                    np.char.add(
                        data[chrom][:,1].astype(str),
                        np.repeat("\n", n)
                    )
                )
            )
            
            file.write(string)

def output_wig_from_table(table, filename):
    '''
    This function will take in a table with three columns: chromosome,
    position, and value, then create a .wig file from that as opposed to a
    dictionary.
    
    Arguments:
        table: The table of data. Column 0 must be the chromosome, column 1
            must be the location, and the column 2 must be the value. Must
            be a numpy array of strings.
        filename: The name of the file to create.
    '''
    with open(filename, "w") as file:
        for chrom in np.unique(table[:,0]):
            file.write("variableStep chrom=" + chrom + " span=1\n")
            selected = table[:,0] == chrom
            n = np.sum(selected)
            string = "".join(
                np.char.add(
                    np.char.add(
                        (table[selected,1].astype(int) + 1).astype(str),
                        np.repeat("\t", n)
                    ),
                    np.char.add(table[selected,2], np.repeat("\n", n))
                )
            )
            file.write(string)

def output_wig_interval_from_table(table, filename, highlight = False):
    '''
    This function will take in a table with three columns: chromosome,
    start position, and stop position, then create a .wig file from that where
    there are 1s placed at every position between the start and stop positions
    (inclusive).
    
    Arguments:
        table: The table of data. Column 0 must be the chromosome, column 1
            must be the start location of the interval to mark, and column 2
            must be the stop location (inclusive) of the interval to mark.
            Must be a numpy array of strings.
        filename: The name of the file to create.
        highlight: If True, will take a location in column 3 of table and make
            the value at that location ten instead of one. Defaults to False.
    '''
    with open(filename, "w") as file:
        for chrom in np.unique(table[:,0]):
            file.write("variableStep chrom=" + chrom + " span=1\n")
            selected = table[:,0] == chrom
            
            positions = table[selected,1][:,np.newaxis].astype(int) + 1
            stops = table[selected,2].astype(int) + 1
            if (highlight):
                highlights = table[selected,3].astype(int) + 1
            
            while (True in (positions[:,-1] < stops)):
                positions = np.column_stack((positions, positions[:,-1] + 1))
            positions = np.unique(positions.flatten())
            
            valString = np.repeat("\t1\n", len(positions)).astype("<U32")
            if (highlight):
                valString[np.isin(positions, highlights)] = "\t10\n"
            txt = "".join(
                np.char.add(
                    positions.astype(str),
                    valString
                )
            )
            file.write(txt)

def output_json(filename, data):
    '''
    This function dumps a data object to a JSON file.
    
    Arguments:
        filename: The name of the JSON file to dump the data into.
        data: The object to save in JSON format. Must be JSON-serializable.
    '''
    with open(filename, "w") as file:
        json.dump(data, file, indent = 4)

def output_fasta(filename, seqs, names = None):
    '''
    This function creates a FASTA file from a list of sequences.
    
    Arguments:
        filename: The name of the file to create.
        seqs: The list of sequences to store.
        names: The names of the sequences. If None, will instead name each
            sequence as "Seq<n>", where n is the index in the list. Defaults
            to None.
    '''
    if (names is None):
        Nms = ["Seq" + str(n) for n in range(len(seqs))]
    else:
        Nms = names
    
    with open(filename, "w", newline = "\n") as file:
        for i in range(len(seqs)):
            file.write(">" + Nms[i] + "\n")
            file.write(seqs[i] + "\n")

#%% Data Processing

def remove_bad_blocks(minus_data, plus_data, blocks):
    '''
    This function will remove all the bad data from the datasets based on the
    supplied dictionary of blocks.
    
    Arguments:
        minus_data: A list of data dictionaries for the minus strand.
        plus_data: A list of data dictionaries for the plus strand.
        blocks: The dictionary of blocks. Each key must be a chromosome, and
            each value must be a numpy array with three columns: block
            starting position, block ending position, and whether or not the
            block is bad.
    '''
    sets = len(minus_data)
    if (len(plus_data) != sets):
        raise Exception(
            "plus_data and minus_data must have same number of entries!"
        )
    
    for chrom in blocks:
        if (1 in blocks[chrom][:,2]):
            badBlocks = blocks[chrom][:,2] == 1
            badStarts = blocks[chrom][badBlocks,0]
            badEnds = blocks[chrom][badBlocks,1]
            for i in range(len(badStarts)):
                for s in range(sets):
                    minus_data[s][chrom] = minus_data[s][chrom][
                        np.logical_or(
                            minus_data[s][chrom][:,0] < badStarts[i],
                            minus_data[s][chrom][:,0] > badEnds[i]
                        ),
                        :
                    ]
                    plus_data[s][chrom] = plus_data[s][chrom][
                        np.logical_or(
                            plus_data[s][chrom][:,0] < badStarts[i],
                            plus_data[s][chrom][:,0] > badEnds[i]
                        ),
                        :
                    ]
    
    return (minus_data, plus_data)

def cleaning_pass(minus_data, plus_data, blocks, stds):
    '''
    This function will take in a list of data for both strands and remove
    outlier regions. The algorithm works by sliding a block of width blockW
    along the genome in intervals of width step and deleting the data in that
    region if the mean of the data in that region exceeds the mean of the
    complete dataset across the whole genome by at least stds standard
    deviations. This is useful for removing repetitive DNA regions from the
    dataset.
    
    Arguments:
        minus_data: A list of data dictionaries for the minus strand.
        plus_data: A list of data dictionaries for the plus strand.
        blocks: The dictionary of blocks. Each key must be a chromosome, and
            each value must be a numpy array with three columns: block
            starting position, block ending position, and whether or not the
            block is bad.
        stds: How many standard deviations by which a region's mean must
            exceed the global mean to be removed.
    
    Note that this function returns a tuple containing the cleaned minus_data,
    the cleaned plus_data, and the updated blocks dictionary.
    '''
    print("~~~ Cleaning Pass ~~~")
    sets = len(minus_data)
    if (len(plus_data) != sets):
        raise Exception(
            "plus_data and minus_data must have same number of entries!"
        )
    # Aggregate data
    print("- Aggregating data...")
    full_data = np.array([])
    for i in range(sets):
        for chrom in blocks:
            full_data = np.concatenate((full_data, minus_data[i][chrom][:,1]))
            full_data = np.concatenate((full_data, plus_data[i][chrom][:,1]))
    mean = np.mean(full_data)
    std = np.std(full_data)
    print("-- Mean: " + str(mean) + "; Std: " + str(std))
    
    # Find bad blocks
    print("- Finding bad blocks...")
    for s in range(sets):
        print("-- Dataset " + str(s + 1) + "/" + str(sets) + "...")
        for chrom in blocks:
            print("--- " + chrom + "...")
            for i in range(np.shape(blocks[chrom])[0]):
                if (blocks[chrom][i,2] == 0):
                    start = blocks[chrom][i,0]
                    end = blocks[chrom][i,1]
                    data = np.vstack(
                        (minus_data[s][chrom], plus_data[s][chrom])
                    )
                    data = data[
                        np.logical_and(data[:,0]>=start,data[:,0]<=end),
                        1
                    ]
                    boxMean = np.mean(data)
                    if (boxMean > mean + stds*std):
                        blocks[chrom][i,2] = 1
                        print("==> Bad box at "+str(start)+": "+str(boxMean))
    
    # Remove bad blocks
    print("- Removing bad blocks...")
    minus_data, plus_data = remove_bad_blocks(minus_data, plus_data, blocks)
    """
    for chrom in blocks:
        if (1 in blocks[chrom][:,2]):
            badBlocks = blocks[chrom][:,2] == 1
            badStarts = blocks[chrom][badBlocks,0]
            badEnds = blocks[chrom][badBlocks,1]
            for i in range(len(badStarts)):
                for s in range(sets):
                    minus_data[s][chrom] = minus_data[s][chrom][
                        np.logical_or(
                            minus_data[s][chrom][:,0] < badStarts[i],
                            minus_data[s][chrom][:,0] > badEnds[i]
                        ),
                        :
                    ]
                    plus_data[s][chrom] = plus_data[s][chrom][
                        np.logical_or(
                            plus_data[s][chrom][:,0] < badStarts[i],
                            plus_data[s][chrom][:,0] > badEnds[i]
                        ),
                        :
                    ]
    """
    
    print("~~~ Cleaning Pass Completed ~~~")
    return (minus_data, plus_data, blocks)

def cull_data(minus_data, plus_data, chromosomes, threshold):
    '''
    This function will go through the data and cull all data points for which
    all the datasets have values below a specificied threshold.
    
    Arguments:
        minus_data: A list of data dictionaries for the minus strand.
        plus_data: A list of data dictionaries for the plus strand.
        chromosomes: A list of chromosomes to cull the data for.
        threshold: The culling threshold.
    
    Note this function returns a tuple containing the culled minus_data and
    culled plus_data.
    '''
    sets = len(minus_data)
    if (len(plus_data) != sets):
        raise Exception(
            "plus_data and minus_data must have same number of entries!"
        )
    
    for chrom in chromosomes:
        all_vals = np.column_stack(
            tuple([minus_data[i][chrom][:,1] for i in range(sets)])
        )
        bad = np.sum(all_vals < threshold, axis = 1) == sets
        for i in range(sets):
            minus_data[i][chrom][np.logical_not[bad],:]
        
        all_vals = np.column_stack(
            tuple([plus_data[i][chrom][:,1] for i in range(sets)])
        )
        bad = np.sum(all_vals < threshold, axis = 1) == sets
        for i in range(sets):
            plus_data[i][chrom][np.logical_not[bad],:]
    
    return (minus_data, plus_data)

def floor_data(minus_data, plus_data, chromosomes, floor):
    '''
    This function will go through the data and set all data points for which
    the value is below a specificied floor value to that floor value.
    
    Arguments:
        minus_data: A list of data dictionaries for the minus strand.
        plus_data: A list of data dictionaries for the plus strand.
        chromosomes: A list of chromosomes to cull the data for.
        floor: The floor value.
    
    Note this function returns a tuple containing the culled minus_data and
    culled plus_data.
    '''
    sets = len(minus_data)
    if (len(plus_data) != sets):
        raise Exception(
            "plus_data and minus_data must have same number of entries!"
        )
    
    for chrom in chromosomes:
        for i in range(sets):
            minus_data[i][chrom][minus_data[i][chrom][:,1] < floor, 1] = floor
            plus_data[i][chrom][plus_data[i][chrom][:,1] < floor, 1] = floor
    
    return (minus_data, plus_data)

def arcsinh_transform_data(minus_data, plus_data, chromosomes):
    '''
    This function will go through the data and arcsinh-transform all the data
    points.
    
    Arguments:
        minus_data: A list of data dictionaries for the minus strand.
        plus_data: A list of data dictionaries for the plus strand.
        chromosomes: A list of chromosomes to transform the data for.
    
    Note this function returns a tuple containing the transformed minus_data
    and transformed plus_data.
    '''
    sets = len(minus_data)
    if (len(plus_data) != sets):
        raise Exception(
            "plus_data and minus_data must have same number of entries!"
        )
    for chrom in chromosomes:
        for i in range(sets):
            minus_data[i][chrom][:,1] = np.arcsinh(minus_data[i][chrom][:,1])
            plus_data[i][chrom][:,1] = np.arcsinh(plus_data[i][chrom][:,1])
    
    return (minus_data, plus_data)

def preprocess_data(
    minus_filenames,
    plus_filenames,
    genome_filename,
    out_dir,
    prefixes,
    passes,
    stds,
    blockW,
    step,
    cull_mode,
    threshold,
    excluded_regions_filename = "",
    excluded_regions_adjustIndex = True,
    transformation = 'none',
    save = True,
    verbose = False
):
    '''
    This function will preprocess raw .wig and .fasta files to much better
    .tsv files, cleaning the data to remove outlier regions and culling/
    flooring the data.
    
    Arguments:
        minus_filenames: List of file names for minus strand .wig files.
        plus_filenames: List of file names for plus strand .wig files.
        genome_filename: The file name for the genome .fasta file.
        out_dir: The output directory.
        prefixes: List of prefixes for the datasets.
        passes: How many cleaning passes to subject the data to.
        stds: How many standard deviations a region's mean must exceed the
            global mean to be considered an outlier.
        cull_mode: 'cull' or 'floor' to cull or floor the data, respectively.
        threshold: Culling threshold or data floor value, based on cull_mode.
        excluded_regions_filename: The name of a .txt file containing genomic
            coordinates of regions to exclude from the analysis. Useful for
            removing marker genes or deletion mutants. Defaults to ''.
        excluded_regions_adjustIndex: If using a .txt file of regions to
            exclude, this variable is used to determine whether or not to
            adjust the index of the coordinates listed in the .txt file. If
            this value is True, it will subtract 1 from the listed coordinates
            to 0-index them. Defaults to True.
        transformation: A string relating to what transformation to subject
            the data to. Defaults to 'none'. Current options are:
            'arcsinh'
        save: True to save the output to files in addition to returning it as
            tuples. Defaults to True.
        verbose: True to output lots more information, such as dipyrimidine
            locations, dinucleotide arrays, etc. Defaults to False.
    
    Note this function returns a tuple containing the reformatted minus_data,
    reformatted plus_data, blocks, genome array, dinucleotides (verbose),
    minus dipyrimidines (verbose), and plus dipyrimidines (verbose).
    
    Also note that the three arguments minus_filenames, plus_filenames, and
    prefixes must be corresponding. So, if you have two datasets, say one for
    naked DNA and one for cellular DNA, you'd need to format the arguments
    as lists like follows:
        minus_filenames = ["naked_data_minus.wig", "cellular_data_minus.wig"]
        plus_filenames = ["naked_data_plus.wig", "cellular_data_plus.wig"]
        prefixes = ["nk", "cl"]
    '''
    print("~~~ Preprocessing Data ~~~")
    sets = len(minus_filenames)
    if (len(plus_filenames) != sets) or (len(prefixes) != sets):
        raise Exception(
            "minus_filenames, plus_filenames, prefixes must be the same len!"
        )
    
    print("- Step 1: Reformatting...")

    # Get genome stuff
    print("-- Genome...")
    print("--- Reading raw genome...")
    str_genome = read_raw_genome(genome_filename)
    print("--- Converting genome...")
    genome = convert_genome_to_array(str_genome)
    print("--- Getting dinucleotides...")
    dinucleotides = get_dinucleotides(genome)
    print("--- Finding dipyrimidines...")
    minus_dipyrimidines = get_dipyrimidine_locations(dinucleotides, "-")
    plus_dipyrimidines = get_dipyrimidine_locations(dinucleotides, "+")
    
    chroms = list(genome.keys())

    print("-- Reading in data...")
    minus_data = []
    plus_data = []
    for i in range(sets):
        print("--- " + minus_filenames[i] + "...")
        raw = read_raw_wig(minus_filenames[i])
        minus_data.append(reformat_wig_data(raw, minus_dipyrimidines))
        
        print("--- " + plus_filenames[i] + "...")
        raw = read_raw_wig(plus_filenames[i])
        plus_data.append(reformat_wig_data(raw, plus_dipyrimidines))
    
    if (excluded_regions_filename != ""):
        excluded_regions = read_intervals(
            excluded_regions_filename,
            excluded_regions_adjustIndex
        )
    else:
        excluded_regions = {}
    
    # Data masking
    print("- Step 2: Cleaning Passes...")
    blocks = {}
    for chrom in genome:
        blockStarts = np.arange(0, len(genome[chrom]) - blockW + step, step)
        blocks[chrom] = np.column_stack(
            (
                blockStarts,
                np.clip(blockStarts + blockW - 1, 0, len(genome[chrom]) - 1),
                np.zeros(len(blockStarts))
            )
        )
    
    if (len(excluded_regions) > 0):
        print("-- Excluding regions...")
        for chrom in excluded_regions:
            blocks[chrom] = np.vstack(
                (
                    blocks[chrom],
                    np.column_stack(
                        (
                            excluded_regions[chrom],
                            np.ones(np.shape(excluded_regions[chrom])[0])
                        )
                    )
                )
            )
        minus_data, plus_data = remove_bad_blocks(minus_data,plus_data,blocks)
    
    for I in range(passes):
        print("-- Pass " + str(I + 1) + "...")
        minus_data, plus_data, blocks = cleaning_pass(
            minus_data,
            plus_data,
            blocks,
            stds
        )

    # Data culling
    print("- Step 3: Culling...")
    if (cull_mode == "cull"):
        print("-- Culling data...")
        minus_data,plus_data=cull_data(minus_data,plus_data,chroms,threshold)
    elif (cull_mode == "floor"):
        print("-- Flooring data...")
        minus_data,plus_data=floor_data(minus_data,plus_data,chroms,threshold)
    else:
        print("-- No valid cull_mode specified. Skipping this step...")

    # Data transformation
    print("- Step 4: Data transformation")
    if (transformation == "arcsinh"):
        print("-- arcsinh transformation...")
        minus_data, plus_data = arcsinh_transform_data(
            minus_data,
            plus_data,
            chroms
        )
    else:
        print("-- No valid transformation specified. Skipping this step...")

    # Save outputs
    if (save):
        print("- Step 5: Saving...")
        if (not os.path.isdir(out_dir)):
            os.mkdir(out_dir)
    
        print("-- Data...")
        for i in range(sets):
            print("--- " + prefixes[i] + "...")
            
            print("= - strand...")
            minus_dir = out_dir + "/" + prefixes[i] + "_minus"
            print("== .tsv file...")
            output_tsv(minus_data[i], minus_dir)
            print("== .wig file...")
            output_wig(minus_data[i], out_dir+"/"+prefixes[i]+"_minus.wig")
            
            print("= + strand...")
            plus_dir = out_dir + "/" + prefixes[i] + "_plus"
            print("== .tsv file...")
            output_tsv(plus_data[i], plus_dir)
            print("== .wig file...")
            output_wig(plus_data[i], out_dir+"/"+prefixes[i]+"_plus.wig")
    
        print("-- Raw Data...")
        for i in range(sets):
            output_wig(
                read_raw_wig(minus_filenames[i]),
                out_dir + "/raw_" + prefixes[i] + "_minus.wig",
                False
            )
            output_wig(
                read_raw_wig(plus_filenames[i]),
                out_dir + "/raw_" + prefixes[i] + "_plus.wig",
                False
            )
    
        print("-- Blocks...")
        output_tsv(blocks, out_dir + "/blocks")
    
        print("-- Genome...")
        if (not os.path.isdir(out_dir + "/genome")):
            os.mkdir(out_dir + "/genome")
        
        with open(out_dir + "/genome/chromosomes.txt", "w") as chromfile:
            for chrom in chroms:
                chromfile.write(chrom + "\n")
        
        output_tsv(genome, out_dir + "/genome/array", fmt = "%s")
        
        # Verbose stuff
        if (verbose):
            print("-- Verbose...")
            output_tsv(dinucleotides,out_dir+"/genome/dinucleotides",fmt="%s")
            
            output_tsv(
                minus_dipyrimidines,
                out_dir + "/genome/dipyrimidines_minus",
                fmt = "%s"
            )
            
            output_tsv(
                plus_dipyrimidines,
                out_dir + "/genome/dipyrimidines_plus",
                fmt = "%s"
            )
    # Final output
    print("~~~ Preprocessed Data ~~~")
    if (verbose):
        return (
            minus_data,
            plus_data,
            blocks,
            genome,
            dinucleotides,
            minus_dipyrimidines,
            plus_dipyrimidines
        )
    else:
        return (minus_data, plus_data, blocks, genome)

def normalize_wig_files(
    minus_files,
    plus_files,
    base_minus,
    base_plus,
    out_name = "data",
    save = True
):
    '''
    This function will take in a list of .wig files and multiply the read
    counts in the input files by a factor to bring the total read counts up
    to the same level as that found in the base dataset.
    
    Arguments:
        minus_files: List of .wig files of minus strand data for replicate
            experiments.
        plus_files: List of .wig files of plus strand data for replicate
            experiments.
        base_minus: A .wig file of minus strand data of the base level to
            normalize the replicate files to.
        base_plus: A .wig file of plus strand data of the base level to
            normalize the replicate files to.
        out_name: The name of the output .wig files. Defaults to 'data'.
        save: Whether or not to save the .wig files. Defaults to True.
    
    Note that this function will output a .wig file that is NOT adjusted
    in any way, and the data output will be as if it was just read in from
    read_raw_wig(). Will append '_norm_<strand>.wig' to the out_name when
    forming the two output files.
    '''
    print("~~~ Normalizing Data ~~~")
    
    if (len(minus_files) != len(plus_files)):
        raise Exception("minus_files and plus_files must be the same length!")

    print("- Reading minus reps...")
    minus_reps = [read_raw_wig(fname) for fname in minus_files]
    print("- Reading plus reps...")
    plus_reps = [read_raw_wig(fname) for fname in plus_files]
    print("- Reading normalized minus...")
    norm_minus = read_raw_wig(base_minus)
    print("- Reading normalized plus...")
    norm_plus = read_raw_wig(base_plus)

    # Combine reps
    print("- Combining reps...")

    # Get chromosomes
    chroms = []
    for rep in minus_reps:
        for chrom in rep:
            if (chrom not in chroms):
                chroms.append(chrom)
    for rep in plus_reps:
        for chrom in rep:
            if (chrom not in chroms):
                chroms.append(chrom)

    # Get locs
    print("-- Getting locations...")
    minus_locs = {}
    for chrom in chroms:
        theseLocs = np.concatenate(
            tuple([rep[chrom][:,0] for rep in minus_reps])
        )
        minus_locs[chrom] = np.sort(np.unique(theseLocs))
    plus_locs = {}
    for chrom in chroms:
        theseLocs = np.concatenate(
            tuple([rep[chrom][:,0] for rep in plus_reps])
        )
        plus_locs[chrom] = np.sort(np.unique(theseLocs))

    # Get damage values
    print("-- Getting damage values...")
    minus_dmg = {}
    for chrom in chroms:
        minus_dmg[chrom] = np.zeros(len(minus_locs[chrom]))
        for rep in minus_reps:
            theseLocs = rep[chrom][:,0]
            theseDmgs = rep[chrom][:,1]
            dummy, theseIndices, locsIndices = np.intersect1d(
                theseLocs,
                minus_locs[chrom],
                return_indices = True
            )
            minus_dmg[chrom][locsIndices] += theseDmgs[theseIndices]
    plus_dmg = {}
    for chrom in plus_reps[0]:
        plus_dmg[chrom] = np.zeros(len(plus_locs[chrom]))
        for rep in plus_reps:
            theseLocs = rep[chrom][:,0]
            theseDmgs = rep[chrom][:,1]
            dummy, theseIndices, locsIndices = np.intersect1d(
                theseLocs,
                plus_locs[chrom],
                return_indices = True
            )
            plus_dmg[chrom][locsIndices] += theseDmgs[theseIndices]

    # Normalization
    print("- Normalization...")
    target_total = 0
    for chrom in norm_minus:
        target_total += np.sum(norm_minus[chrom][:,1])
    for chrom in norm_plus:
        target_total += np.sum(norm_plus[chrom][:,1])

    this_total = 0
    for rep in minus_reps:
        for chrom in rep:
            this_total += np.sum(rep[chrom][:,1])
    for rep in plus_reps:
        for chrom in rep:
            this_total += np.sum(rep[chrom][:,1])

    factor = target_total / this_total

    for chrom in minus_dmg:
        minus_dmg[chrom] *= factor
    for chrom in plus_dmg:
        plus_dmg[chrom] *= factor

    # Compile to final state
    print("- Compiling...")
    minus = {}
    plus = {}
    for chrom in chroms:
        minus[chrom] = np.column_stack((minus_locs[chrom], minus_dmg[chrom]))
        plus[chrom] = np.column_stack((plus_locs[chrom], plus_dmg[chrom]))

    # Output
    print("- Saving...")
    if (save):
        output_wig_raw(minus, out_name + "_norm_minus.wig")
        output_wig_raw(plus, out_name + "_norm_plus.wig")
    
    print("~~~ Normalized Data ~~~")
    return (minus, plus)

def get_bad_locations(locations, blocks):
    '''
    This function returns a dictionary of true/false arrays each of length
    equal to the length of the array stored in locations at that same key.
    True means the location is inside one of the bad blocks.
    
    Arguments:
        locations: A dictionary of arrays holding genomic coordinates
        blocks: The blocks dictionary of block locations
    '''
    returnDict = {}
    for chrom in locations:
        badBlocks = blocks[chrom][blocks[chrom][:,2] == 1,0:2]
        
        badLocs = np.zeros(len(locations[chrom])).astype(bool)
        if (np.shape(badBlocks)[0] > 0):
            for i in range(np.shape(badBlocks)[0]):
                badLocs = np.logical_or(
                    badLocs,
                    np.logical_and(
                        locations[chrom] >= badBlocks[i,0],
                        locations[chrom] <= badBlocks[i,1]
                    )
                )
        
        returnDict[chrom] = badLocs
    
    return returnDict

def get_bad_intervals(intervals, blocks):
    '''
    This function returns a dictionary of true/false arrays each of length
    equal to the length of the array stored in intervals at that same key.
    True means the location is inside one of the bad blocks.
    
    Arguments:
        intervals: A dictionary of arrays holding start and end coordinates
        blocks: The blocks dictionary of block locations
    '''
    returnDict = {}
    
    starts = {}
    ends = {}
    
    for chrom in intervals:
        starts[chrom] = intervals[chrom][:,0]
        ends[chrom] = intervals[chrom][:,1]
    
    badStarts = get_bad_locations(starts, blocks)
    badEnds = get_bad_locations(ends, blocks)
    badInsides = {}
    for chrom in intervals:
        badBlocks = blocks[chrom][blocks[chrom][:,2] == 1,0:2]
        
        badInts = np.zeros(np.shape(intervals[chrom])[0]).astype(bool)
        if (np.shape(badBlocks)[0] > 0):
            for i in range(np.shape(badBlocks)[0]):
                badInts = np.logical_or(
                    badInts,
                    np.logical_and(
                        starts[chrom] <= badBlocks[i,0],
                        ends[chrom] >= badBlocks[i,1]
                    )
                )
        badInsides[chrom] = badInts
    
    for chrom in intervals:
        returnDict[chrom] = np.logical_or(
            badStarts[chrom],
            np.logical_or(
                badEnds[chrom],
                badInsides[chrom]
            )
        )
    
    return returnDict

def make_parallel_dataset(
    minus_wigs,
    plus_wigs,
    prefixes,
    minus_dipyrimidines,
    plus_dipyrimidines,
    blocks,
    out_dir
):
    '''
    This function will take in a separate dataset and make a 'parallel'
    dataset, one where the same blocks have been removed.
    
    Arguments:
        minus_wigs: A list of .wig files for the minus-stranded data.
        plus_wigs: A list of .wig files for the plus-stranded data.
        prefixes: A list of prefixes to name the outputs.
        minus_dipyrimidines: the dipyrimidine locations dictionary for the
            minus strand.
        plus_dipyrimidines: The dipyrimidine locations dictionary for the plus
            strand.
        blocks: The blocks dictionary to remove data points based on.
        out_dir: The name of the directory to put the results into.
    
    Note that this function returns a tuple containing (minus, plus) lists of
    parallel datasets.
    '''
    print("~~~ Making Parallel Dataset ~~~")
    N = len(plus_wigs)
    if (len(minus_wigs) != N):
        raise Exception("minus_wigs must have same length as plus_wigs!")
    if (len(prefixes) != N):
        raise Exception("Must have as many prefixes as datasets!")
    
    print("- Loading datasets...")
    plus_data = [read_raw_wig(plus_wig) for plus_wig in plus_wigs]
    minus_data = [read_raw_wig(minus_wig) for minus_wig in minus_wigs]
    
    print("- Reformatting datasets...")
    plus_data = [
        reformat_wig_data(raw, plus_dipyrimidines)
        for raw in plus_data
    ]
    minus_data = [
        reformat_wig_data(raw, minus_dipyrimidines)
        for raw in minus_data
    ]
    
    print("- Removing bad positions...")
    plus_locs = [{} for i in range(N)]
    minus_locs = [{} for i in range(N)]
    for chrom in plus_dipyrimidines:
        for i in range(N):
            plus_locs[i][chrom] = plus_data[i][chrom][:,0]
            minus_locs[i][chrom] = minus_data[i][chrom][:,0]
    
    bad_plus = [get_bad_locations(loc, blocks) for loc in plus_locs]
    bad_minus = [get_bad_locations(loc, blocks) for loc in minus_locs]
    
    for chrom in plus_dipyrimidines:
        for i in range(N):
            plus_data[i][chrom] = plus_data[i][chrom][
                np.logical_not(bad_plus[i][chrom]),
                :
            ]
            minus_data[i][chrom] = minus_data[i][chrom][
                np.logical_not(bad_minus[i][chrom]),
                :
            ]
    
    print("- Outputting...")
    for i in range(N):
        tsv_dir = out_dir + "/" + prefixes[i] + "_plus"
        output_tsv(plus_data[i], tsv_dir)
        output_wig(plus_data[i], out_dir + "/" + prefixes[i] + "_plus.wig")
        
        tsv_dir = out_dir + "/" + prefixes[i] + "_minus"
        output_tsv(minus_data[i], tsv_dir)
        output_wig(minus_data[i], out_dir + "/" + prefixes[i] + "_minus.wig")
    
    print("~~~ Made Parallel Dataset ~~~")
    return (minus_data, plus_data)

def preprocess_supplementary_wigs(wigFiles, blocks, out_names):
    '''
    This function will take in a list of .wig files of supplementary data and
    convert that data to a zero-indexed dictionary of numpy arrays and output
    corresponding .tsv and .wig files with bad blocks removed.
    
    Arguments:
        wigFiles: A list of names of .wig files of supplementary data.
        blocks: The array of blocks.
        out_names: A list of names that will serve as both the directory names
            for outputting the .tsv files and the filenames of the .wig files.
    '''
    N = len(wigFiles)
    if (len(out_names) != N):
        raise Exception("Number of out_dirs and wigFiles must be the same!")
    
    returnData = []
        
    for i in range(N):
        wig_data = read_raw_wig(wigFiles[i])
        locs = {}
        for chrom in wig_data:
            wig_data[chrom][:,0] -= 1 # Zero-index
            locs[chrom] = wig_data[chrom][:,0]
        badLocs = get_bad_locations(locs, blocks)
        
        for chrom in wig_data:
            wig_data[chrom]=wig_data[chrom][np.logical_not(badLocs[chrom]),:]
        
        output_tsv(wig_data, out_names[i], chroms = True)
        output_wig(wig_data, out_names[i] + ".wig")
        
        returnData.append(wig_data)
    
    return returnData