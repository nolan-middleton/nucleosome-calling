##############################################################################
# This file contains library functions for more advanced data analysis       #
# methods.                                                                   #
##############################################################################

#%% Setup

# Imports
import numpy as np
import scipy.cluster as cluster
import scipy.optimize as optimize
import scipy.integrate as integrate
import sklearn.metrics as metrics
import os as os
from .DataIO import output_tsv, output_wig

#%% Miscellaneous Analysis Functions

def transform_dataset(minus_data, plus_data, function):
    '''
    This function will take in a dataset and transform it according to a given
    function. Note that the function must take in a 1-dimensional numpy array
    as its only argument and return a 1-dimensional numpy array of the same
    length. Returns a tuple (transformed minus data, transformed plus data).
    
    Arguments:
        minus_data: The minus-stranded data.
        plus_data: The plus-stranded data.
        function: The function to transform the data values by.
    '''
    transf_minus = {}
    for chrom in minus_data:
        transf_minus[chrom] = np.column_stack(
            (minus_data[chrom][:,0], function(minus_data[chrom][:,1]))
        )
    
    transf_plus = {}
    for chrom in plus_data:
        transf_plus[chrom] = np.column_stack(
            (plus_data[chrom][:,0], function(plus_data[chrom][:,1]))
        )
    
    return (transf_minus, transf_plus)

def make_transformed_dataset(minus_data, plus_data, function, prefix, out_dir):
    '''
    This function will take in a dataset and transform it according to a given
    function, then output the results to a directory. Note that the function
    must take in a 1-dimensional numpy array as its only argument and return a
    1-dimensional numpy array of the same length. Returns a tuple (transformed
    minus data, transformed plus data).
    
    Arguments:
        minus_data: The minus-stranded data.
        plus_data: The plus-stranded data.
        function: The function to transform the data values by.
        prefix: The prefix to label the outputs with.
        out_dir: The directory to put the outputs into.
    '''
    print("~~~ Making Transformed Dataset ~~~")
    
    print("- Getting transformed data...")
    transf_minus,transf_plus=transform_dataset(minus_data,plus_data,function)
    
    print("- Outputting...")
    output_tsv(transf_minus, out_dir + "/" + prefix + "_minus")
    output_tsv(transf_plus, out_dir + "/" + prefix + "_plus")
    output_wig(transf_minus, out_dir + "/" + prefix + "_minus.wig")
    output_wig(transf_plus, out_dir + "/" + prefix + "_plus.wig")
    
    print("~~~ Made Transformed Dataset ~~~")
    return (transf_minus, transf_plus)

def find_interval_overlaps(intervals, regions):
    '''
    This function will return the length of the longest overlap between each
    interval and any of the regions.
    
    Arguments:
        intervals: A dictionary of intervals, column 0 = start, column 1 = end
        regions: A dictionary of regions, column 0 = start, column 1 = end
    '''
    overlaps = {}
    
    for chrom in intervals:
        overlaps[chrom] = np.zeros(np.shape(intervals[chrom])[0])
        
        if (chrom in regions):
            a = intervals[chrom][:,0]
            b = intervals[chrom][:,1]
            for i in range(np.shape(regions[chrom])[0]):
                c = regions[chrom][i,0]
                d = regions[chrom][i,1]
                
                overlaps[chrom] = np.max(
                    np.column_stack(
                        (
                            overlaps[chrom],
                            np.clip(
                                np.min(np.column_stack((b-c+1,d-a+1)),axis=1),
                                0,
                                b - a + 1 # I'm assuming and inclusive endpt
                            )
                        )
                    ),
                    axis = 1
                )
        
    return overlaps

def get_inbetween_positions(intervals):
    '''
    This function takes in a dictionary of of start and stop positions and
    returns another dictionary that contains every inbetween position between
    the start and stop positions. Note that the intervals must all be evenly
    spaced. For instance, if there's an interval [1, 5], this function will
    expand that interval to a row containing [1.5, 2.5, 3.5, 4.5].
    
    Arguments:
        intervals: An array of evenly-spaced intervals.
    '''
    if (len(intervals.keys()) == 0):
        raise Exception("You can't give an empty dictionary!")
    
    firstTable = intervals[list(intervals.keys())[0]]
    if (np.shape(firstTable)[0] == 0):
        raise Exception("Empty entry in first key!")
    
    span = firstTable[0,1] - firstTable[0,0]
    
    for chrom in intervals:
        if (False in (intervals[chrom][:,1] - intervals[chrom][:,0] == span)):
            raise Exception("Intervals must be evenly-spaced!")
    
    offsets = np.array(range(span)) + 0.5
    m = len(offsets)
    
    inbetweens = {}
    for chrom in intervals:
        inbetweens[chrom] = np.repeat(
            intervals[chrom][:,0][:,np.newaxis],
            m,
            axis = 1
        ) + np.repeat(
            offsets[np.newaxis,:],
            np.shape(intervals[chrom])[0],
            axis = 0
        )
    
    return inbetweens

def scaled_integral(
    objective_function,
    integrand,
    parameters,
    bounds
):
    '''
    This function will find a scaled integral over a set of bounds using an
    objective function to scale the integral and make it tractable via a log
    factor. Returns (value of the scaled integral, log factor).
    
    Arguments:
        objective_function: The objective function to minimize to rescale the
            integral. The first argument must be the parameter we're
            integrating over, the next arguments must be the arguments in
            parameters in the order they appear in parameters.
        integrand: The integrand function. The first argument must be the
            parameter we're integrating over, the next arguments must be
            the arguments in parameters in the order they appear in parameters,
            and the last argument must be log_factor.
        parameters: A LIST of other arguments for the objective_function and
            integrand that go after the variable we're integrating over and, in
            the case of integrand, before log_factor.
        bounds: A LIST of two elements: the lower and upper bounds for
            the integral.
    '''
    log_factor = -optimize.minimize(
        objective_function,
        x0 = [(bounds[0] + bounds[1]) / 2],
        args = tuple(parameters),
        method = "Powell",
        bounds = [bounds]
    ).fun
    
    scaled_integral = integrate.quad(
        integrand,
        bounds[0],
        bounds[1],
        args = tuple(parameters + [log_factor])
    )[0]
    
    return (scaled_integral, log_factor)

def get_distance_matrix(locs1, locs2, signed = False):
    '''
    This function returns a pairwise distance matrix dictionary between two
    dictionaries of genomic locations.
    
    Arguments:
        locs1: A dictionary of genomic locations. Each key must contain a 1d
            numpy array.
        locs2: A dictionary of genomic locations. Each key must contain a 1d
            numpy array.
        signed: If True, will return the signed distance of arr1 - arr2. Else
            returns the magnitude, the true distance |arr1 - arr2|. Defaults to
            False.
    
    Note that if any chromosome is not shared between the two dictionaries,
    then it will be excluded from the resulting distance matrix dictionary.
    '''
    D = {}
    for chrom in locs1:
        if (chrom in locs2):
            A = np.repeat(locs1[chrom][:,np.newaxis],len(locs2[chrom]),axis=1)
            B = np.repeat(locs2[chrom][np.newaxis,:],len(locs1[chrom]),axis=0)
            
            if (not signed):
                D[chrom] = np.abs(A - B)
            else:
                D[chrom] = A - B
    return D

def get_consecutive_stretches(scores):
    '''
    This function takes in a dictionary of numpy arrays, each with the first
    column being genomic coordinates. It returns a dictionary of python lists,
    where each entry is a continuous stretch of scores.
    
    Arguments:
        scores: A dictionary of numpy arrays. Each entry must have at least 2
            columns, and the first column must be the location.
    '''
    sequences = {}
    for chrom in scores:
        sequences[chrom] = []
        
        theseScores = np.copy(scores[chrom])
        diffs = theseScores[1:,0] - theseScores[:-1,0]
        diffs = np.append([0], diffs)
        
        if (True not in (diffs > 1)):
            sequences[chrom].append(theseScores)
        else:
            while (True in (diffs > 1)):
                stop = np.argmax(diffs > 1)
                sequences[chrom].append(theseScores[:stop,:])
                diffs = diffs[stop:]
                diffs[0] = 0
                theseScores = theseScores[stop:,:]
            sequences[chrom].append(theseScores)
    
    return sequences

#%% Statistics and Classification

def normal(x, mu, sigma):
    '''
    This function will return the value of the normal pdf with mean mu and
    standard devation sigma at the value x.
    
    Arguments:
        x: The value to evaluate the pdf at.
        mu: The mean of the normal distribution.
        sigma: The standard deviation of the normal distribution.
    '''
    return (1/np.sqrt(2*np.pi*(sigma**2)))*np.exp(-((x-mu)**2)/(2*(sigma**2)))

def log_normal(x, mu, sigma):
    '''
    This function will return the value of the log of the normal pdf with mean
    mu and standard devation sigma at the value x.
    
    Arguments:
        x: The value to evaluate the pdf at.
        mu: The mean of the normal distribution.
        sigma: The standard deviation of the normal distribution.
    '''
    return -0.5*np.log(2*np.pi*(sigma**2)) - ((x-mu)**2)/(2*(sigma**2))

def clustering(data, iterations = 5, maxCluster = 7):
    '''
    This function will perform the k-means clustering for the points in the
    damage_matrix, returning the set of labels with the lowest
    silhouette score across all possible k-values over 100 iterations each.
    
    Arguments:
        data: The matrix of data points in space to cluster.
        iterations: How many iterations of the k-means clustering to do.
        maxCluster: Try forming k=2 to this many clusters.
    
    Returns a tuple containing (k, silhouette score, labels)
    '''
    print("~~~ Clustering ~~~")
    scores = [-100 for k in range(maxCluster - 1)]
    all_labels=np.zeros((np.shape(data)[0], maxCluster - 1), dtype = "int32")
    
    for k in range(2, maxCluster + 1):
        print("- k = " + str(k) + "...")
        
        score = -100
        labels = np.zeros(np.shape(data)[0], dtype = "int32")
        for i in range(iterations):
            print("-- Iteration " + str(i + 1) + "...")
            
            print("--- Clustering...")
            dummy, theseLabs = cluster.vq.kmeans2(data, k)
            print("--- Scoring...")
            thisScore = metrics.silhouette_score(data, theseLabs)
            print("--- Silhouette score: " + str(thisScore))
            
            if (thisScore > score):
                score = thisScore
                labels = theseLabs
        all_labels[:,k - 2] = labels
        scores[k - 2] = score
    
    best_k = int(np.argmax(scores) + 2)
    
    print("~~~ Done Clustering ~~~")
    return (best_k, scores[best_k - 2], all_labels[:,best_k - 2])

def accuracy(C):
    '''
    This function takes in a dictionary corresponding to the confusion matrix
    and returns the accuracy.
    
    Arguments:
        C: A dictionary containing the keys "TP", "FP", "TN",
            and "FN" corresponding to True Positives, False Positives, True
            Negatives, and False Negatives, respectively.
    '''
    return (C["TP"] + C["TN"]) / (C["TP"] + C["FP"] + C["TN"] + C["FN"])

def sensitivity(C):
    '''
    This function takes in a dictionary corresponding to the confusion matrix
    and returns the sensitivity (True Positive Rate, i.e. correctly called
    positives over all actual positives). AKA recall.
    
    Arguments:
        C: A dictionary containing the keys "TP", "FP", "TN",
            and "FN" corresponding to True Positives, False Positives, True
            Negatives, and False Negatives, respectively.
    '''
    return C["TP"] / (C["TP"] + C["FN"])

recall = sensitivity
TPR = sensitivity

def specificity(C):
    '''
    This function takes in a dictionary corresponding to the confusion matrix
    and returns the specificity (True Negative Rate, i.e. correctly called
    negatives over all actual negatives).
    
    Arguments:
        C: A dictionary containing the keys "TP", "FP", "TN",
            and "FN" corresponding to True Positives, False Positives, True
            Negatives, and False Negatives, respectively.
    '''
    return C["TN"] / (C["TN"] + C["FP"])

TNR = specificity

def FPR(C):
    '''
    This function takes in a dictionary corresponding to the confusion matrix
    and returns 1 - the specificity (False Positive Rate, i.e. incorrectly
    called positives over all actual negatives).
    
    Arguments:
        C: A dictionary containing the keys "TP", "FP", "TN",
            and "FN" corresponding to True Positives, False Positives, True
            Negatives, and False Negatives, respectively.
    '''
    return 1 - specificity(C)
    
def precision(C):
    '''
    This function takes in a dictionary corresponding to the confusion matrix
    and returns the Precision (Positive Predictive Value, i.e. correctly called
    positives over all guessed positives).
    
    Arguments:
        C: A dictionary containing the keys "TP", "FP", "TN",
            and "FN" corresponding to True Positives, False Positives, True
            Negatives, and False Negatives, respectively.
    '''
    return C["TP"] / (C["TP"] + C["FP"])

PPV = precision

def NPV(C):
    '''
    This function takes in a dictionary corresponding to the confusion matrix
    and returns the Negative Predictive Value, i.e. correctly called negatives
    over all guessed negatives).
    
    Arguments:
        C: A dictionary containing the keys "TP", "FP", "TN",
            and "FN" corresponding to True Positives, False Positives, True
            Negatives, and False Negatives, respectively.
    '''
    return C["TN"] / (C["TN"] + C["FN"])

def phi_coefficient(C):
    '''
    This function takes in a dictionary corresponding to the confusion matrix
    and returns the phi coefficient (or, as the computer scientists INCORRECTLY
    call it, the Matthews Correlation Coefficient).
    
    Arguments:
        C: A dictionary containing the keys "TP", "FP", "TN",
            and "FN" corresponding to True Positives, False Positives, True
            Negatives, and False Negatives, respectively.
    '''
    return (C["TP"]*C["TN"] + C["FP"]*C["FN"]) / np.sqrt(
        (C["TP"]+C["FP"])*(C["TP"]+C["FN"])*(C["TN"]+C["FP"])*(C["TN"]+C["FN"])
    )

# I'm intentionally not making this alias.
# MCC = phi_coefficient

def AUC(false_positive_rates, true_positive_rates):
    '''
    This function takes in two arrays of equal lengths, one for the FPR, or
    1 - specificity, and another for the TPR, or sensitivity, and computes the
    area under the ROC curve.
    
    Arguments:
        false_positive_rates: An array of false positive rates, the x-values of
            the ROC curve.
        true_positive_rates: An array of true positive rates, the y-values of
            the ROC curve.
    '''
    if (len(false_positive_rates) != len(true_positive_rates)):
        raise Exception(
            "false_positive_rates and true_positive_rates must be same length!"
        )
    
    order = np.argsort(false_positive_rates)
    X = false_positive_rates[order]
    Y = true_positive_rates[order]
    
    A = 0
    for i in range(len(X) - 1):
        A += (X[i+1] - X[i]) * Y[i]
    
    return A

#%% Dictionary Handling

def column_split(array_dict):
    '''
    This function will take in a dictionary of arrays and return a tuple of
    dictionaries, where each dictionary has a column from the original
    dictionary instead of the entire array.
    
    Arguments: 
        array_dict: A dictionary of arrays.
    
    Note that the arrays in the dictionary must all have the same number of
    columns.
    '''
    ncols = np.shape(array_dict[list(array_dict.keys())[0]])[1]
    return_list = [{} for i in range(ncols)]
    for chrom in array_dict:
        if (np.shape(array_dict[chrom])[1] != ncols):
            raise Exception("Each entry must have the same number of columns!")
        
        for i in range(ncols):
            return_list[i][chrom] = array_dict[chrom][:,i]
    
    return tuple(return_list)

def column_join(dicts):
    '''
    This function will take in a list of dictionaries of arrays and return a
    tuple dictionary, where each key in the dictionary has the columns from
    the same key from each dictionary in the list stacked together.
    
    Arguments: 
        dicts: A list of dictionaries of arrays.
    
    Note that the arrays in each dictionary must all have the same number of
    rows and the same keys.
    '''
    if (len(dicts) < 1):
        raise Exception("Must supply at least one dictionary!")
    return_dict = {}
    for key in dicts[0]:
        return_dict[key] = np.column_stack(tuple([d[key] for d in dicts]))
    
    return return_dict

def make_type(D, dtype):
    '''
    This function will take in a dictionary of arrays and return the same
    arrays in another dictionary, where each array has had the data type
    changed to dtype.
    
    Arguments: 
        D: A dictionary of arrays.
        dtype: The data type to change the arrays in D to.
    
    Note that the arrays in each dictionary must all be able to be converted
    to dtype.
    '''
    return_dict = {}
    for key in D:
        return_dict[key] = D[key].astype(dtype)
    
    return return_dict

def table_to_dict(table, dtype = str):
    '''
    This function takes in a table of strings and converts it to a dictionary.
    The first column of the table must be the key value for the dict. The rest
    of the columns are data columns and will be condensed into the dictionary
    values.
    
    Arguments:
        table: The table of data. Must be table of strings where the first
            column is the key for the dictionary entry.
        dtype: The data type to convert the rest of the columns to when
            putting them into the dictionary. Defaults to str.
    '''
    return_dict = {}
    for key in np.unique(table[:,0]):
        selected = table[:,0] == key
        return_dict[str(key)] = table[selected,1:].astype(dtype)
    
    return return_dict

def dict_stack(dictionary):
    '''
    This function takes in a dictionary of numpy arrays and returns one giant
    numpy table that is all the entries for each key in the dictionary stacked
    on top of one another.
    
    Arguments:
        dictionary: The dictionary to stack the entries of.
    
    Note that the arrays in each spot of the dictionary must all have the same
    number of columns.
    '''
    if (len(dictionary.keys()) == 0):
        raise Exception("You can't supply an empty dictionary!")
    colnum = np.shape(dictionary[list(dictionary.keys())[0]])[1]
    for key in dictionary:
        if (np.shape(dictionary[key])[1] != colnum):
            raise Exception(
                "All arrays in dictionary must have same number of columns!"
            )
    return np.vstack(tuple([dictionary[key] for key in dictionary]))

def dict_concat(dictionary):
    '''
    This function takes in a dictionary of numpy arrays and returns one giant
    numpy array that is all the entries for each key in the dictionary
    concatenated together.
    
    Arguments:
        dictionary: The dictionary to concatenate the entries of.
    
    Note that the arrays in each spot of the dictionary must all be 1d.
    '''
    if (len(dictionary.keys()) == 0):
        raise Exception("You can't supply an empty dictionary!")
    if (len(np.shape(dictionary[list(dictionary.keys())[0]])) != 1):
        raise Exception("Dictionary entries must be 1d arrays!")
    for key in dictionary:
        if (len(np.shape(dictionary[key])) != 1):
            raise Exception("All arrays in dictionary must be 1d!")
    return np.concatenate(tuple([dictionary[key] for key in dictionary]))

def dict_min(dictionary, axis = None):
    '''
    This function takes in a dictionary of numpy arrays and finds the min of
    each array.
    
    Arguments:
        dictionary: The dictionary of numpy arrays.
        axis: The axis to take the min along. Defaults to None.
    
    Note that all entries in the dictionary must have at least as many axes as
    you specify.
    '''
    return_dict = {}
    for chrom in dictionary:
        return_dict[chrom] = np.min(dictionary[chrom], axis = axis)
    
    return return_dict

def dict_argmin(dictionary, axis = None):
    '''
    This function takes in a dictionary of numpy arrays and finds the argmin of
    each array.
    
    Arguments:
        dictionary: The dictionary of numpy arrays.
        axis: The axis to take the argmin along. Defaults to None.
    
    Note that all entries in the dictionary must have at least as many axes as
    you specify.
    '''
    return_dict = {}
    for chrom in dictionary:
        return_dict[chrom] = np.argmin(dictionary[chrom], axis = axis)
    
    return return_dict

def dict_index(dictionary, index_dict):
    '''
    This function takes in a dictionary of numpy arrays and a dictionary of
    indices and indexes each array of dictionary by the index_dict array.
    
    Arguments:
        dictionary: The dictionary of numpy arrays.
        index_dict: A dictionary of numpy arrays to index dictionary with. The
            keys in index_dict must be a subset of the keys in dictionary, and
            the entries must be boolean numpy arrays or integer numpy arrays to
            index the entry of dictionary at the same key with.
    '''
    return_dict = {}
    for chrom in index_dict:
        return_dict[chrom] = dictionary[chrom][index_dict[chrom]]
    
    return return_dict

def dict_sort(dictionary, axis = -1, kind = None, order = None):
    '''
    This function takes in a dictionary of numpy arrays and sorts each one.
    
    Arguments:
        dictionary: The dictionary of numpy arrays to sort.
        axis: The axis to sort each of the entries along. Defaults to -1.
        kind: The kind of sorting algorithm to use. Defaults to None.
        order: The order of the sorting algorithm. See numpy.sort for details.
            Defaults to None.
    '''
    return_dict = {}
    for chrom in dictionary:
        return_dict[chrom] = np.sort(
            dictionary[chrom],
            axis = axis,
            kind = kind,
            order = order
        )
    
    return return_dict

def dict_argsort(dictionary, axis = -1, kind = None, order = None):
    '''
    This function takes in a dictionary of numpy arrays and returns the argsort
    dictionary for each entry.
    
    Arguments:
        dictionary: The dictionary of numpy arrays to argsort.
        axis: The axis to sort each of the entries along. Defaults to -1.
        kind: The kind of sorting algorithm to use. Defaults to None.
        order: The order of the sorting algorithm. See numpy.sort for details.
            Defaults to None.
    '''
    return_dict = {}
    for chrom in dictionary:
        return_dict[chrom] = np.argsort(
            dictionary[chrom],
            axis = axis,
            kind = kind,
            order = order
        )
    
    return return_dict

def dict_abs(dictionary):
    '''
    This function takes the absolute value of every entry in a dictionary of
    numpy arrays.
    
    Arguments:
        dictionary: The dictionary of numpy arrays to take the absolute values
            of.
    '''
    return_dict = {}
    for chrom in dictionary:
        return_dict[chrom] = np.abs(dictionary[chrom])
    
    return return_dict

def dict_overall_mean(dictionary):
    '''
    This function takes in a dictionary of numpy arrays and returns the global
    dictionary-wide mean of all the numpy arrays.
    
    Arguments: dictionary: A dictionary of numpy arrays to find the overall
        mean of.
    '''
    total_array = np.zeros(0)
    for chrom in dictionary:
        total_array.concatenate(dictionary[chrom].flatten())
    
    return np.mean(total_array)

def dict_floor(dictionary):
    '''
    This function takes the floor of every entry in a dictionary of numpy
    arrays.
    
    Arguments:
        dictionary: The dictionary of numpy arrays to take the floor of.
    '''
    return_dict = {}
    for chrom in dictionary:
        return_dict[chrom] = np.floor(dictionary[chrom])
    
    return return_dict

def dict_compare(dictionary, value, comparison = "=="):
    '''
    This function takes in a dictionary of numpy arrays and a value, then
    returns a dictionary of boolean numpy arrays that are True wherever the
    comparison between the array and value is True.
    
    Arguments:
        dictionary: The dictionary of numpy arrays to search.
        value: The value to look for.
        comparison: String representing the type of comparison to perform. Must
            be one of '==', '>=', '<=', '>', '<', '!=', 'in', or 'not in'.
            Defaults to '=='.
    
    Note that the dictionary entries and value must be compatible with the
    comparison type.
    '''
    if (comparison not in ["==", ">=", "<=", ">", "<", "!=", "in", "not in"]):
        raise Exception("Invalid Comparison Type!")
    
    return_dict = {}
    if (comparison == "=="):
        for chrom in dictionary:
            return_dict[chrom] = dictionary[chrom] == value
    elif (comparison == ">="):
        for chrom in dictionary:
            return_dict[chrom] = dictionary[chrom] >= value
    elif (comparison == "<="):
        for chrom in dictionary:
            return_dict[chrom] = dictionary[chrom] <= value
    elif (comparison == ">"):
        for chrom in dictionary:
            return_dict[chrom] = dictionary[chrom] > value
    elif (comparison == "<"):
        for chrom in dictionary:
            return_dict[chrom] = dictionary[chrom] < value
    elif (comparison == "!="):
        for chrom in dictionary:
            return_dict[chrom] = dictionary[chrom] != value
    elif (comparison == "in"):
        for chrom in dictionary:
            return_dict[chrom] = np.isin(dictionary[chrom], value)
    elif (comparison == "not in"):
        for chrom in dictionary:
            return_dict[chrom] = ~np.isin(dictionary[chrom], value)
    
    return return_dict

def dict_negate(dictionary):
    '''
    This function takes in a dictionary of Boolean numpy arrays and returns the
    logical negation of all the entries.
    
    Arguments:
        dictionary: The dictionary of Boolean numpy arrays to negate.
    '''
    return_dict = {}
    for chrom in dictionary:
        return_dict[chrom] = ~dictionary[chrom]
    
    return return_dict

def dict_sum(dictionaries):
    '''
    This function takes in a list of dictionaries of numpy arrays and returns
    a dictionary that adds the numpy arrays across keys.
    
    Arguments:
        dictionaries: A list of dictionaries of numpy arrays. Each dictionary
            in the list must have the exact same set of keys, and for any key,
            the numpy arrays in each dictionary that correspond to that key
            must have the same shape across all dictionaries in the list.
    '''
    if (len(dictionaries) == 0):
        raise Exception("You must supply at least one dictionary!")
    
    return_dict = {}
    for chrom in dictionaries[0]:
        return_dict[chrom] = np.zeros(np.shape(dictionaries[0][chrom]))
        for D in dictionaries:
            return_dict[chrom] += D[chrom]
    
    return return_dict

#%% Comparisons

def absolute_difference(data_a_minus, data_a_plus, data_b_minus, data_b_plus):
    '''
    This function takes in two datasets (which MUST be reformatted and have
    the same locations) and returns a new dataset containing the difference
    between them.
    
    Arguments:
        data_a_minus: The minus strand data for one dataset, the minuend.
        data_a_plus: The plus strand data for one dataset, the minuend.
        data_b_minus: The minus strand data for the other dataset, the
            subtrahend.
        data_b_plus: The plus strand data for the other dataset, the
            subtrahend.
    
    Note this function returns a tuple containing (minus_diffs, plus_diffs).
    '''
    diffs_plus = {}
    diffs_minus = {}
    for chrom in data_a_plus:
        diffs_plus[chrom] = np.column_stack(
            (
                data_a_plus[chrom][:,0],
                data_a_plus[chrom][:,1] - data_b_plus[chrom][:,1]
            )
        )
        diffs_minus[chrom] = np.column_stack(
            (
                data_a_minus[chrom][:,0],
                data_a_minus[chrom][:,1] - data_b_minus[chrom][:,1]
            )
        )
    
    return (diffs_minus, diffs_plus)

def log2_fold_ratio(
    data_a_minus,
    data_a_plus,
    data_b_minus,
    data_b_plus,
    zeroHandlingMode = "none",
    zeroHandlingValue = 0
):
    '''
    This function takes in two datasets (which MUST be reformatted and have
    the same locations) and returns a new dataset containing the log ratio of
    the values in one dataset over the values in the other dataset.
    
    Arguments:
        data_a_minus: The minus strand data for one dataset, the numerator.
        data_a_plus: The plus strand data for one dataset, the numerator.
        data_b_minus: The minus strand data for the other dataset, the
            denominator.
        data_b_plus: The plus strand data for the other dataset, the
            denominator.
        zeroHandlingMode: Which method to use to handle zero values. Must be
            'none', 'floor', or 'pseudocount'. Defaults to 'none'. If 'none',
            may throw errors if zero values are present in the dataset.
        zeroHandlingValue: The value to use in the zero-handling mode.
    
    Note this function returns a tuple containing (minus, plus) and that it
    expects nonnegative values. Negative values could introduce zero values if
    the zero-handling mode is 'psuedocount', and the 'floor' mode will raise
    all negative values to the floor value.
    '''
    diffs_plus = {}
    diffs_minus = {}
        
    for chrom in data_a_plus:
        plus_a_vals = data_a_plus[chrom][:,1]
        plus_b_vals = data_b_plus[chrom][:,1]
        
        minus_a_vals = data_a_minus[chrom][:,1]
        minus_b_vals = data_b_minus[chrom][:,1]
        
        if (zeroHandlingMode == "floor"):
            plus_a_vals[plus_a_vals <= 0] = zeroHandlingValue
            plus_b_vals[plus_b_vals <= 0] = zeroHandlingValue
            minus_a_vals[minus_a_vals <= 0] = zeroHandlingValue
            minus_b_vals[minus_b_vals <= 0] = zeroHandlingValue
        elif (zeroHandlingMode == "pseudocount"):
            plus_a_vals += zeroHandlingValue
            plus_b_vals += zeroHandlingValue
            minus_a_vals += zeroHandlingValue
            minus_b_vals += zeroHandlingValue
        
        diffs_plus[chrom] = np.column_stack(
            (data_a_plus[chrom][:,0], np.log2(plus_a_vals / plus_b_vals))
        )
        diffs_minus[chrom] = np.column_stack(
            (data_a_minus[chrom][:,0], np.log2(minus_a_vals / minus_b_vals))
        )
    
    return (diffs_minus, diffs_plus)

def fold_ratio(
    data_a_minus,
    data_a_plus,
    data_b_minus,
    data_b_plus,
    zeroHandlingMode = "none",
    zeroHandlingValue = 0
):
    '''
    This function takes in two datasets (which MUST be reformatted and have
    the same locations) and returns a new dataset containing the ratio of
    the values in one dataset over the values in the other dataset.
    
    Arguments:
        data_a_minus: The minus strand data for one dataset, the numerator.
        data_a_plus: The plus strand data for one dataset, the numerator.
        data_b_minus: The minus strand data for the other dataset, the
            denominator.
        data_b_plus: The plus strand data for the other dataset, the
            denominator.
        zeroHandlingMode: Which method to use to handle zero values in the
            denominator. Must be 'none', 'floor', or 'pseudocount'. Defaults
            to 'none'. If 'none', may throw errors if zero values are present
            in the b dataset.
        zeroHandlingValue: The value to use in the zero-handling mode.
    
    Note this function returns a tuple containing (minus, plus) and that it
    expects nonnegative values. Negative values could introduce zero values if
    the zero-handling mode is 'psuedocount', and the 'floor' mode will raise
    all negative values to the floor value.
    '''
    diffs_plus = {}
    diffs_minus = {}
        
    for chrom in data_a_plus:
        plus_a_vals = data_a_plus[chrom][:,1]
        plus_b_vals = data_b_plus[chrom][:,1]
        
        minus_a_vals = data_a_minus[chrom][:,1]
        minus_b_vals = data_b_minus[chrom][:,1]
        
        if (zeroHandlingMode == "floor"):
            plus_b_vals[plus_b_vals <= 0] = zeroHandlingValue
            minus_b_vals[minus_b_vals <= 0] = zeroHandlingValue
        elif (zeroHandlingMode == "pseudocount"):
            plus_b_vals += zeroHandlingValue
            minus_b_vals += zeroHandlingValue
        
        diffs_plus[chrom] = np.column_stack(
            (data_a_plus[chrom][:,0], plus_a_vals / plus_b_vals)
        )
        diffs_minus[chrom] = np.column_stack(
            (data_a_minus[chrom][:,0], minus_a_vals / minus_b_vals)
        )
    
    return (diffs_minus, diffs_plus)

def make_comparison(
    minus_a_data,
    plus_a_data,
    minus_b_data,
    plus_b_data,
    prefix,
    out_dir,
    mode = "diff",
    zeroHandlingMode = "none",
    zeroHandlingValue = 0
):
    '''
    This function takes in two datasets and produces a comparison between them.
    
    Arguments:
        minus_a_data: The minus strand data from one dataset.
        plus_a_data: The plus strand data from one dataset.
        minus_b_data: The minus strand data from the other dataset.
        plus_b_data: The plus strand data from the other dataset.
        prefix: The prefix to add to the filenames.
        out_dir: The directory to store the results in.
        mode: 'diff', 'log2', or 'ratio'. To use the absolute difference, log2
            fold ratio, or untransformed fold ratio, respectively. Defaults to
            'diff'.
        zeroHandlingMode: Passed onto 'zeroHandlingMode' in log2_fold_ratio()
            or fold_ratio() if mode is 'log2' or 'ratio'. Defaults to 'none'.
        zeroHandlingValue: Passed onto 'zeroHandlingValue' in
            log2_fold_ratio() or fold_ratio() if mode is 'log2' or 'ratio'.
            Defaults to 'none'.
    
    Note that the out_dir directory is NOT created.
    '''
    print("~~~ Comparing Datasets ~~~")
    if (mode not in ["diff", "log2", "ratio"]):
        raise Exception("Unsupported comparison mode!")
    
    if (mode == "diff"):
        print("- Calculating difference...")
        diff_minus, diff_plus = absolute_difference(
            minus_a_data,
            plus_a_data,
            minus_b_data,
            plus_b_data
        )
    elif (mode == "log2"):
        print("- Calculating log2 fold ratio...")
        diff_minus, diff_plus = log2_fold_ratio(
            minus_a_data,
            plus_a_data,
            minus_b_data,
            plus_b_data,
            zeroHandlingMode,
            zeroHandlingValue
        )
    elif (mode == 'ratio'):
        print("- Calculating fold ratio...")
        diff_minus, diff_plus = fold_ratio(
            minus_a_data,
            plus_a_data,
            minus_b_data,
            plus_b_data,
            zeroHandlingMode,
            zeroHandlingValue
        )
    print("- Outputting...")
    tsv_dir = out_dir + "/" + prefix + "_plus"
    if (not os.path.isdir(tsv_dir)):
        os.mkdir(tsv_dir)
    output_tsv(diff_plus, tsv_dir)
    output_wig(diff_plus, out_dir + "/" + prefix + "_plus.wig")
    
    tsv_dir = out_dir + "/" + prefix + "_minus"
    if (not os.path.isdir(tsv_dir)):
        os.mkdir(tsv_dir)
    output_tsv(diff_minus, tsv_dir)
    output_wig(diff_minus, out_dir + "/" + prefix + "_minus.wig")
    
    print("~~~ Compared Datasets ~~~")
    return (diff_minus, diff_plus)