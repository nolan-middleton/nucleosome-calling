##############################################################################
# This file contains library functions for handling nucleosome locations and #
# damage patterns within CPD-seq data.                                       #
##############################################################################

#%% Setup

# Imports
import numpy as np
import Analysis as an
import DataIO as io
import statsmodels.api as sm
import scipy as sp

#%% Obtaining Nucleosome Damage Patterns

def get_nucleosome_strand_pattern(nucPositions, data):
    '''
    This function will produce the nucleosome damage pattern from a given set
    of data for one strand. Will return (values, selected, pattern), where
    values contains all of the data values, selected contains zeros and ones
    indicating whether there actually was data at that position. pattern
    has four rows: mean, var, and n where each column corresponds to a
    different position along the nucleosome.
    
    Arguments:
        nucPositions: A dictionary of nucleosome positions
        data: The CPD data for one strand.
    '''
    print("~~~ Getting Nucleosome Strand Pattern ~~~")
    if (len(nucPositions.keys()) == 0):
        raise Exception("You must supply some nucleosomes!")
    cnums = np.shape(nucPositions[list(nucPositions.keys())[0]])[1]
    for chrom in nucPositions:
        if (np.shape(nucPositions[chrom])[1] != cnums):
            raise Exception(
                "Nucleosome positions must have the same number of positions!"
            )
    
    print("- Finding data points...")
    selected = {}
    values = {}
    for chrom in nucPositions:
        print("-- " + chrom + "...")
        thesePositions = nucPositions[chrom].flatten()
        shape = np.shape(nucPositions[chrom])
        theseValues = np.zeros(shape).flatten()
        theseSelected = np.isin(thesePositions, data[chrom][:,0])
        
        theseDataSelected = np.isin(data[chrom][:,0], thesePositions)
        selectedPositions = data[chrom][theseDataSelected,0]
        selectedValues = data[chrom][theseDataSelected,1]
        
        for i in range(np.sum(theseDataSelected)):
            selectedPosition=selectedPositions[i]
            theseValues[thesePositions==selectedPosition] = selectedValues[i]
        
        values[chrom] = np.reshape(theseValues, shape)
        selected[chrom] = np.reshape(theseSelected, shape)
    
    print("- Collecting data...")
    table = np.zeros((0,cnums))
    selected_table = np.zeros((0,cnums), dtype = bool)
    for chrom in nucPositions:
        table = np.vstack((table, values[chrom]))
        selected_table = np.vstack((selected_table, selected[chrom]))
    
    print("- Calculating statistics...")
    n = np.sum(selected_table, axis = 0)
    mean = np.sum(table, axis = 0) / n
    var = np.sum(
        (
            table - np.repeat(
                mean[np.newaxis,:],
                np.shape(table)[0],
                axis = 0
            )*selected_table
        )**2,
        axis = 0
    ) / (n - 1)
    
    print("~~~ Got Nucleosome Strand Pattern ~~~")
    return (values, selected, np.vstack((mean, var, n)))

def get_nucleosome_pattern(nucPositions, data_minus, data_plus):
    '''
    This function will produce the nucleosome damage pattern from a given set
    of data. Will return a tuple of values as follows: (values_minus,
    values_plus,selected_minus,selected_plus,pattern_minus,pattern_plus),
    where the pattern is an array with 4 rows: mean, var, and n, where each
    column corresponds to a different position along the nucleosome.
    
    Arguments:
        nucPositions: A dictionary of exact nucleosome positions.
        data_minus: The minus-stranded CPD data.
        data_plus: The plus-stranded CPD data.
    '''
    # Stranded Patterns
    plus_results = get_nucleosome_strand_pattern(nucPositions, data_plus)
    minus_results = get_nucleosome_strand_pattern(nucPositions, data_minus)
    
    return (
        minus_results[0],
        plus_results[0],
        minus_results[1],
        plus_results[1],
        minus_results[2],
        plus_results[2]
    )

def make_nucleosome_pattern(
    nucleosomes,
    data_minus,
    data_plus,
    out_dir,
    N = 20,
    m = -73,
    M = 73
):
    '''
    This function will take in a dictionary of nucleosome positions and CPD
    data and will output the nucleosome pattern to a directory and return the
    damage pattern. Will NOT create out_dir. Will return a tuple of values as
    follows: (values_minus, values_plus, selected_minus, selected_plus,
    pattern), where pattern is a numpy array with 9 rows: the detrended mean,
    composite trended mean, composite variance, minus mean, plus mean, minus
    variance, plus variance, minus n, plus n.
    
    Arguments:
        nucleosomes: The dictionary of nucleosome positions. Must contain
            arrays with three columns: start, stop, dyad.
        data_minus: The dictionary of minus-stranded CPD data.
        data_plus: The dictionary of plus-stranded CPD data.
        out_dir: The directory to push the outputs into.
        N: The number of base pairs for the smoothing window. Defaults to 20.
    '''
    print("~~~ Making Nucleosome Pattern ~~~")
    
    print("- Getting nucleosome positions...")
    nucleosome_positions = an.get_inbetween_positions(
        an.make_type(an.column_join(an.column_split(nucleosomes)[0:2]), int)
    )

    print("- Finding damage pattern...")
    results = get_nucleosome_pattern(nucleosome_positions,data_minus,data_plus)
    minus_vals, plus_vals, minus_selected, plus_selected = results[0:4]
    minus_pattern, plus_pattern = results[4:]
    
    minus_means = minus_pattern[0,:][::-1]
    plus_means = plus_pattern[0,:]
    minus_vars = minus_pattern[1,:][::-1]
    plus_vars = plus_pattern[1,:]
    minus_n = minus_pattern[2,:][::-1]
    plus_n = plus_pattern[2,:]
    
    all_mean = (minus_means*minus_n + plus_means*plus_n) / (minus_n + plus_n)
    X = np.arange(-73, 73) + 0.5
    trend = sm.nonparametric.lowess(all_mean, X, N/len(X))
    if (np.shape(trend)[0] != len(X)):
        raise Exception("LOWESS dropped some values!")
    flat_mean = all_mean - trend[:,1]
    
    all_var = np.zeros(len(X))
    for chrom in nucleosome_positions:
        for i in range(np.shape(nucleosome_positions[chrom])[0]):
            S = minus_selected[chrom][i,:]
            all_var[::-1][S] += (minus_vals[chrom][i,S] - all_mean[S])**2
            
            S = plus_selected[chrom][i,:]
            all_var[S] += (plus_vals[chrom][i,S] - all_mean[S])**2
    all_var /= minus_n + plus_n
    
    master_pattern = np.vstack(
        (
            flat_mean,
            all_mean,
            all_var,
            minus_means[::-1], # We'll re-reverse the minus stranded data
            plus_means,
            minus_vars[::-1],
            plus_vars,
            minus_n[::-1],
            plus_n
        )
    )

    print("- Outputting results...")
    io.output_tsv(minus_vals, out_dir + "/minus_vals")
    io.output_tsv(plus_vals, out_dir + "/plus_vals")

    io.output_tsv(minus_selected, out_dir + "/minus_selected", fmt = "%d")
    io.output_tsv(plus_selected, out_dir + "/plus_selected", fmt = "%d")

    np.savetxt(
        out_dir + "/damagePattern.tsv",
        master_pattern,
        delimiter = "\t"
    )
    
    print("~~~ Made Nucleosome Pattern ~~~")
    return (minus_vals,plus_vals,minus_selected,plus_selected,master_pattern)

def load_nucleosome_pattern_dataset(directory, chromosomes):
    '''
    This function will load in all of the data associated with the nucleosome
    pattern, including the pattern array itself, the values dictionaries, and
    the dictionaries of selected locations. pattern is a numpy array with 8
    rows: the detrended mean, composite trended mean, minus mean, plus mean,
    minus variance, plus variance, minus n, plus n
    
    Arguments:
        directory: The directory that contains the damage pattern data. Must
            contain the damagePattern_minus.tsv and damagePattern_plus.tsv
            files and the minus_vals, plus_vals, minus_selected, and
            plus_selected directories.
        chromosomes: The chromosomes of the nucleosomes used to generate the
            pattern dataset. Used in io.load_tsv for the values and selected
            dictionaries.
    
    Note that this function returns a tuple containing (pattern, values_minus,
    values_plus, selected_minus, selected_plus).
    '''
    pattern = np.loadtxt(
        directory + "/damagePattern.tsv",
        delimiter = "\t"
    )
    
    values_minus = io.load_tsv(directory + "/minus_vals", chromosomes)
    values_plus = io.load_tsv(directory + "/plus_vals", chromosomes)
    selected_minus = an.make_type(
        io.load_tsv(directory + "/minus_selected", chromosomes),
        bool
    )
    selected_plus = an.make_type(
        io.load_tsv(directory + "/plus_selected", chromosomes),
        bool
    )
    
    return (
        pattern,
        values_minus,
        values_plus,
        selected_minus,
        selected_plus
    )

def get_N():
    '''
    This function will find the number of random datasets compiled, stored in
    the file RandomPattern/datasets.txt as the only data in the file.
    '''
    with open("RandomPattern/datasets.txt") as file:
        for line in file:
            N = int(line)
    
    return N

#%% Statistics

def normal_log_likelihood(
    minus_values,
    plus_values,
    minus_dipys,
    plus_dipys,
    means,
    stds
):
    '''
    This function will take in a set of data for a stretch of 146 inbetween
    positions and a set of means and standard deviations and return the ln of
    the likelihood under the model N(x | mean, std). All these lists must be
    the same length.
    
    Arguments:
        minus_values: A list of the data values for the minus strand for one
            stretch of 146 inbetween positions.
        plus_values: A list of the data values for the plus strand for one
            stretch of 146 inbetween positions.
        minus_dipys: A boolean array for which values in minus_values to use.
        plus_dipys: A boolean for which values in plus_values to use.
        means: A list of means for each of the data positions.
        stds: A list of standard deviations for each of the data positions.
    '''
    total = np.sum(
        an.log_normal(
            minus_values[minus_dipys],
            means[::-1][minus_dipys],
            stds[::-1][minus_dipys]
        )
    )
    total += np.sum(
        an.log_normal(
            plus_values[plus_dipys],
            means[plus_dipys],
            stds[plus_dipys]
        )
    )
    return total

def normal_pooled_log_likelihood(
    minus_values,
    plus_values,
    minus_dipys,
    plus_dipys,
    means,
    stds
):
    '''
    This function will take in a set of data for many stretches of 146
    inbetween positions and a set of means and standard deviations and return
    the ln of the likelihood under the model N(x | mean, std). All these lists
    must be the same length.
    
    Arguments:
        minus_values: A table of all the data values for the minus strand,
            where each row corresponds to a different stretch of 146
            inbetween positions.
        plus_values: A table of all the data values for the plus strand,
            where each row corresponds to a different stretch of 146
            inbetween positions.
        minus_dipys: A boolean array for which values in minus_values to use.
        plus_dipys: A boolean for which values in plus_values to use.
        means: A list of means for each of the data positions.
        stds: A list of standard deviations for each of the data positions.
    '''
    return np.sum(
        normal_log_likelihood(
            minus_values.flatten(),
            plus_values.flatten(),
            minus_dipys.flatten(),
            plus_dipys.flatten(),
            np.tile(means, np.shape(minus_values)[0]),
            np.tile(stds, np.shape(minus_values)[0])
        )
    )

#%% Maximum Likelihood Estimation

def nucleosome_flat_mean(A):
    '''
    This function returns the means for a particular value of a under the model
    that the underlying mean ~ a + offsets
    
    Arguments:
        A: A list containing [a, offsets, X]
    '''
    return A[0] + A[1]

def nucleosome_linear_mean(A):
    '''
    This function returns the means for a particular value of a under the model
    that the underlying mean ~ a + b*X + offsets
    
    Arguments:
        A: A list containing [a, b, offsets, X]
    '''
    return A[0] + A[1]*A[3] + A[2]

def nucleosome_quadratic_mean(A):
    '''
    This function returns the means for a particular value of a under the model
    that the underlying mean ~ a + p*(X - h)**2 + offsets
    
    Arguments:
        A: A list containing [a, p, h, offsets, X]
    '''
    return A[0] + A[1]*(A[4] - A[2])**2 + A[3]

def pooled_objective(
    A,
    B,
    minus_values,
    plus_values,
    minus_dipys,
    plus_dipys,
    mean_fun,
    mean_params,
    std_fun,
    std_params
):
    '''
    This function will return the -ln of the pooled likelihood for a normal
    model. This function is used in scipy.optimize.minimize to find MLEs.
    
    Arguments:
        A: The list of parameters to minimize.
        B: The number of entries in A that belong to the mean_fun. A[:B] is
            passed through to mean_fun, while A[B:] is passed onto std_fun.
        minus_values: A table of minus CPD values, where each row is another
            stretch of 146 inbetween positions.
        plus_values: A table of plus CPD values, where each row is another
            stretch of 146 inbetween positions.
        minus_dipys: A table of boolean values corresponding to which values
            in minus_values are actually valid.
        plus_dipys: A table of boolean values corresponding to which values
            in plus_values are actually valid.
        mean_fun: The function that is used to determine the mean at each
            position in the stretch of 146 inbetween positions. Must take in
            only one argument, which is a list. The first B entries in that
            list must be the parameters you want to minimize, while any
            arguments after that are supplied by mean_params. The composite
            list A[:B] + mean_params is passed on to mean_fun.
        mean_params: The list of the other entries in the list of arguments
            for mean_fun, after the first B arguments supplied by A.
        std_fun: The function that is used to determine the standard deviation
            at each position in the stretch of 146 inbetween positions based
            on the means. Must take in two arguments. The first must be a
            list. The first len(A)-B entries in that list must be the
            parameters you want to minimize supplied by the elements of A
            after B, while any arguments after that are supplied by
            std_params. The composite list A[B:] + std_params is passed on to
            std_fun as the first argument. The second argument of std_fun must
            be the means.
        std_params: The list of the other entries in the list of parameters
            for the first argument of std_fun, after the first len(A)-B
            parameters supplied by A[B:].
    '''
    means = mean_fun(list(A[:B]) + list(mean_params))
    stds = std_fun(list(A[B:]) + list(std_params), means)
    
    return -normal_pooled_log_likelihood(
        minus_values,
        plus_values,
        minus_dipys,
        plus_dipys,
        means,
        stds
    )

def pooled_constraint(A, B, mean_fun, mean_params, std_fun, std_params):
    '''
    This function will return the minimum of the standard deviations
    calculated for a normal model. This is used as an inequality constraint
    for scipy.optimize.minimize, since the standard deviations must all be
    positive.
    
    Arguments:
        A: The list of parameters to minimize.
        B: The number of entries in A that belong to the mean_fun. A[:B] is
            passed through to mean_fun, while A[B:] is passed onto std_fun.
        mean_fun: The function that is used to determine the mean at each
            position in the stretch of 146 inbetween positions. Must take in
            only one argument, which is a list. The first B entries in that
            list must be the parameters you want to minimize, while any
            arguments after that are supplied by mean_params. The composite
            list A[:B] + mean_params is passed on to mean_fun.
        mean_params: The list of the other entries in the list of arguments
            for mean_fun, after the first B arguments supplied by A.
        std_fun: The function that is used to determine the standard deviation
            at each position in the stretch of 146 inbetween positions based
            on the means. Must take in two arguments. The first must be a
            list. The first len(A)-B entries in that list must be the
            parameters you want to minimize supplied by the elements of A
            after B, while any arguments after that are supplied by
            std_params. The composite list A[B:] + std_params is passed on to
            std_fun as the first argument. The second argument of std_fun must
            be the means.
        std_params: The list of the other entries in the list of parameters
            for the first argument of std_fun, after the first len(A)-B
            parameters supplied by A[B:].
    '''
    means = mean_fun(list(A[:B]) + list(mean_params))
    stds = std_fun(list(A[B:]) + list(std_params), means)
    
    return np.min(stds)

def pooled_MLE(
    initial_guess,
    B,
    minus_values,
    plus_values,
    minus_dipys,
    plus_dipys,
    mean_fun,
    mean_params,
    std_fun,
    std_params = [],
    additional_constraints = [],
    options = {},
    method = "COBYLA",
    verbose = False
):
    '''
    This function performs the minimization for finding maximum likelihood
    estimates for a normal model under the pooled dataset. Uses the COBYLA
    method for scipy.optimize.minimize by default. Refer to that documentation
    also.
    
    Arguments:
        initial_guess: The initial guess for the parameters you want to
            minimize.
        B: The number of entries in A that belong to the mean_fun. A[:B] is
            passed through to mean_fun, while A[B:] is passed onto std_fun.
        minus_values: A table of minus CPD values, where each row is another
            stretch of 146 inbetween positions.
        plus_values: A table of plus CPD values, where each row is another
            stretch of 146 inbetween positions.
        minus_dipys: A table of boolean values corresponding to which values
            in minus_values are actually valid.
        plus_dipys: A table of boolean values corresponding to which values
            in plus_values are actually valid.
        mean_fun: The function that is used to determine the mean at each
            position in the stretch of 146 inbetween positions. Must take in
            only one argument, which is a list. The first B entries in that
            list must be the parameters you want to minimize, while any
            arguments after that are supplied by mean_params. The composite
            list A[:B] + mean_params is passed on to mean_fun.
        mean_params: The list of the other entries in the list of arguments
            for mean_fun, after the first B arguments supplied by A.
        std_fun: The function that is used to determine the standard
            deviation at each position in the stretch of 146 inbetween
            positions based on the means. Must take in two arguments. The first
            must be a list. The first len(A)-B entries in that list must be the
            parameters you want to minimize supplied by the elements of A
            after B, while any arguments after that are supplied by
            std_params. The composite list A[B:] + std_params is passed on to
            std_fun as the first argument. The second argument of std_fun must
            be the means.
        std_params: The list of the other entries in the list of parameters
            for the first argument of std_fun, after the first len(A)-B
            parameters supplied by A[B:]. Defaults to [].
        additional_constraints: Extra constraints to impose on the minimization
            minimization. Read the scipy.optimize.minimize documentation for
            the format. Essentially a list of dictionaries.
        options: A dictionary of options to pass through to
            scipy.optimize.minimize. Defaults to {}.
        method: The minimization method to use. Defaults to COBYLA.
        verbose: If True, will return the whole OptimizeResult object from
            scipy.optimize.minimize. If not, will return just the list of
            parameters. Defaults to False.
    '''
    results = sp.optimize.minimize(
        fun = pooled_objective,
        x0 = initial_guess,
        args = (
            B,
            minus_values,
            plus_values,
            minus_dipys,
            plus_dipys,
            mean_fun,
            mean_params,
            std_fun,
            std_params
        ),
        method = method,
        constraints = [
            {
                "type": "ineq",
                "fun": pooled_constraint,
                "args": (
                    B,
                    mean_fun,
                    mean_params,
                    std_fun,
                    std_params
                )
            }
        ] + additional_constraints,
        options = options
    )
    
    if (not results.success):
        print("/!\\ Minimization did not succeed!")
    
    if (verbose):
        return results
    else:
        return results.x

def single_objective(
    A,
    minus_values,
    plus_values,
    minus_dipys,
    plus_dipys,
    mean_fun,
    mean_params,
    std_fun,
    std_params
):
    '''
    This function will return the -ln of the pooled likelihood for a normal
    model. This function is used in scipy.optimize.minimize to find MLEs.
    
    Arguments:
        A: The list of parameters to minimize.
        minus_values: An array of minus CPD values for one stretch of 146
            inbetween positions.
        plus_values: An array of minus CPD values for one stretch of 146
            inbetween positions.
        minus_dipys: An array of boolean values corresponding to which values
            in minus_values are actually valid.
        plus_dipys: An array of boolean values corresponding to which values
            in plus_values are actually valid.
        mean_fun: The function that is used to determine the mean at each
            position in the stretch of 146 inbetween positions. Must take in
            only one argument, which is a list. The first entry in that list
            must be the parameter you want to minimize, while any arguments
            after that are supplied by mean_params. The composite list
            [A[0]] + mean_params is passed on to mean_fun.
        mean_params: The list of the other entries in the list of arguments
            for mean_fun, after the first argument supplied by A.
        std_fun: The function that is used to determine the standard deviation
            at each position in the stretch of 146 inbetween positions based
            on the means. Must take in two arguments. The first must be a
            list. The second argument of std_fun must be the means.
        std_params: The list of parameters for the first argument of std_fun.
    '''
    means = mean_fun([A[0]] + list(mean_params))
    stds = std_fun(std_params, means)
    
    return -normal_log_likelihood(
        minus_values,
        plus_values,
        minus_dipys,
        plus_dipys,
        means,
        stds
    )

def single_constraint(A, mean_fun, mean_params, std_fun, std_params):
    '''
    This function will return the minimum of the standard deviations
    calculated for a normal model. This is used as an inequality constraint
    for scipy.optimize.minimize, since the standard deviations must all be
    positive.
    
    Arguments:
        A: The list of parameters to minimize.
        minus_values: An array of minus CPD values for one stretch of 146
            inbetween positions.
        plus_values: An array of minus CPD values for one stretch of 146
            inbetween positions.
        minus_dipys: An array of boolean values corresponding to which values
            in minus_values are actually valid.
        plus_dipys: An array of boolean values corresponding to which values
            in plus_values are actually valid.
        mean_fun: The function that is used to determine the mean at each
            position in the stretch of 146 inbetween positions. Must take in
            only one argument, which is a list. The first entry in that list
            must be the parameter you want to minimize, while any arguments
            after that are supplied by mean_params. The composite list
            [A[0]] + mean_params is passed on to mean_fun.
        mean_params: The list of the other entries in the list of arguments
            for mean_fun, after the first argument supplied by A.
        std_fun: The function that is used to determine the standard deviation
            at each position in the stretch of 146 inbetween positions based
            on the means. Must take in two arguments. The first must be a
            list. The second argument of std_fun must be the means.
        std_params: The list of parameters for the first argument of std_fun.
    '''
    means = mean_fun([A[0]] + list(mean_params))
    stds = std_fun(std_params, means)
    
    return np.min(stds)

def single_MLE(
    initial_guess,
    minus_values,
    plus_values,
    minus_dipys,
    plus_dipys,
    mean_fun,
    mean_params,
    std_fun,
    std_params,
    additional_constraints = [],
    options = {},
    verbose = False
):
    '''
    This function performs the minimization for finding maximum likelihood
    estimates for a normal model for a single stretch of 146 inbetween
    positions. Uses the COBYLA method for scipy.optimize.minimize. Refer to
    that documentation, also.
    
    Arguments:
        initial_guess: The initial guess for the parameter you want to
            minimize.
        minus_values: An array of minus CPD values for one stretch of 146
            inbetween positions.
        plus_values: An array of minus CPD values for one stretch of 146
            inbetween positions.
        minus_dipys: An array of boolean values corresponding to which values
            in minus_values are actually valid.
        plus_dipys: An array of boolean values corresponding to which values
            in plus_values are actually valid.
        mean_fun: The function that is used to determine the mean at each
            position in the stretch of 146 inbetween positions. Must take in
            only one argument, which is a list. The first entry in that list
            must be the parameter you want to minimize, while any arguments
            after that are supplied by mean_params. The composite list
            [a] + mean_params is passed on to mean_fun.
        mean_params: The list of the other entries in the list of arguments
            for mean_fun, after the first argument, which is a.
        std_fun: The function that is used to determine the standard deviation
            at each position in the stretch of 146 inbetween positions based
            on the means. Must take in two arguments. The first must be a
            list. The second argument of std_fun must be the means.
        std_params: The list of parameters for the first argument of std_fun.
        additional_constraints: Extra constraints to impose on the COBYLA
            minimization. Read the scipy.optimize.minimize documentation for
            the format. Essentially a list of dictionaries. This uses COBYLA,
            so all constraints must be inequality constraints.
        options: A dictionary of options to pass on to scipy.optimize.minimize.
            Defaults to {}.
        verbose: If True, will return the whole OptimizeResult object from
            scipy.optimize.minimize. If not, will return just the list of
            parameters. Defaults to False.
    '''
    results = sp.optimize.minimize(
        fun = single_objective,
        x0 = [initial_guess],
        args = (
            minus_values,
            plus_values,
            minus_dipys,
            plus_dipys,
            mean_fun,
            mean_params,
            std_fun,
            std_params
        ),
        method = "COBYLA",
        constraints = [
            {
                "type": "ineq",
                "fun": single_constraint,
                "args": (mean_fun, mean_params, std_fun, std_params)
            }
        ] + additional_constraints,
        options = options
    )
    
    if (not results.success):
        print("/!\\ Minimization did not succeed!")
    
    if (verbose):
        return results
    else:
        return results.x[0]

def prior_objective(M, a_hats):
    '''
    This function will return the -ln of the prior distribution (being
    treated as a "likelihood" of sorts) to estimate the hyperparameters. Used
    in scipy.optimize.minimize.
    
    Arguments:
        M: A list of parameters you want to minimze. M[0] must be mu, the mean,
            while M[1] must be tau, the standard deviation.
        a_hats: The MLE estimates for the parameter you want to estimate the
            prior distribution for.
    '''
    total = 0
    for chrom in a_hats:
        total -= np.sum(an.log_normal(a_hats[chrom], M[0], M[1]))
    
    return total

def hyperparameter_MLE(
    initial_guess,
    minus_values,
    plus_values,
    minus_dipys,
    plus_dipys,
    mean_fun,
    mean_params,
    std_fun,
    std_params,
    additional_constraints = [],
    cobyla_options = {},
    powell_options = {},
    verbose = False
):
    '''
    This function will find the MLE for the hyperparameters based on the MLE
    for the parameter for every data point. Uses the Powell method for
    scipy.optimize.minimize. Read that documentation, also.
    
    Arguments:
        initial_guess: The initial guess for the value of the parameter you
            want to find the hyperparameters for. The MLE usually works pretty
            good here.
        minus_values: A dictionary of arrays of minus-stranded CPD data, where
            each row in the array is a different stretch of 146 inbetween
            positions. Each column then corresponds to a different position
            along the stretch. For each key k in minus_values, row j in the
            array paired with k must represent the same stretch as row j in
            the array paired with key k in plus_values and in minus_dipys.
        plus_values: A dictionary of arrays of plus-stranded CPD data, where
            each row in the array is a different stretch of 146 inbetween
            positions. Each column then corresponds to a different position
            along the stretch. For each key k in plus_values, row j in the
            array paired with k must represent the same stretch as row j in
            the array paired with key k in minus_values and in plus_dipys.
        minus_dipys: A dictionary of boolean arrays, where each row in the
            array is a different stretch of 146 inbetween positions. The
            entries in these array must correspond to whether or not the CPD
            data in the index in minus_values is valid (i.e. there is actually
            a dipyrimidine at that location).
        plus_dipys: A dictionary of boolean arrays, where each row in the
            array is a different stretch of 146 inbetween positions. The
            entries in these array must correspond to whether or not the CPD
            data in the index in plus_values is valid (i.e. there is actually
            a dipyrimidine at that location).
        mean_fun: The function that is used to determine the mean at each
            position in the stretch of 146 inbetween positions. Must take in
            only one argument, which is a list. The first entry in that list
            must be the parameter you want to minimize, while any arguments
            after that are supplied by mean_params. The composite list
            [a] + mean_params is passed on to mean_fun.
        mean_params: The list of the other entries in the list of arguments
            for mean_fun, after the first argument, which is a.
        std_fun: The function that is used to determine the standard deviation
            at each position in the stretch of 146 inbetween positions based
            on the means. Must take in two arguments. The first must be a
            list. The second argument of std_fun must be the means.
        std_params: The list of parameters for the first argument of std_fun.
        additional_constraints: Extra constraints to impose on the COBYLA
            minimization for the individual MLE estimates. Read the
            scipy.optimize.minimize documentation for the format. Essentially
            a list of dictionaries. This will be passed on to single_MLE, which
            uses COBYLA, so all constraints must be inequality constraints.
        cobyla_options: A dictionary of options to pass on to the first
            scipy.optimize.minimize for the COBYLA minimization. Defaults to
            {}.
        powell_options: A dictionary of options to pass on to the second
            scipy.optimize.minimize for the Powell minimization. Defaults to
            {}.
        verbose: If True, will return the whole OptimizeResult object from
            scipy.optimize.minimize. If not, will return just a dictionary of
            parameters. Defaults to False.
    
    Note that this function returns a dictionary of results. The key "mu" is
    the mean, the key "tau" is the standard deviation, and the key "a_hat"
    contains all of the MLEs for the a parameter at every position. If verbose
    is set to True, it will instead return the entire scipy.optimize.minimize
    object for the hyperparameter MLE and the list of a MLEs.
    '''
    print("~~~ Finding Hyperparameter MLEs ~~~")
    a_hats = {}
    for chrom in minus_values:
        print("- " + chrom + "...")
        J = np.shape(minus_values[chrom])[0]
        a_hats[chrom] = np.zeros(J)
        for j in range(J):
            a_hats[chrom][j] = single_MLE(
                initial_guess,
                minus_values[chrom][j,:],
                plus_values[chrom][j,:],
                minus_dipys[chrom][j,:],
                plus_dipys[chrom][j,:],
                mean_fun,
                mean_params,
                std_fun,
                std_params,
                additional_constraints,
                cobyla_options
            )

    hyper_results = sp.optimize.minimize(
        fun = prior_objective,
        x0 = [initial_guess, 1],
        args = a_hats,
        method = "Powell",
        bounds = [(-np.inf, np.inf), (0, np.inf)],
        options = powell_options
    )
    
    if (not hyper_results.success):
        print("/!\\ Minimization did not succeed!")
    
    print("~~~ Found Hyperparameter MLEs ~~~")
    if (not verbose):
        return {
            "mu": hyper_results.x[0],
            "tau": hyper_results.x[1],
            "a_hat": a_hats
        }
    else:
        return (hyper_results, a_hats)

def output_MLEs(out_name, pooled_MLE, hyper_MLE):
    '''
    This function will output all of the MLEs for a specific model in a single,
    combined dictionary.
    
    Arguments:
        out_name: The name of the output file to save the data to.
        pooled_MLE: The results of the pooled MLEs.
        hyper_MLE: The results of the hyperparameter MLEs.
    '''
    output_a_hats = {}
    for chrom in hyper_MLE["a_hat"]:
        output_a_hats[chrom] = list(hyper_MLE["a_hat"][chrom])
    
    out_dict = {}
    for key in pooled_MLE:
        out_dict[key] = pooled_MLE[key]
    out_dict["mu"] = hyper_MLE["mu"]
    out_dict["tau"] = hyper_MLE["tau"]
    out_dict["a_hat"] = output_a_hats

    io.output_json(out_name, out_dict)

#%% Integration

def nucleosome_flat_model(a, MLEs, offsets, X):
    '''
    This function returns the expected values for the nucleosomes under the
    flat model.
    
    Arguments:
        a: The a parameter, generally being integrated over.
        MLEs: The MLEs for the other parameters. Unused in this function.
        offsets: The nucleosome offsets.
        X: The X-values. Should be -72.5 to 72.5
    '''
    return a + offsets

def nucleosome_linear_model(a, MLEs, offsets, X):
    '''
    This function returns the expected values for the nucleosomes under the
    linear model.
    
    Arguments:
        a: The a parameter, generally being integrated over.
        MLEs: The MLEs for the other parameters. Must have the "b" key.
        offsets: The nucleosome offsets.
        X: The X-values. Should be -72.5 to 72.5
    '''
    return a + MLEs["b"]*X + offsets

def nucleosome_quadratic_model(a, MLEs, offsets, X):
    '''
    This function returns the expected values for the nucleosomes under the
    quadratic model.
    
    Arguments:
        a: The a parameter, generally being integrated over.
        MLEs: The MLEs for the other parameters. Must have "p" and "h" keys.
        offsets: The nucleosome offsets.
        X: The X-values. Should be -72.5 to 72.5
    '''
    return a + MLEs["p"]*(X - MLEs["h"])**2 + offsets

def nucleosome_bounds(
    MLEs,
    mean_fun,
    offsets,
    X,
    std_relation,
    granularity = 1000,
    R = 4
):
    '''
    This function will search for valid bounds for the nucleosome model.
    
    Arguments:
        MLEs: A dictionary of MLEs for the nucleosome model. Must contain the
            "mu" and "tau" keys plus whatever keys are needed by mean_fun and
            std_relation.
        mean_fun: A function to determine the means. Must take a as the first
            argument, MLEs as the second argument, offsets as the third
            argument, and X as the fourth and final argument.
        offsets: The periodic nucleosome offsets.
        X: The inbetween positions.
        std_relation: The function to derive the standard deviations from the
            means. Must take MLEs as the first argument and the means as the
            second argument.
        granularity: The number of divisions to break the search area up into.
            Defaults to 1000.
        R: The number of standard deviations to integrate over. Defaults to 4.
    '''
    A = np.linspace(
        MLEs["mu"] - R*MLEs["tau"],
        MLEs["mu"] + R*MLEs["tau"],
        granularity
    )
    
    means = [mean_fun(a, MLEs, offsets, X) for a in A]
    stds = [std_relation(MLEs, mean) for mean in means]
    
    valid = np.array([False not in (std > 0) for std in stds])
    
    if (True not in valid):
        raise Exception("Bounds totally invalid!")
    
    if (False in valid):
        S = np.split(valid, np.arange(granularity)[~valid])
        I = np.split(np.arange(granularity), np.arange(granularity)[~valid])
        I = [I[i][S[i]] for i in range(len(I))]
        
        C = I[np.argmax([len(a) for a in I])]
        
        return (A[C[0]], A[C[-1]])
    else:
        return (MLEs["mu"] - R*MLEs["tau"], MLEs["mu"] + R*MLEs["tau"])

def prior_int(bounds, MLE):
    '''
    This function takes in the bounds of the integral and the MLEs and returns
    the value of the integral of the prior function over integration bounds.
    Necessary for normalizing the probability distributions.
    
    Arguments:
        bounds: The bounds of the integration function.
        MLE: The maximum likelihood estimates. Must contain the "mu" and "tau"
            keys.
    
    Note that this function returns a tuple: (integral value, log_factor)
    '''
    return an.scaled_integral(
        lambda a, mu, tau : -an.log_normal(a, mu, tau),
        lambda a,mu,tau,log_factor:np.exp(an.log_normal(a,mu,tau)-log_factor),
        [MLE["mu"], MLE["tau"]],
        bounds
    )

def nucleosome_integrand(
    a,
    values_minus,
    values_plus,
    dipys_minus,
    dipys_plus,
    offsets,
    X,
    MLEs,
    mean_fun,
    std_relation,
    log_factor
):
    '''
    This is the integrand function for the nucleosome model. To be used with
    scipy.integrate.quad, or similar.
    
    Arguments:
        a: The a parameter of the model to integrate over.
        values_minus: The values for the minus strand of one stretch of 146
            inbetween positions.
        values_plus: The values for the plus strand of one stretch of 146
            inbetween positions.
        dipys_minus: A boolean array for the minus strand of one stretch of 146
            inbetween positions, containing True for every position that is a
            valid CPD.
        dipys_plus: A boolean array for the plus strand of one stretch of 146
            inbetween positions, containing True for every position that is a
            valid CPD.
        offsets: The nucleosomal period offsets.
        X: The inbetween position values, should be -72.5 to 72.5.
        MLEs: A dictionary of MLEs for the other parameters. Must have the "mu"
            and "tau" keys, plus any keys necessary for mean_fun and
            std_relation.
        mean_fun: The function to derive the means. Must take the a parameter
            as its first argument, the MLE dictionary as its second argument,
            the offsets as its third argument, and the X values as its fourth
            and final argument.
        std_relation: The function to derive the standard deviations from the
            means. Must take in the MLE dictionary as its first argument and
            the means as its second argument.
        log_factor: A log factor used to scale the integral and make it
            tractable.
    '''
    means = mean_fun(a, MLEs, offsets, X)
    stds = std_relation(MLEs, means)
    
    log_likelihood = normal_log_likelihood(
        values_minus,
        values_plus,
        dipys_minus,
        dipys_plus,
        means,
        stds
    )
    
    log_prior = an.log_normal(a, MLEs["mu"], MLEs["tau"])
    
    return np.exp(np.sum(log_likelihood) + np.sum(log_prior) - log_factor)

def nucleosome_integrand_objective(
    a,
    values_minus,
    values_plus,
    dipys_minus,
    dipys_plus,
    offsets,
    X,
    MLEs,
    mean_fun,
    std_relation
):
    '''
    This is the objective function used to maximize the integrand function for
    the nucleosome model. To be used with an.scaled_integral, or similar. This
    function effectively returns the -log of the integrand.
    
    Arguments:
        a: The a parameter of the model to integrate over.
        values_minus: The values for the minus strand of one stretch of 146
            inbetween positions.
        values_plus: The values for the plus strand of one stretch of 146
            inbetween positions.
        dipys_minus: A boolean array for the minus strand of one stretch of 146
            inbetween positions, containing True for every position that is a
            valid CPD.
        dipys_plus: A boolean array for the plus strand of one stretch of 146
            inbetween positions, containing True for every position that is a
            valid CPD.
        offsets: The nucleosomal period offsets.
        X: The inbetween position values, should be -72.5 to 72.5.
        MLEs: A dictionary of MLEs for the other parameters. Must have the "mu"
            and "tau" keys, plus any keys necessary for mean_fun and
            std_relation.
        mean_fun: The function to derive the means. Must take the a parameter
            as its first argument, the MLE dictionary as its second argument,
            the offsets as its third argument, and the X values as its fourth
            and final argument.
        std_relation: The function to derive the standard deviations from the
            means. Must take in the MLE dictionary as its first argument and
            the means as its second argument.
    '''
    means = mean_fun(a, MLEs, offsets, X)
    stds = std_relation(MLEs, means)
    
    log_likelihood = normal_log_likelihood(
        values_minus,
        values_plus,
        dipys_minus,
        dipys_plus,
        means,
        stds
    )
    
    log_prior = an.log_normal(a, MLEs["mu"], MLEs["tau"])
    
    return -np.sum(log_likelihood) - np.sum(log_prior)

def random_bounds(
    MLEs,
    std_relation,
    granularity = 1000,
    R = 4
):
    '''
    This function will search for valid bounds for the random model.
    
    Arguments:
        MLEs: A dictionary of MLEs for the nucleosome model. Must contain the
            "mu" and "tau" keys plus whatever keys are needed by std_relation.
        std_relation: The function to derive the standard deviations from the
            means. Must take MLEs as the first argument and the means as the
            second argument.
        granularity: The number of divisions to break the search area up into.
            Defaults to 1000.
        R: The number of standard deviations to integrate over. Defaults to 4.
    '''
    A = np.linspace(
        MLEs["mu"] - R*MLEs["tau"],
        MLEs["mu"] + R*MLEs["tau"],
        granularity
    )
    
    means = A
    stds = std_relation(MLEs, means)
    
    valid = np.array(stds > 0)
    
    if (True not in valid):
        raise Exception("Bounds totally invalid!")
    
    if (False in valid):
        S = np.split(valid, np.arange(granularity)[~valid])
        I = np.split(np.arange(granularity), np.arange(granularity)[~valid])
        I = [I[i][S[i]] for i in range(len(I))]
        
        C = I[np.argmax([len(a) for a in I])]
        
        return (A[C[0]], A[C[-1]])
    else:
        return (MLEs["mu"] - R*MLEs["tau"], MLEs["mu"] + R*MLEs["tau"])

def random_integrand(
    a,
    values_minus,
    values_plus,
    dipys_minus,
    dipys_plus,
    MLEs,
    std_relation,
    log_factor
):
    '''
    This is the integrand function for the random model. To be used with
    scipy.integrate.quad, or similar.
    
    Arguments:
        a: The a parameter of the model to integrate over.
        values_minus: The values for the minus strand of one stretch of 146
            inbetween positions.
        values_plus: The values for the plus strand of one stretch of 146
            inbetween positions.
        dipys_minus: A boolean array for the minus strand of one stretch of 146
            inbetween positions, containing True for every position that is a
            valid CPD.
        dipys_plus: A boolean array for the plus strand of one stretch of 146
            inbetween positions, containing True for every position that is a
            valid CPD.
        MLEs: A dictionary of MLEs for the other parameters. Must have the "mu"
            and "tau" keys, plus any keys necessary for std_relation.
        std_relation: The function to derive the standard deviations from the
            means. Must take in the MLE dictionary as its first argument and
            the means as its second argument.
        log_factor: A log factor used to scale the integral and make it
            tractable.
    '''
    means = a + np.zeros(len(values_minus))
    stds = std_relation(MLEs, means)
    
    log_likelihood = normal_log_likelihood(
        values_minus,
        values_plus,
        dipys_minus,
        dipys_plus,
        means,
        stds
    )
    
    log_prior = an.log_normal(a, MLEs["mu"], MLEs["tau"])
    
    return np.exp(np.sum(log_likelihood) + np.sum(log_prior) - log_factor)

def random_integrand_objective(
    a,
    values_minus,
    values_plus,
    dipys_minus,
    dipys_plus,
    MLEs,
    std_relation
):
    '''
    This is the objective function used to maximize the integrand function for
    the nucleosome model. To be used with an.scaled_integral, or similar. This
    function effectively returns the -log of the integrand.
    
    Arguments:
        a: The a parameter of the model to integrate over.
        values_minus: The values for the minus strand of one stretch of 146
            inbetween positions.
        values_plus: The values for the plus strand of one stretch of 146
            inbetween positions.
        dipys_minus: A boolean array for the minus strand of one stretch of 146
            inbetween positions, containing True for every position that is a
            valid CPD.
        dipys_plus: A boolean array for the plus strand of one stretch of 146
            inbetween positions, containing True for every position that is a
            valid CPD.
        MLEs: A dictionary of MLEs for the other parameters. Must have the "mu"
            and "tau" keys, plus any keys necessary for std_relation.
        std_relation: The function to derive the standard deviations from the
            means. Must take in the MLE dictionary as its first argument and
            the means as its second argument.
    '''
    
    means = a + np.zeros(len(values_minus))
    stds = std_relation(MLEs, means)
    
    log_likelihood = normal_log_likelihood(
        values_minus,
        values_plus,
        dipys_minus,
        dipys_plus,
        means,
        stds
    )
    
    log_prior = an.log_normal(a, MLEs["mu"], MLEs["tau"])
    
    return -np.sum(log_likelihood) - np.sum(log_prior)

def score_genome(
    genome,
    blocks,
    chroms,
    data_minus,
    data_plus,
    offsets,
    nucleosome_model,
    nucleosome_std_relation,
    nucleosome_MLE,
    random_std_relation,
    random_MLE
):
    '''
    This function will score the entire genome according to one model. Returns
    a tuple of dicts: (bayes factors, nucleosome integrals, random integrals,
    nucleosome prior normalization factor, random prior normalization factor).
    
    Arguments:
        genome: The genome array dictionary.
        blocks: The blocks dictionary detailing which blocks are invalid.
        chroms: The chromosomes of the genome to score.
        data_minus: The minus-stranded CPD data.
        data_plus: The plus-stranded CPD data.
        offsets: The periodic offsets of the nucleosome model.
        nucleosome_model: The nucleosome mean function, passed on to
            nucleosome_integrand_objective and nucleosome_integrand. Must take
            the a parameter as its first argument, the MLE dictionary for its
            second argument, the offsets as its third argument, and the X
            values as its fourth argument.
        nucleosome_std_relation: The nucleosome function to derive the
            standard deviation from the means, passed on to
            nucleosome_integrand_objective, and nucleosome_integrand. Must take
            in the MLE dictionary as its first argument and the means as its
            second argument.
        nucleosome_MLE: The nucleosome MLE values for the other parameters.
            Must have keys for all of the parameters for the nucleosome_model
            and nucleosome_std_relation and for the prior distribution (mu and
            tau).
        random_std_relation: The random function to derive the standard
            deviation from the means, passed on to random_integrand_objective
            and random_integrand. Must take in the MLE dictionary as its first
            argument and the means as its second argument.
        random_MLE: The random MLE values for the other parameters. Must have
            keys for all of the parameters for random_std_relation and for the
            prior distribution (mu and tau).
    '''
    print("~~~ Scoring Genome ~~~")
    # Initial values
    D = len(offsets) #"diameter"
    R = int(D/2) #"radius"
    X = np.arange(-R, R) + 0.5
    
    # Find bounds
    print("- Finding bounds...")
    random_int_bounds = random_bounds(random_MLE, random_std_relation)
    nucleosome_int_bounds = nucleosome_bounds(
        nucleosome_MLE,
        nucleosome_model,
        offsets,
        X,
        nucleosome_std_relation
    )
    
    # Find normalization factors
    print("- Integrating priors...")
    rand_P_int = prior_int(random_int_bounds, random_MLE)
    nuc_P_int = prior_int(nucleosome_int_bounds, nucleosome_MLE)
    
    # Main loop...
    nucleosome_results = {}
    random_results = {}
    bayes_factors = {}
    for chrom in chroms:
        print("- " + chrom + "...")
        badBlocks = blocks[chrom][blocks[chrom][:,2].astype(bool),0:2]
        
        print("-- Reformatting data...")
        positions = np.arange(0, len(genome[chrom]) - 1) + 0.5
        
        minus_values = np.zeros(len(positions))
        minus_dipys = np.zeros(len(positions), dtype = bool)
        dummy, positionIndices, dataIndices = np.intersect1d(
            positions,
            data_minus[chrom][:,0],
            return_indices = True
        )
        minus_values[positionIndices] = data_minus[chrom][dataIndices,1]
        minus_dipys[positionIndices] = True
        
        plus_values = np.zeros(len(positions))
        plus_dipys = np.zeros(len(positions), dtype = bool)
        dummy, positionIndices, dataIndices = np.intersect1d(
            positions,
            data_plus[chrom][:,0],
            return_indices = True
        )
        plus_values[positionIndices] = data_plus[chrom][dataIndices,1]
        plus_dipys[positionIndices] = True
        
        print("-- Filtering out invalid positions...")
        valid_genome_locs = np.zeros(len(genome[chrom]), dtype = bool)
        valid_genome_locs[R:-R] = True
        for i in range(np.shape(badBlocks)[0]):
            valid_genome_locs[(badBlocks[i,0]-R):(badBlocks[i,1]+R)] = False
        valid_genome_locs = np.arange(len(genome[chrom]))[valid_genome_locs]
        
        test_indices = np.arange(len(positions))[
            np.isin((positions + R - 0.5).astype(int), valid_genome_locs)
        ]
        
        print("-- Commencing integration...")
        random_results[chrom] = np.zeros((len(test_indices), 3))
        nucleosome_results[chrom] = np.zeros((len(test_indices), 3))
        
        for i in range(len(test_indices)):
            index = test_indices[i]
            dyad = int(positions[index] + R - 0.5)
            if (dyad % 10000 == 0):
                print("--- "+str(dyad)+"/"+str(len(genome[chrom]))+"...")
            
            theseMinusVals = minus_values[index:(index + D)]
            thesePlusVals = plus_values[index:(index + D)]
            theseMinusDipys = minus_dipys[index:(index + D)]
            thesePlusDipys = plus_dipys[index:(index + D)]
            
            theseRandomIntResults = an.scaled_integral(
                random_integrand_objective,
                random_integrand,
                [
                    theseMinusVals,
                    thesePlusVals,
                    theseMinusDipys,
                    thesePlusDipys,
                    random_MLE,
                    random_std_relation
                ],
                random_int_bounds
            )
            random_results[chrom][i,1]=theseRandomIntResults[0]/rand_P_int[0]
            random_results[chrom][i,2]=theseRandomIntResults[1]-rand_P_int[1]
            random_results[chrom][i,0] = dyad
            
            theseNucIntResults = an.scaled_integral(
                nucleosome_integrand_objective,
                nucleosome_integrand,
                [
                    theseMinusVals,
                    thesePlusVals,
                    theseMinusDipys,
                    thesePlusDipys,
                    offsets,
                    X,
                    nucleosome_MLE,
                    nucleosome_model,
                    nucleosome_std_relation
                ],
                nucleosome_int_bounds
            )
            nucleosome_results[chrom][i,1]=theseNucIntResults[0]/nuc_P_int[0]
            nucleosome_results[chrom][i,2]=theseNucIntResults[1]-nuc_P_int[1]
            nucleosome_results[chrom][i,0] = dyad
        
        print("-- Bayes factors...")
        bayes_factors[chrom] = np.zeros(
            (np.shape(nucleosome_results[chrom])[0], 2)
        )
        bayes_factors[chrom][:,1] = np.log(
            nucleosome_results[chrom][:,1]
        ) + nucleosome_results[chrom][:,2] - np.log(
            random_results[chrom][:,1]
        ) - random_results[chrom][:,2]
        bayes_factors[chrom][:,0] = nucleosome_results[chrom][:,0]
    
    print("~~~ Scored Genome ~~~")
    return (
        bayes_factors,
        nucleosome_results,
        random_results,
        nuc_P_int,
        rand_P_int
    )

def score_genome_shuffled(
    genome,
    blocks,
    chroms,
    data_minus,
    data_plus,
    offsets,
    nucleosome_model,
    nucleosome_std_relation,
    nucleosome_MLE,
    random_std_relation,
    random_MLE,
    rng
):
    '''
    This function will score the entire genome according to one model. Returns
    a tuple of dicts: (bayes factors, nucleosome integrals, random integrals,
    nucleosome prior normalization factor, random prior normalization factor).
    
    This function shuffles the values around first as a control.
    
    Arguments:
        genome: The genome array dictionary.
        blocks: The blocks dictionary detailing which blocks are invalid.
        chroms: The chromosomes of the genome to score.
        data_minus: The minus-stranded CPD data.
        data_plus: The plus-stranded CPD data.
        offsets: The periodic offsets of the nucleosome model.
        nucleosome_model: The nucleosome mean function, passed on to
            nucleosome_integrand_objective and nucleosome_integrand. Must take
            the a parameter as its first argument, the MLE dictionary for its
            second argument, the offsets as its third argument, and the X
            values as its fourth argument.
        nucleosome_std_relation: The nucleosome function to derive the
            standard deviation from the means, passed on to
            nucleosome_integrand_objective, and nucleosome_integrand. Must take
            in the MLE dictionary as its first argument and the means as its
            second argument.
        nucleosome_MLE: The nucleosome MLE values for the other parameters.
            Must have keys for all of the parameters for the nucleosome_model
            and nucleosome_std_relation and for the prior distribution (mu and
            tau).
        random_std_relation: The random function to derive the standard
            deviation from the means, passed on to random_integrand_objective
            and random_integrand. Must take in the MLE dictionary as its first
            argument and the means as its second argument.
        random_MLE: The random MLE values for the other parameters. Must have
            keys for all of the parameters for random_std_relation and for the
            prior distribution (mu and tau).
        rng: The random number generator object to use.
    '''
    print("~~~ Scoring Genome ~~~")
    # Initial values
    D = len(offsets) #"diameter"
    R = int(D/2) #"radius"
    X = np.arange(-R, R) + 0.5
    
    # Find bounds
    print("- Finding bounds...")
    random_int_bounds = random_bounds(random_MLE, random_std_relation)
    nucleosome_int_bounds = nucleosome_bounds(
        nucleosome_MLE,
        nucleosome_model,
        offsets,
        X,
        nucleosome_std_relation
    )
    
    # Find normalization factors
    print("- Integrating priors...")
    rand_P_int = prior_int(random_int_bounds, random_MLE)
    nuc_P_int = prior_int(nucleosome_int_bounds, nucleosome_MLE)
    
    # Main loop...
    nucleosome_results = {}
    random_results = {}
    bayes_factors = {}
    for chrom in chroms:
        print("- " + chrom + "...")
        badBlocks = blocks[chrom][blocks[chrom][:,2].astype(bool),0:2]
        
        print("-- Reformatting data...")
        positions = np.arange(0, len(genome[chrom]) - 1) + 0.5
        
        minus_values = np.zeros(len(positions))
        minus_dipys = np.zeros(len(positions), dtype = bool)
        dummy, positionIndices, dataIndices = np.intersect1d(
            positions,
            data_minus[chrom][:,0],
            return_indices = True
        )
        minus_values[positionIndices] = data_minus[chrom][dataIndices,1]
        minus_dipys[positionIndices] = True
        
        plus_values = np.zeros(len(positions))
        plus_dipys = np.zeros(len(positions), dtype = bool)
        dummy, positionIndices, dataIndices = np.intersect1d(
            positions,
            data_plus[chrom][:,0],
            return_indices = True
        )
        plus_values[positionIndices] = data_plus[chrom][dataIndices,1]
        plus_dipys[positionIndices] = True
        
        print("-- Filtering out invalid positions...")
        valid_genome_locs = np.zeros(len(genome[chrom]), dtype = bool)
        valid_genome_locs[R:-R] = True
        for i in range(np.shape(badBlocks)[0]):
            valid_genome_locs[(badBlocks[i,0]-R):(badBlocks[i,1]+R)] = False
        valid_genome_locs = np.arange(len(genome[chrom]))[valid_genome_locs]
        
        test_indices = np.arange(len(positions))[
            np.isin((positions + R - 0.5).astype(int), valid_genome_locs)
        ]
        
        print("-- Commencing integration...")
        random_results[chrom] = np.zeros((len(test_indices), 3))
        nucleosome_results[chrom] = np.zeros((len(test_indices), 3))
        
        for i in range(len(test_indices)):
            index = test_indices[i]
            dyad = int(positions[index] + R - 0.5)
            if (dyad % 10000 == 0):
                print("--- "+str(dyad)+"/"+str(len(genome[chrom]))+"...")
            
            theseMinusVals = minus_values[index:(index + D)]
            thesePlusVals = plus_values[index:(index + D)]
            theseMinusDipys = minus_dipys[index:(index + D)]
            thesePlusDipys = plus_dipys[index:(index + D)]
            
            # Shuffling
            theseMinusVals[theseMinusDipys] = rng.choice(
                theseMinusVals[theseMinusDipys],
                np.sum(theseMinusDipys),
                replace = False
            )
            thesePlusVals[thesePlusDipys] = rng.choice(
                thesePlusVals[thesePlusDipys],
                np.sum(thesePlusDipys),
                replace = False
            )
            
            theseRandomIntResults = an.scaled_integral(
                random_integrand_objective,
                random_integrand,
                [
                    theseMinusVals,
                    thesePlusVals,
                    theseMinusDipys,
                    thesePlusDipys,
                    random_MLE,
                    random_std_relation
                ],
                random_int_bounds
            )
            random_results[chrom][i,1]=theseRandomIntResults[0]/rand_P_int[0]
            random_results[chrom][i,2]=theseRandomIntResults[1]-rand_P_int[1]
            random_results[chrom][i,0] = dyad
            
            theseNucIntResults = an.scaled_integral(
                nucleosome_integrand_objective,
                nucleosome_integrand,
                [
                    theseMinusVals,
                    thesePlusVals,
                    theseMinusDipys,
                    thesePlusDipys,
                    offsets,
                    X,
                    nucleosome_MLE,
                    nucleosome_model,
                    nucleosome_std_relation
                ],
                nucleosome_int_bounds
            )
            nucleosome_results[chrom][i,1]=theseNucIntResults[0]/nuc_P_int[0]
            nucleosome_results[chrom][i,2]=theseNucIntResults[1]-nuc_P_int[1]
            nucleosome_results[chrom][i,0] = dyad
        
        print("-- Bayes factors...")
        bayes_factors[chrom] = np.zeros(
            (np.shape(nucleosome_results[chrom])[0], 2)
        )
        bayes_factors[chrom][:,1] = np.log(
            nucleosome_results[chrom][:,1]
        ) + nucleosome_results[chrom][:,2] - np.log(
            random_results[chrom][:,1]
        ) - random_results[chrom][:,2]
        bayes_factors[chrom][:,0] = nucleosome_results[chrom][:,0]
    
    print("~~~ Scored Genome ~~~")
    return (
        bayes_factors,
        nucleosome_results,
        random_results,
        nuc_P_int,
        rand_P_int
    )

#%% Mapping

def greedy_placement(
    scores,
    exclusion_radius = 117,
    threshold = 0,
    maxiter = 100000,
    verbose = True
):
    '''
    This function will take in a dictionary of scores across the genome and
    place the nucleosomes in a 'greedy algorithmic' fashion, where a nucleosome
    dyad is placed at the point with the highest score, a region of width is
    then excluded from consideration, and the process is repeated.
    
    Arguments:
        scores: The dictionary of scores. Must be a dictionary of arrays, each
            with two columns: the location and the score.
        exclusion_radius: The distance on either side of called dyad to
            exclude before another dyad can be called. A region of width
            2*exclusion_radius+1, centered on the called dyad, is then 'masked
            off.' This is basically the minimum acceptable distance between two
            adjacent dyads. Defaults to 117, resulting in a 117bp exclusion
            window around each dyad, as in the Hartemink paper (Zhong et al.
            PMID: 26772197).
        threshold: The threshold of scores to consider. Will only consider
            scores greater than threshold to put a nucleosome in. Defaults to
            -np.inf.
        maxiter: The maximum number of nucleosome placement iterations to go
            through. Safety variable. Defaults to 100000.
        verbose: If True, will print out detailed status messages. Defaults to
            True.
    '''
    print("~~~ Greedily Placing Nucleosomes ~~~")
    dyads = {}
    for chrom in scores:
        print("- " + chrom + "...")
        dyads[chrom] = np.zeros(np.shape(scores[chrom])[0], dtype = int)
        n_dyads = 0
        
        potential_dyads = scores[chrom][:,0].astype(int)
        potential_scores = scores[chrom][:,1]
        
        potential_dyads = potential_dyads[potential_scores > threshold]
        potential_scores = potential_scores[potential_scores > threshold]
        
        i = 0
        while (len(potential_dyads) > 0) and (i < maxiter):
            if (verbose) and (i % 1000 == 0):
                print("-- Iteration: " + str(i) + "...")
                print("--- Potential dyads: "+str(len(potential_dyads))+"...")
            i += 1
            
            thisDyad = potential_dyads[np.argmax(potential_scores)]
            dyads[chrom][n_dyads] = thisDyad
            n_dyads += 1
            
            still_valid = np.logical_or(
                potential_dyads < thisDyad - exclusion_radius,
                potential_dyads > thisDyad + exclusion_radius
            )
            potential_dyads = potential_dyads[still_valid]
            potential_scores = potential_scores[still_valid]
        
        if (i >= maxiter):
            print("/!\\ maxiter reached")
        
        dyads[chrom] = dyads[chrom][0:n_dyads]
        if (verbose):
            print("-- Found " + str(n_dyads) + " dyads...")
    
    print("~~~ Greedily Placed Nucleosomes ~~~")
    return dyads

def dynamic_P_G(S, N, tau = 1, beta = 0.5):
    '''
    This function is the forward process for the dynamic nucleosome placement
    algorithm that will produce a probability distribution for nucleosome
    placement. G(i) = sum(all scaled weights of configurations of subsequence
    0...i). Sourced from Field et al. 2008.
    
    Arguments:
        S: A list of ONE-INDEXED scores (scores[0] should be NaN, etc. It will
            not be read.) of CONSECUTIVE positions along a sequence of
            nucleosome START POSITIONS, (**NOT** DYAD POSITIONS). scores[i]
            should be the score of a nucleosome starting at position i in the
            genome (1-indexed).
        N: The number of positions in this stretch.
        tau: An apparent nucleosome density. Defaults to 1.
        beta: An apparent inverse temperature parameter. Defaults to 0.5.
    
    Note that G is 1-indexed. This function also returns a tuple containing (G,
    scale_factors).
    '''
    G = np.zeros(N + 1)
    scale_factors = np.ones(N + 1)
    
    G[0] = 1
    i = 1
    while (i <= N):
        if (i <= 146):
            G[i] = G[i-1]
        else:
            G[i] = G[i-1] + np.prod(
                scale_factors[(i-146):i]
            )*G[i-147]*tau*np.exp(beta*S[i-146])
        
        if (G[i] > 1e9):
            scale_factors[i] = 0.9
            G[i] = G[i]*0.9
        i += 1
    
    return G, scale_factors

def dynamic_P_K(S, scale_factors, N, tau = 1, beta = 0.5):
    '''
    This function is the scaled reverse process for the dynamic nucleosome
    placement algorithm that will produce a probability distribution for
    nucleosome placement. K(i) = sum(all scaled weights of configurations of
    subsequence i...N). Sourced from Field et al. 2008.
    
    Arguments:
        S: A list of ONE-INDEXED scores (scores[0] should be NaN, etc. It will
            not be read.) of CONSECUTIVE positions along a sequence of
            nucleosome START POSITIONS, (**NOT** DYAD POSITIONS). scores[i]
            should be the score of a nucleosome starting at position i in the
            genome (1-indexed). The scores should be subtracted by
            ln(scale_factor).
        scale_factors: An array of factors used to scale the results at each
            step.
        N: The number of positions in this stretch.
        tau: An apparent nucleosome density. Defaults to 1.
        beta: An apparent inverse temperature parameter. Defaults to 0.5.
    
    Note that K is 1-indexed.
    '''
    K = np.zeros(N + 2)
    
    K[N+1] = 1
    i = N
    while (i > 0):
        if (i >= N - 145):
            K[i] = scale_factors[i]*K[i+1]
        else:
            K[i] = scale_factors[i]*K[i+1] + np.prod(
                scale_factors[i:(i+146+1)]
            )*K[i+147]*tau*np.exp(beta*S[i])
        
        i -= 1
    
    return K

def dynamic_Q(S, N, tau = 1, beta = 0.5):
    '''
    This function will take in a set of scores for one consecutive block of the
    genome and return the scaled probability distribution for a nucleosome
    being present at each position on that block.
    
    Arguments:
        S: A list of ONE-INDEXED scores (scores[0] should be NaN, etc. It will
            not be read.) of CONSECUTIVE positions along a sequence of
            nucleosome START POSITIONS, (**NOT** DYAD POSITIONS). scores[i]
            should be the score of a nucleosome starting at position i in the
            genome (1-indexed).
        N: The number of positions in this stretch.
        tau: An apparent nucleosome density. Defaults to 1.
        beta: An apparent inverse temperature parameter. Defaults to 0.5.
        
        Note that Q is 1-indexed. This function also returns a tuple containing
        (Q, scale factors).
    '''
    G, scale_factors = dynamic_P_G(S, N, tau, beta)
    K = dynamic_P_K(S, scale_factors, N, tau, beta)
    
    Q = np.zeros(N + 1)
    
    for i in range(1, N - 146 + 1):
        Q[i] = G[i-1]*tau*np.exp(beta*S[i])*K[i+147] / K[1]
    
    return (Q, scale_factors)

def get_dynamic_probability_distribtion(scores, tau = 1, beta = 0.5):
    '''
    This function takes in a dictionary of scores and returns the probability
    distribution of a nucleosome being found at every location.
    
    Arguments:
        scores: The dictionary of dyad position scores across the whole genome.
        tau: An apparent nucleosome density. Defaults to 1.
        beta: An apparent inverse temperature parameter. Defaults to 0.5.
    '''
    print("~~~ Getting Dynamic Probability Distribution ~~~")
    
    # Find consecutive sequences
    print("- Finding consecutive stretches...")
    sequences = an.get_consecutive_stretches(scores)
    
    # Get distributions
    print("- Getting distributions...")
    probability_dist = {}
    for chrom in sequences:
        print("-- " + chrom + "...")
        probability_dist[chrom] = []
        for i in range(len(sequences[chrom])):
            S = np.append([np.nan], sequences[chrom][i][:,1])
            N = len(S) - 1 + 146
            
            Q, scale_factors = dynamic_Q(S, N, tau, beta)
            Q = Q[:len(S)]
            P = np.zeros(len(Q))
            for j in range(1, len(S)):
                P[j] = np.prod(scale_factors[j:(j+147)])*Q[j]
            P = P[1:]
            
            probability_dist[chrom].append(
                np.column_stack((sequences[chrom][i][:,0], P))
            )
    
    # Reformat return dictionary
    print("- Reformatting...")
    for chrom in probability_dist:
        probability_dist[chrom] = np.vstack(tuple(probability_dist[chrom]))
    
    print("~~~ Got Dynamic Probability Distribution ~~~")
    return probability_dist

def viterbi_placement(scores, tau = 1, beta = 0.5):
    '''
    This function uses the Viterbi algorithm (as detailed in Rabiner, 1989) to
    place the nucleosomes.
    
    Arguments:
        scores: The dictionary of dyad position scores across the whole genome.
        tau: An apparent nucleosome density. Defaults to 1.
        beta: An apparent inverse temperature parameter. Defaults to 0.5.
    '''
    print("~~~ Using Viterbi's Method to Place Nucleosomes ~~~")
    
    print("- Finding consecutive stretches...")
    sequences = an.get_consecutive_stretches(scores)
    
    print("- Doing Viterbi stuff...")
    locations = {}
    for chrom in sequences:
        print("-- " + chrom + "...")
        locations[chrom] = np.zeros(0, dtype = int)
        for I in range(len(sequences[chrom])):
            print(
                "--- Sequence "+str(I+1)+"/"+str(len(sequences[chrom]))+"..."
            )
            S = np.append([np.nan], sequences[chrom][I][:,1])
            N = len(S) - 1 + 146
            
            print("= Finding epsilon and phi...")
            epsilon = np.zeros(N + 1)
            phi = np.zeros(N + 1, dtype = int)
            for i in range(1,N+1):
                if (i <= 146):
                    epsilon[i] = epsilon[i-1]
                    phi[i] = 0
                else:
                    options = (
                        epsilon[i-1],
                        epsilon[i-147] + np.log(tau) + beta*S[i-146]
                    )
                    phi[i] = np.argmax(options)
                    epsilon[i] = max(options)
            
            print("= Backtracing...")
            i = N
            L = np.zeros(0, dtype = int)
            while (i > 0):
                if (phi[i] == 1):
                    L = np.append(L, i - 146)
                    i -= 147
                else:
                    i -= 1
            
            L -= 1
            
            locations[chrom] = np.concatenate(
                (locations[chrom], sequences[chrom][I][L,0].astype(int))
            )
    
    print("- Sorting...")
    locations = an.dict_sort(locations)
    
    print("~~~ Used Viterbi's Method to Place Nucleosomes ~~~")
    return locations

def viterbi_placement_with_overlap(scores, tau = 1, beta = 0.5, overlap = 40):
    '''
    This function uses the Viterbi algorithm (as detailed in Rabiner, 1989) to
    place the nucleosomes.
    
    Arguments:
        scores: The dictionary of dyad position scores across the whole genome.
        tau: An apparent nucleosome density. Defaults to 1.
        beta: An apparent inverse temperature parameter. Defaults to 0.5.
        overlap: The number of base pairs to allow nucleosomes to overlap.
    '''
    print("~~~ Using Viterbi's Method to Place Nucleosomes ~~~")
    
    print("- Finding consecutive stretches...")
    sequences = an.get_consecutive_stretches(scores)
    length = 147 - overlap
    print("- Doing Viterbi stuff...")
    locations = {}
    for chrom in sequences:
        print("-- " + chrom + "...")
        locations[chrom] = np.zeros(0, dtype = int)
        for I in range(len(sequences[chrom])):
            print(
                "--- Sequence "+str(I+1)+"/"+str(len(sequences[chrom]))+"..."
            )
            S = np.append([np.nan], sequences[chrom][I][:,1])
            N = len(S) - 1 + 146
            
            print("= Finding epsilon and phi...")
            epsilon = np.zeros(N + 1)
            phi = np.zeros(N + 1, dtype = int)
            for i in range(1,N+1):
                if (i <= length - 1):
                    epsilon[i] = epsilon[i-1]
                    phi[i] = 0
                else:
                    options = (
                        epsilon[i-1],
                        # We need to still reference the score 146 bp back!
                        epsilon[i-length] + np.log(tau) + beta*S[i-146]
                    )
                    phi[i] = np.argmax(options)
                    epsilon[i] = max(options)
            
            print("= Backtracing...")
            i = N
            L = np.zeros(0, dtype = int)
            while (i > 0):
                if (phi[i] == 1):
                    L = np.append(L, i - 146) # Nucleosome placed 146bp back!
                    i -= length
                else:
                    i -= 1
            
            L -= 1
            
            locations[chrom] = np.concatenate(
                (locations[chrom], sequences[chrom][I][L,0].astype(int))
            )
    
    print("- Sorting...")
    locations = an.dict_sort(locations)
    
    print("~~~ Used Viterbi's Method to Place Nucleosomes ~~~")
    return locations

def output_placements(placements, name):
    '''
    This function takes in a dictionary of nucleosome locations and outputs
    both a .tsv file and a .wig file for easy visualization.
    
    Arguments:
        placements: The dictionary of nucleosome dyad positions.
        name: The name of the directory/file. Should include the path to the
            file. Leave off the file extension.
    '''
    print("~~~ Outputting Nucleosome Placements ~~~")
    
    print("- Outputting .tsv file...")
    io.output_tsv(placements, name, fmt = "%d", chroms = True)
    
    print("- Outputting .wig file...")
    table = np.zeros((0, 4), dtype = "<U32")
    for chrom in placements:
        table = np.vstack(
            (
                table,
                np.column_stack(
                    (
                        np.repeat(chrom, np.shape(placements[chrom])[0]),
                        placements[chrom] - 73,
                        placements[chrom] + 73,
                        placements[chrom]
                    )
                )
            )
        )

    io.output_wig_interval_from_table(table, name + ".wig", highlight = True)
    
    print("~~~ Output Nucleosome Placements ~~~")

#%% Evaluation

def get_free_region_centers(reference_dyads, genome, blocks, padding = 73):
    '''
    This function takes in a set of reference 'true' dyad positions and finds
    regions of the genome that are devoid of any nucleosomes, returning a
    dictionary of 'free region centers,' where a 'free region' is a region of
    the genome at least 147bp wide that is devoid of nucleosome dyads with a
    'padding' area of width padding on either end that is also devoid of
    nucleosome dyads. These 'free regions' are then divided into 147bp windows,
    and the centers of these windows are returned.
    
    Arguments:
        reference_dyads: A dictionary of SORTED locations of dyads for some
            reference dataset.
        genome: The genome array dictionary.
        blocks: The array of blocks outlining which regions of the genome are
            excluded from analysis.
        padding: The width of the padding region on either end of a 'free
            region' that must be devoid of nucleosomes. Defaults to 73.
    '''
    bad_blocks = {}
    for chrom in blocks:
        if (np.sum(blocks[chrom][:,2]) > 0):
            bad_blocks[chrom] = blocks[chrom][
                blocks[chrom][:,2].astype(bool),
                0:2
            ]
    
    free_regions = {}
    for chrom in reference_dyads:
        potential_free_region_ends = np.append(
            reference_dyads[chrom],
            len(genome[chrom])
        )
        potential_free_region_starts = np.append([0], reference_dyads[chrom])
        dists = potential_free_region_ends - potential_free_region_starts
        valid_regions = dists - 2*padding > 147
        
        free_region_ends = potential_free_region_ends[valid_regions] - padding
        free_region_starts=potential_free_region_starts[valid_regions]+padding
        
        free_regions[chrom] = np.column_stack(
            (free_region_starts, free_region_ends)
        )
        
    overlaps = an.find_interval_overlaps(free_regions, bad_blocks)
    for chrom in overlaps:
        free_regions[chrom] = free_regions[chrom][overlaps[chrom] == 0,:]
        if (np.shape(free_regions[chrom])[0] == 0):
            free_regions.remove(chrom)
    
    free_centers = {}
    for chrom in free_regions:
        centers = []
        for i in range(np.shape(free_regions[chrom])[0]):
            start = free_regions[chrom][i,0]
            end = free_regions[chrom][i,1]
            thisCtr = start + 73
            centers.append(thisCtr)
            while (thisCtr + 147 < end - 73):
                thisCtr += 147
                centers.append(thisCtr)
        
        free_centers[chrom] = np.array(centers)
    
    return free_centers

def evaluate_genomic_nucleosome_placement(
    putative_dyads,
    reference_dyads,
    reference_free_centers,
    P,
    N,
    distance_threshold = 73
):
    '''
    This function takes in a set of putative dyad calls and a set of 'true'
    reference dyad calls and a set of 'true' centers of nucleosome free regions
    and returns the confusion matrix for the nucleosome placement. Iterating
    over every reference dyad: a "TP" occurs if a putative dyad was placed
    within distance_threshold of the true dyad; a "FN" occurs if there is no
    putative dyad placed within distance_threshold of the true dyad. Then,
    iterating over every free center: a "FP" occurs if there is a putative dyad
    placed within distance_threshold of the center; a "TN" occurs if there is
    no putative dyad placed within distance_threshold of the center.
    
    Arguments:
        putative_dyads: A dictionary of putative dyad locations.
        reference_dyads: A dictionary of reference 'true' dyad locations.
        reference_free_centers: A dictionary of reference 'true' non-nucleosome
            locations.
        P: The number of reference dyads there are in all (the number of all
            "positives").
        N: The number of free centers there are in all (the number of all
            "negatives").
        distance_threshold: The maximum distance a putative dyad and reference
            dyad/free region center can be to count as being a nucleosome call.
    '''
    TP = 0
    D = an.get_distance_matrix(putative_dyads, reference_dyads)
    for chrom in D:
        TP += np.sum(np.min(D[chrom], axis = 0) <= distance_threshold)
    FN = P - TP
    
    FP = 0
    D = an.get_distance_matrix(putative_dyads, reference_free_centers)
    for chrom in D:
        FP += np.sum(np.min(D[chrom], axis = 0) <= distance_threshold)
    TN = N - FP
    
    return {"TP": TP, "FN": FN, "FP": FP, "TN": TN}

def get_nucleosome_calls(
    reference_dyads,
    putative_dyads,
    putative_dyad_scores,
    distance_threshold = 73
):
    '''
    This function takes in a set of reference dyads and, for each one,
    generates a "call," which is the maximum score of all the putative dyads
    that are within distance_threshold. If no putative dyads are within
    distance_threshold of the reference dyad, a score of -np.inf is used.
    
    Arguments:
        reference_dyads: A dictionary of reference dyad positions.
        putative_dyads: A dictionary of putative dyad positions.
        putative_dyad_scores: A dictionary of scores corresponding to the dyad
            positions.
        distance_threshold: The distance threshold to consider a dyad call.
            Defaults to 73.
    '''
    D = an.get_distance_matrix(reference_dyads, putative_dyads)
    
    calls = {}
    for chrom in reference_dyads:
        calls[chrom] = np.zeros(len(reference_dyads[chrom])) - np.inf
    
    for chrom in D:
        for j in range(np.shape(D[chrom])[0]):
            theseScores = putative_dyad_scores[chrom][
                D[chrom][j,:] <= distance_threshold
            ]
            if (len(theseScores) > 0):
                calls[chrom][j] = np.max(theseScores)
    
    return calls

def get_confusion_matrix_from_calls(
    positive_call_list,
    negative_call_list,
    score_thresholds = None
):
    '''
    This function takes in two calls, one for positive positions and one for
    negative positions, and returns the confusion matrix at every score
    threshold.
    
    Arguments:
        positive_call_list: An array of positive calls, the scores. NOT a
            dictionary.
        negative_call_list: An array of negative calls, the scores. NOT a
            dictionary.
        score_thresholds: The array of score thresholds to test for. If None,
            will test every unique score value present in the calls. Defaults
            to None.
    
    Note that this function returns a tuple, (confusion matrix, thresholds).
    '''
    if (score_thresholds is None):
        thresh = np.concatenate((positive_call_list, negative_call_list))
        thresh = np.sort(np.unique(thresh))
        thresh = thresh[thresh != -np.inf]
    else:
        thresh = score_thresholds
    
    T = len(thresh)
    C = {
        "TP": np.zeros(T, dtype = int),
        "TN": np.zeros(T, dtype = int),
        "FP": np.zeros(T, dtype = int),
        "FN": np.zeros(T, dtype = int)
    }
    
    P = len(positive_call_list)
    N = len(negative_call_list)
    
    for i in range(T):
        C["TP"][i] = np.sum(positive_call_list >= thresh[i])
        C["FN"][i] = P - C["TP"][i]
        C["TN"][i] = np.sum(negative_call_list < thresh[i])
        C["FP"][i] = N - C["TN"][i]
    
    return (C, thresh)

def evaluate_genomic_placement_as_score_varies(
    reference_dyads,
    putative_dyads,
    putative_dyad_scores,
    genome,
    blocks,
    distance_threshold = 73,
    score_thresholds = None
):
    '''
    This function takes in a set of reference dyads, putative dyads, scores,
    and a distance threshold and a set of score thresholds and computes the
    confusion matrix for every score threshold.
    
    Arguments:
        reference_dyads: The dictionary of reference dyad positions.
        putative_dyads: The dictionary of putative dyad positions.
        putative_dyad_scores: The dictionary of scores corresponding to the
            putative dyads.
        genome: The genome array dictionary.
        blocks: The dictionary of blocks outlining which regions of the genome
            to exclude from analysis.
        distance_threshold: The maximum center-center distance between a
            reference dyad and putative dyad to be considered "correct."
            Defaults to 73.
        score_thresholds: The list of score thresholds to compute the confusion
            matrices for. Defaults to None, which will use every possible score
            threshold.
    '''
    print("~~~ Evaluating Nucleosome Placement Across Multiple Scores ~~~")
    ref_dyads = an.dict_sort(reference_dyads) 
    
    print("- Getting free region centers...")
    free_centers = get_free_region_centers(ref_dyads, genome, blocks)
    
    print("- Getting calls...")
    positive_calls = an.dict_concat(
        get_nucleosome_calls(
            ref_dyads,
            putative_dyads,
            putative_dyad_scores,
            distance_threshold = distance_threshold
        )
    )
    negative_calls = an.dict_concat(
        get_nucleosome_calls(
            free_centers,
            putative_dyads,
            putative_dyad_scores,
            distance_threshold = distance_threshold
        )
    )
    
    print("- Getting confusion matrices...")
    C = get_confusion_matrix_from_calls(
        positive_calls,
        negative_calls,
        score_thresholds = score_thresholds
    )
    
    print("~~~ Evaluated Nucleosome Placement Across Multiple Scores ~~~")
    return C

def get_dyad_scores(all_scores, dyad_locations):
    '''
    This function takes in a dictionary of scores across the whole genome and
    a SORTED list of genomic positions for the dyad locations and returns a
    dictionary of scores, where the scores correspond to dyad_locations.
    
    Arguments:
        all_scores: The dictionary of genome scores. The first column must be
            the genome location, while the second column must be the score.
        dyad_locations: The dictionary of dyad locations. MUST BE SORTED.
    '''
    scores = {}
    for chrom in dyad_locations:
        theseDyads = dyad_locations[chrom]
        theseScores = all_scores[chrom][np.argsort(all_scores[chrom][:,0]),:]
        
        scores[chrom] = theseScores[np.isin(theseScores[:,0], theseDyads),1]
    
    return scores

def get_distance_freq_table(putative_dyads, reference_dyads, R = (-100,100)):
    '''
    This function takes in putative and reference dyad positions and returns
    a frequency table of all the distances.
    
    Arguments:
        putative_dyads: The putative dyad positions.
        reference_dyads: The reference_dyad_positions.
        R: The range of distances to consider. Defaults to (-100, 100).
    '''
    freqTable = np.zeros((R[1]- R[0] + 1, 2), dtype = int)
    freqTable[:,0] = np.arange(R[0], R[1] + 1)
    for chrom in reference_dyads:
        for i in range(len(reference_dyads[chrom])):
            theseDists = putative_dyads[chrom] - reference_dyads[chrom][i]
            theseDists=theseDists[(theseDists >= R[0]) & (theseDists <= R[1])]
            for d in theseDists:
                freqTable[d - R[0],1] += 1
    
    return freqTable

def get_distance_freq_table_from_distances(distances):
    '''
    This function takes in a dictionary of center-center distances and returns
    a frequency table of all the distances.
    
    Arguments:
        distances: A dictionary of center-center distances between nucleosome
            positions.
    '''
    Min = min([np.min(distances[chrom]) for chrom in distances])
    Max = max([np.max(distances[chrom]) for chrom in distances])
    this_table = np.column_stack(
        (np.arange(Min, Max + 1), np.zeros(Max + 1 - Min, dtype = int))
    )
    
    for chrom in distances:
        for D in np.unique(distances[chrom]):
            this_table[this_table[:,0]==D,1] += np.sum(distances[chrom]==D)
    
    return this_table

def get_periodic_proportion(
        distance_frequencies,
        minor_out_positions,
        distance_threshold = 65,
        distance_tolerance = 1
    ):
    '''
    This function takes in a table of distance frequencies and returns the
    proportion of distances that line up with the period.
    
    Arguments:
        distance_frequencies: The table of frequencies of nucleosome center-
            center distances.
        minor_out_positions: An array listing the minor out positions along
            the length of the nucleosome.
        distance_threshold: Only consider nucleosomes placed at a distance
            equal to this or closer away from the true dyad.
        distance_tolerance: The tolerance, in terms of numbers of base pairs,
            of width around each periodic offset to consider "lining up" with
            the period.
    
    Note that this function returns a tuple containing (successes, total,
    expected_value).
    '''
    thisTable = distance_frequencies[
        np.abs(distance_frequencies[:,0]) <= distance_threshold,
        :
    ]
    
    total = 0
    for spot in minor_out_positions:
        if (abs(spot) <= distance_threshold):
            total += np.sum(
                thisTable[
                    thisTable[:,0] == spot,
                    1
                ]
            )

    return (total, np.sum(thisTable[:,1]))
