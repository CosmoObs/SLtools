"""Module to calculate the magnitude limit of a survey"""


##@file maglimit
#
# ============================================================================
# Authors:
# Bruno Moraes - bruno.a.l.moraes@gmail.com - 01/Nov/2011 - Update:19/Feb/2014
# ============================================================================

from __future__ import division
import numpy as np
from scipy.optimize import curve_fit


def create_counts_data(mag_data, binsize):
    """
    Bin 1d magnitude data and calculate number counts and Poisson errors.

    Input:
     - mag_data        np.array : Magnitude data array
     - binsize           float : Bin size
    Output:
     - counts_data    np.array : [bin values, counts, Poisson errors]
    """

    # Define bin limits to be closest lower/higher int to mag_data min/max

    bininf = int(mag_data.min())
    binsup = int(mag_data.max()+1) # int rounds to lower integer
    bins = np.arange(bininf, binsup, binsize)

    counts, binedges = np.histogram(mag_data, bins=bins)

    bincenters = binedges[0:-1] + binsize/2
    counts_data = np.array([bincenters, counts, np.sqrt(counts)])

    return counts_data


def counts_model(mag, a, b):
    """
    Define the analytical form of galaxy/star number counts per pointing.

    The expected form, for theoretical reasons, is a*10^(b*magnitude),
    where magnitude is in the AB system and b is close to 0.3.

    Input:
     - mag        float : Magnitude value
     - a          float : Parameter
     - b          float : Parameter
    Output:
     - f          float : Count model value
    """

    return a*10**(b*mag)


def _fit_counts_data(counts_data, initparams=(1, 0.3), minf=None, msup=None):

    """
    Fit count data to analytical model in a selected magnitude range.

    If minf and msup are equal to None, the fitting range is automatically
    calculated using the peak of the counts distribution as a reference
    point.

    Input:
     - counts_data    np.array : Counts in the format [bins, counts, errors]
     - initparams        tuple : Initial parameter guesses
     - minf              float : Fitting range inferior limit
     - msup              float : Fitting range superior limit
    Output:
     - fitparams      np.array : Fitted parameters in format [afit, bfit]
    """

    if minf is None:
        minf = counts_data[0][counts_data[1].argmax()] - 3

    if msup is None:
        msup = counts_data[0][counts_data[1].argmax()] - 1 

    # Bug if the user doesn't pass a valid int or float
    mask = (counts_data[0] > float(minf)) & (counts_data[0] < float(msup))
    counts_data = counts_data[:, mask]

    # Perform fit

    mag = counts_data[0]
    counts = counts_data[1]
    counts_errors = counts_data[2]

    fitparams, __ = curve_fit(counts_model, mag, counts, initparams,
                               counts_errors)

    return fitparams


def _calc_rel_diff(counts_data, fitparams):

    """
    Calculate the relative difference between number count data and model

    If data_i > model_i, values are negative. Else, values are positive.

    Input:
     - counts_data  np.array : Counts in the format [bins, counts, errors]
     - fitparams    np.array : Fitted parameters in format [afit, bfit]
    Output:
     - rel_diff     np.array : Relative difference per magnitude bin
    """

    a = fitparams[0]
    b = fitparams[1]

    rel_diff = 1 - counts_data[1]/counts_model(counts_data[0], a, b)

    return rel_diff


def get_mag_limit_counts(mag_data, binsize, initparams=(1, 0.3), minf=None,
                         msup=None, comp_limit=0.8):
    """
    Get magnitude limit from the number counts in magnitude bins.

    For large field imaging surveys, there exists an empirical formula
    for the number of objects detected as a function of magnitude,
    given by N = a*10^{b*mag}, where a and b are constants and mag
    is the magnitude in the AB system. Deviation from this law at faint
    magnitudes indicates a loss of completeness in the observations. This
    function creates number count data by binning a numpy vector of all
    object magnitudes in a given field, fits the model to a selected
    magnitude range, compares the fit to the data at the faint end and
    obtains the magnitude value at which completeness drops below a given
    threshold. It returns this magnitude value.

    Input:
     - mag_data      np.array : Magnitude data array
     - binsize          float : Bin size
     - initparams       tuple : Initial parameter guesses
     - fitparams     np.array : Fitted parameters in format [afit, bfit]
     - minf             float : Fitting range inferior limit
     - msup             float : Fitting range superior limit
     - comp_limit       float : Completeness limit
    Output:
     - maglim           float : Magnitude limit
    """

    # Defend against bad input - Add code later

    counts_data = create_counts_data(mag_data, binsize)

    kwargsfit = {'initparams': initparams}
    if minf is not None:
        kwargsfit['minf'] = minf
    if msup is not None:
        kwargsfit['msup'] = msup

    fitparams = _fit_counts_data(counts_data, **kwargsfit)

    rel_diff = _calc_rel_diff(counts_data, fitparams)

    mask = ( rel_diff < (1 - comp_limit) )

    # Bug if all mask entries are false

    maglim = counts_data[0][mask].max()

    #if plot is not None:

    return maglim

