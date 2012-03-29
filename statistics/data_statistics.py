from __future__ import division
import numpy as np

import pyper

def get_quantile_errors(data,quantiles):

    """
    Gets arbitrary quantile values for given data

    This function uses the PypeR library to call commands from the R language
    used to calculate an arbitrary number of quantile intervals of a numpy
    array. In the case of bi-dimensional numpy arrays, it will calculate
    these statistics along the first axis, meaning that each line corresponds
    to an independent data set.


    Input:
     - data        numpy.ndarray : data input, each line being a data set
                                   (ndim=2,dtype=float)
     - quantiles   numpy.ndarray : percentages (ndim=1,dtype=float)
                                   Ex: quartiles : array([0.25,0.5,0.75])

    Output:
     - qtl_errors  numpy.ndarray : quantile error values (ndim=2,dtype=float) 
     
    ---
    """

    myR = pyper.R()

    myR['data'] = np.transpose(data)
    
    myR['quantiles'] = quantiles

    # Calculate the quantiles for each data line
    
    myR("""qtls <- t(sapply(as.data.frame(data),function(x) quantile(x,
                     quantiles,names=FALSE)))""")

    qtl_errors = myR['qtls']

    del myR

    return qtl_errors

