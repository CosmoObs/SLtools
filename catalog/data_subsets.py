#!/usr/bin/env python
# ====================================================
# Authors:
# Bruno Moraes - bruno.a.l.moraes@gmail.com - 04/Apr/2011
# ====================================================


from __future__ import division


def create_data_subsets(data,col_name,cuts):

    """
    Creates subsets from data, cutting in numerical value ranges of a given field

    This function takes different sets of input value ranges for a given field in
    the data and uses them to create numpy masks that are then applied to the 
    original data in order to create subsets of it. These subsets and the 
    corresponding masks are returned.

    An arbitrary number of cuts can be performed, but they must be of one of the
    following three formats:

        - ['.',b] --> performs the cut data.field(col_name) < b
        - [a,b] ----> performs the cut a <= data.field(col_name) < b
        - [a,'.'] --> performs the cut a <= data.field(col_name)

    There's currently no support for cuts in string fields or to change from
    semi-open to open or closed intervals.


    Input:
     - data                        hdr : Input data
     - col_name                    str : Name of the column chosen to perform
                                         the cuts
     - cuts  [[a_1,b_1],...,[a_n,b_n]] : Cuts for different subsets

    Output:
     - cut_data         [pyfits_bla_1,...,pyfits_bla_n] : List of data subsets
                                            (each data subset: ndim=2, dtype=)
     - masks    [ndarray_1,...,ndarray_n] : List of masks used to create the 
                                            data subsets 
                                            (each mask: ndim=1,dtype=bool) 
    """

    # Try/Except if the field doesn't exist in the catalog?    

    cut_data = []
    masks = []

    for i in range(len(cuts)):

        if cuts[i][0] == '.':

            masks.append(data.field(col_name) < cuts[i][1])

        elif cuts[i][1] == '.':
    
            masks.append(cuts[i][0] <= data.field(col_name))

        else:

            masks.append((cuts[i][0] <= data.field(col_name)) & (data.field(col_name) < cuts[i][1]))


        cut_data.append(data[masks[i]])

    return cut_data, masks


def create_matched_subsets(data_A,data_B,col_name,cuts):

    """
    Creates subsets of data A and B, cutting in numerical value ranges of a 
    given field in data A

    This function creates subsets of matched data A and B according to cuts
    performed in a given field of data A, returning two lists, containing the 
    subsets of A and B. Input data A and B must have the same number of rows,
    ideally coming from a matching procedure.

    For details of the cutting procedure, check the documentation of the
    function create_data_subsets.


    Input:
     - data_A                      hdr : Input data set A
     - data_B                      hdr : Input data set B (same size as A)
     - col_name                    str : Name of the column chosen to perform
                                         the cuts
     - cuts  [[a_1,b_1],...,[a_n,b_n]] : Cuts for different subsets

    Output:
     - cut_data_A     [pyfits_bla_1,...,pyfits_bla_n] : List of data subsets
                                            (each data subset: ndim=2, dtype=)
     - cut_data_B     [pyfits_bla_1,...,pyfits_bla_n] : List of data subsets
                                            (each data subset: ndim=2, dtype=)
    """

    cut_data_A, masks_A = create_data_subsets(data_A,col_name,cuts)

    cut_data_B = []

    for i in range(len(cut_data_A)):
        
        cut_data_B.append(data_B[masks_A[i]])

    return cut_data_A, cut_data_B
