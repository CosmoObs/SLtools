#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Determines the redshift bins in which the sources will be placed. """

##@package get_source_redshifts
# 
#
# Returns a distribution of source redshifts and redshift bin widths


#=======================================================================================================
def get_source_redshifts (zl, cosmological_background_parameters=(), halo_model=()):
	"""
	Returns a distribution of source redshifts and redshift bin widths

	Input:
	 - zl  <float> : lens redshift

	Output:
	 - zs_bins        <list> : List with the redshift bins (center) where the sources should be 
				   projected
	 - zs_bins_width  <list> : List with the width of each source redshift bin 

	"""
    return  [1.5*zl, 2*zl , 2.5*zl, 3*zl], [0.5*zl, 0.5*zl,0.5*zl,0.5*zl] 


