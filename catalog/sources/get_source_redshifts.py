
##@package get_source_redshifts [formerly zsdist] 
# Returns a distribution of source redshifts and redshift bin widths
#
#@param halo_model [defines the halo, for analytic models it is a list containing the halo parameters]
#@param zl
#@param cosmological_background_parameters
#@return a tuple with two lists, the first list contains the source redshifts, the second the bin size
def get_source_redshifts (zl, cosmological_background_parameters=(), halo_model=()):
    return  [1.5*zl, 2*zl , 2.5*zl, 3*zl], [0.5*zl, 0.5*zl,0.5*zl,0.5*zl] 


