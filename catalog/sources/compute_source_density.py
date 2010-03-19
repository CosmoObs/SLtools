#---------------------------------------------------------------------------------------------------------------------
# Maria's compute_source_density: I changed the way she uses the X variable (she used to name the catalog with argv; I define X inside the function with the name of the UDF catalog)
##@package compute_source_density COMPUTES THE SOURCE DENSITY BASED ON THE UDF CATALOG
# Based on the Hubble Ultra Deep Field (HUDF) catalog, computes the source density for a given redshift and redshift bin. 
#
#@param zs
#@param delta_zs
#@param zl
#@param source_selection_criteria (to get the catalog path and enhance_nsource, that artificially increases the source density by a multiplicative factor)
#
#@return Return the source density (increased by enhance_nsource)
def compute_source_density(zS, delta_zS, source_selection_criteria): 

	UDF = source_selection_criteria['source_catalog']
	enhance_nsource = float(source_selection_criteria['enhance_nsource']);
	X = open('%s' % (UDF)).readlines()
	z0 = []
	zs = []
	for i in range(len(X)):
		z0.append(float(X[i].split()[7]))              # redshifts from the UDF catalog
		if z0[i]>= (zS-delta_zS/2.0) and z0[i]< (zS+delta_zS/2.0): 
			zs.append(z0[i])
	N_area = len(zs)/43092.0
	#print min(zs)                               # or 11.97.0 arcmin**2
	return N_area * enhance_nsource
