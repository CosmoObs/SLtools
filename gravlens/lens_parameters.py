##@package lens_parameters
# Creates 2 strings that are used as input for gravlens
#
#@param lens_model (kappas, xsl,el,thetal) currently a dictionary with the lens model parameters 
#@param gravlens_params (set of gravlens parameters: gridhi1, maxlev, etc)
#@return inputlens, setlens
def lens_parameters(lens_model, gravlens_params): # lens_model= [kappas, xsl, el, thetal]

	setlens = 'nfw %f 0 0 %f %f 0 0 %f 0 0 '% ( float(lens_model['kappas']), float(lens_model['el']), float(lens_model['thetal']), float(lens_model['xsl']) )  # (kappas, el, thetal, xsl)
	inputlens = """set gridhi1=%0.12f
set ngrid1 = %d
set ngrid2 = %d
set xtol = %0.12f 
set crittol = %0.12f 
set inttol = %0.12f 
set maxlev = %d
set gallev = %d 
set imglev = %d 
startup 1 1
%s
0 0 0 0 0 0 0 0 0 0\n""" % ( float(gravlens_params['gridhi1']),
                             int(gravlens_params['ngrid1']),
                             int(gravlens_params['ngrid2']),
                             float(gravlens_params['xtol']),
                             float(gravlens_params['crittol']),
                             float(gravlens_params['inttol']),
                             int(gravlens_params['maxlev']),
                             int(gravlens_params['gallev']),
                             int(gravlens_params['imglev']), setlens) 
# xtol: Tolerance on image positions in numerical root finding.
	return inputlens, setlens
