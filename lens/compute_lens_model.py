#!/usr/bin/env python
# ==================================
# Authors:
# Pedro Ferreira - pferreira@dfte.ufrn.br
# ==================================

""" Module computes NFW parameters \f$ \kappa_s \f$ e  \f$ r_s \f$"""

##@package compute_lens_model
#
# Computes NFW parameters \f$ \kappa_s \f$ e  \f$ r_s \f$


import compute_nfw_lens_parameters as KSmod # import KSmod (Gabriel's modules); 
import logging


#=======================================================================================================
def compute_lens_model (zl, zs, halo_model, cosmological_background_params):
	"""
	Computes NFW parameters \f$ \kappa_s \f$ e  \f$ r_s \f$.

	See documentation at http://twiki.linea.gov.br/bin/view/StrongLensing/KappaLens for details of
	the calculations.

	Input:
	 - zl                             <float> : lens redshift
	 - zs                             <float> : source redshift
	 - halo_model                      <dict> : (cluster) lens parameters: 'm200', 'c200', 'el' and
						    'thetal' meaning the lens mass (\f$ 10^{14}M_{\odot} 
						    \f$), concentration parameter, elipticity and 
						    position angle, respectively. 
	 - cosmological_background_params  <dict> : Cosmological parameters 'omega_m' (\f$ \Omega_m \f$)
						    and 'omega_l' (\f$ \Omega_{\Lambda} \f$) 

	Output: 
	 - <dict> : NFW parameters. The keys are 'kappas', 'xsl' (in arcsec),'el' and 'thetal'

	"""

        m200 = float(halo_model['m200']);
        c200 = float(halo_model['c200']);
        el = float(halo_model['el']);
        thetal = float(halo_model['thetal']);
        omega_m = float(cosmological_background_params['omega_m']);
        omega_l = float(cosmological_background_params['omega_l']);
        # kappas and xsl
        #------------------------------------------------------------------------------
        kappas = KSmod.kappa_s(zl,zs,m200,c200,omega_m,omega_l) # NFW kappas calculated by Gabriel's module
        xsl = 60* KSmod.x_s(zl,m200,c200,omega_m,omega_l) # NFW xsl calculated by Gabriel's module. Converted from arcmin to arcsec
        #------------------------------------------------------------------------------
        logging.debug('Obtained kappas = %f and xsl = %f from Gabriel\'s module' % (kappas, xsl) )
        # creates a dictionary 
        return {'kappas':kappas, 'xsl':xsl,'el':el,'thetal':thetal} 


