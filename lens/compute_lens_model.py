
##@package compute_lens_model [formerly kappalens]
# converts halo model to lens model. if halo is a surface density map, the lens model will be a kappa map. 
#
# here, the lens model is the list of projected NFW parameters kappas and xsl 
#@param zl 
#@param zs
#@param halo_model (m200, c200,el,thetal) 
#@param cosmological_background_params (Omega_m, Omega_l)
#@return Return a dictionary with the lens model parameters
def compute_lens_model (zl, zs, halo_model, cosmological_background_params):
	from arcgeneratormodules import log
        import compute_NFW_lens_parameters as KSmod # import KSmod (Gabriel's modules); 
        m200=halo_model['m200'] 
        c200=halo_model['c200']
        el=halo_model['el'] 
        thetal=halo_model['thetal'] 
        omega_m = cosmological_background_params['omega_m'] 
        omega_l = cosmological_background_params['omega_l'] 
        # kappas and xsl
        #------------------------------------------------------------------------------
        kappas = KSmod.kappa_s(zl,zs,m200,c200,omega_m,omega_l) # NFW kappas calculated by Gabriel's module
        xsl = 60* KSmod.x_s(zl,m200,c200,omega_m,omega_l) # NFW xsl calculated by Gabriel's module. Converted from arcmin to arcsec
        #------------------------------------------------------------------------------
        log('Obtained kappas = %f and xsl = %f from Gabriel\'s module\n' % (kappas, xsl) )
        # creates a dictionary 
        return {'kappas':kappas, 'xsl':xsl,'el':el,'thetal':thetal} 


