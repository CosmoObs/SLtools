

##@package compute_lens_model [formerly kappalens]
# converts halo model to lens model. if halo is a surface density map, the lens model will be a kappa map. 
#
# here, the lens model is the list of projected NFW parameters kappas and xsL 
#@param zl 
#@param zs
#@param halo_model (m200, c200,el,thetal) 
#@param cosmological_background_parameters (Omega_m, Omega_l)
#@return Return a dictionary with the lens model parameters
def compute_lens_model (zl, zs, halo_model, cosmological_background_parameters):
    import ksmod as KSmod # import KSmod (Gabriel's modules); 
    m200=halo_model['m200'] 
    c200=halo_model['c200']
    el=halo_model['el'] 
    thetal=halo_model['thetal'] 
    Omega_m = cosmological_background_parameters['Omega_m'] 
    Omega_l = cosmological_background_parameters['Omega_l'] 
    # kappas and xsL
    #------------------------------------------------------------------------------
    kappas = KSmod.kappa_s(zl,zs,m200,c200,Omega_m,Omega_l) # NFW kappas calculated by Gabriel's module
    xsL = 60* KSmod.x_s(zl,m200,c200,Omega_m,Omega_l) # NFW xsL calculated by Gabriel's module. Converted from arcmin to arcsec
    #------------------------------------------------------------------------------
    log('Obtained kappas = %f and xsL = %f from Gabriel\'s module\n' % (kappas, xsL) )
    # creates a dictionary 
    return({'kappas':kappas, 'xsL':xsL,'el':el,'thetal':thetal})   


