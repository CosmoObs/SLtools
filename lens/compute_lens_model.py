import compute_nfw_lens_parameters as KSmod # import KSmod (Gabriel's modules); 
import logging

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


