#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
Package to read predefined catalogs
"""

import time

import astropy.io.fits as pyfits

import socket

def path_stuff(key):
    """
    Function to return the path for several stuff. Useful when using different
    computers
    Input
     - key str : string to identify the "stuff"
    Output
     - path
    """

    pc_gcaminha_dic = { \
                        "cs82_redmapper_members" : \
                        "/home/gcaminha/utils/catalogs/" + \
                        "cs82_redmapper_v3.21_catalog_members.fit", \
                        #-------------------------------------------------------
                        "cs82_redmapper" : \
                        "/home/gcaminha/utils/catalogs/" + \
                        "cs82_redmapper_v3.21_catalog.fit", \
                        #-------------------------------------------------------
                        "sogras_cris" : \
                       "/home/gcaminha/utils/catalogs/sogras_observed_cris.txt",
                        #-------------------------------------------------------
                        "bcg_r200_old" : \
                        "/home/gcaminha/utils/catalogs/bcg_db_r200_full.fit",
                        #-------------------------------------------------------
                        "ecgmm_old" : \
                        "/home/gcaminha/utils/catalogs/Coadd_ecgmm_belended_2" \
                        "com.fit",
                        #-------------------------------------------------------
                        "gmbcg_dr7" : \
                        "/home/gcaminha/utils/catalogs/GMBCG_SDSS_DR7_PUB.fit",
                        #-------------------------------------------------------
                        "gmb_s82" : \
                        "/home/gcaminha/utils/catalogs/" \
                        "gmb_stripe82_clusters_toplevel_v1.fits"
                      }

    virgulino_dic = { \
                        "cs82_redmapper_members" : \
                        "/home/gbcaminha/utils/catalogs/" + \
                        "cs82_redmapper_v3.21_catalog_members.fit", \
                        #-------------------------------------------------------
                        "cs82_redmapper" : \
                        "/home/gbcaminha/utils/catalogs/" + \
                        "cs82_redmapper_v3.21_catalog.fit", \
                        #-------------------------------------------------------
                        "sogras_cris" : \
                      "/home/gbcaminha/utils/catalogs/sogras_observed_cris.txt",
                        #-------------------------------------------------------
                        "bcg_r200_old" : \
                        "/home/gbcaminha/utils/catalogs/bcg_db_r200_full.fit",
                        #-------------------------------------------------------
                        "ecgmm_old" : \
                        "/home/gbcaminha/utils/catalogs/Coadd_ecgmm_belended_"\
                        "2com.fit",
                        #-------------------------------------------------------
                        "gmbcg_dr7" : \
                        "/home/gbcaminha/utils/catalogs/GMBCG_SDSS_DR7_PUB.fit",
                        #-------------------------------------------------------
                        "gmb_s82" : \
                        "/home/gbcaminha/utils/catalogs/" \
                        "gmb_stripe82_clusters_toplevel_v1.fits"
                      }




    dic_all = {"pc-gcaminha" : pc_gcaminha_dic, "virgulino" : virgulino_dic}
    host_name = socket.gethostname()

    if (host_name in dic_all.keys()):
        #-----------------------------------------------------------------------
        if (key in dic_all[host_name].keys()):
            return dic_all[host_name][key]
        else:
            print "ERROR: path_stuff - \"" + key + "\" not found in \"" + \
                  host_name + "\""
        #-----------------------------------------------------------------------
    else:
        print "ERROR: path_stuff - hostname \"" + host_name + "\" is not known"
###############################################################################
def read_catalog(key, return_hdu = False):
    """
    Function to read catalogs
    Input
     - key str : string to identify the cataloge
    Output
     - data, names
    """

    if type(key) != str:
        print "ERROR: read_catalog - key must me a str"

    #---------------------------------------------------------------------------
    if key == "cs82_redmapper_members":
        redm = pyfits.open( path_stuff("cs82_redmapper_members") )
        redm_data = redm[1].data
        redm_names = redm[1].columns.names
        if return_hdu:
            return redm_data, redm_names, redm
        else:
           return redm_data, redm_names 
    #---------------------------------------------------------------------------
    elif key == "cs82_redmapper":
        red = pyfits.open( path_stuff("cs82_redmapper") )
        red_data = red[1].data
        red_names = red[1].columns.names
        if return_hdu:
            return red_data, red_names, red
        else:
            return red_data, red_names
    #---------------------------------------------------------------------------
    elif key == "sogras_cris":
        file_cris = open( path_stuff("sogras_cris"), "r" )
        id_cris = []
        ra_cris = []
        dec_cris = []
        for line in file_cris:
            id_cris.append(line.split()[0])
            ra_cris.append( float(line.split()[1]) )
            dec_cris.append( float(line.split()[2]) )
        dic_cris = {"ID" : id_cris, "RA" : ra_cris, "DEC" : dec_cris}
        return dic_cris, dic_cris.keys()
    #---------------------------------------------------------------------------
    elif key == "bcg_r200_old":
        fits_bcg = pyfits.open( path_stuff("bcg_r200_old") )
        bcg_data = fits_bcg[1].data
        bcg_names = fits_bcg[1].columns.names
        if return_hdu:
            return bcg_data, bcg_names, fits_bcg
        else:
            return bcg_data, bcg_names
    #---------------------------------------------------------------------------
    elif key == "ecgmm_old":
        fits_ecg = pyfits.open( path_stuff("ecgmm_old") )
        ecg_data = fits_ecg[1].data
        ecg_names = fits_ecg[1].columns.names
        if return_hdu:
            return ecg_data, ecg_names, fits_ecg
        else:
            return ecg_data, ecg_names
    elif key == "gmbcg_dr7":
        fits_gmbcg = pyfits.open( path_stuff("gmbcg_dr7") )
        gmbcg_data = fits_gmbcg[1].data
        gmbcg_names = fits_gmbcg[1].columns.names
        if return_hdu:
            return gmbcg_data, gmbcg_names, fits_gmbcg
        else:
            return gmbcg_data, gmbcg_names
    elif key == "gmb_s82":
        fits_gmbs82 = pyfits.open( path_stuff("gmb_s82") )
        gmbs82_data = fits_gmbs82[1].data
        gmbs82_names = fits_gmbs82[1].columns.names
        if return_hdu:
            return gmbs82_data, gmbs82_names, fits_gmbs82
        else:
            return gmbs82_data, gmbs82_names
    else:
        print "ERROR: read_catalog - key \"" + key + "\" is not known"
###############################################################################
def main_read_catalogs():
    """
    Main function
    Input
     -
    Output
     -
    """
    print "Hello world! from main_read_catalogs()"
###############################################################################
if __name__ == '__main__':

    print "Hello world!"
    main_read_catalogs()
