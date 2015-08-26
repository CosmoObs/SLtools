#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package description
"""
import numpy as np
import read_catalogs as rcat
import functions as func
import matplotlib.pylab as plt

def main_catalog_analysis():
    """
    Main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    print "Hello world! from main_catalog_analysis()"
    red_data, red_names = rcat.read_catalog("cs82_redmapper")
    func.s82_ra(red_data["RA"])

    gmbs82_data, gmbs82_names = rcat.read_catalog("gmb_s82")
    func.s82_ra(gmbs82_data["RA"])

    ecg_data, ecg_names = rcat.read_catalog("ecgmm_old")
    func.s82_ra(ecg_data["RA"])

    red_data, red_names = rcat.read_catalog("cs82_redmapper")
    func.s82_ra(red_data["RA"])

    sog_data, sog_names = rcat.read_catalog("sogras_cris")
    func.s82_ra(sog_data["RA"])

    print "gmb_stripe82_clusters_toplevel_v1.fits, %i objects:" % \
           len(gmbs82_data["RA"])
    print gmbs82_names


    print "\n Coadd_ecgmm_belended_2com.fit, %i objects:" % \
          len(ecg_data["RA"])
    print ecg_names



    print "\ncs82_redmapper_v3.21_catalog.fit, %i objects:" % \
           len(red_data["RA"])
    print red_names

    """
    plt.plot(red_data["RA"], red_data["DEC"], ".")
    plt.plot(gmbs82_data["RA"], gmbs82_data["DEC"], ".")
    plt.plot(sog_data["RA"], sog_data["DEC"], "o", alpha = 0.5)
    plt.show()
    plt.plot(ecg_data["RA"], ecg_data["DEC"], ".")
    plt.show()
    """
    msk = (ecg_data["GM_SCALED_NGALS"] > 20)
    print len( ecg_data["GM_SCALED_NGALS"][msk] )
    """
    plt.hist(np.log10(red_data["LAMBDA_CHISQ"]))
    plt.title("cs82_redmapper_v3.21_catalog.fit")
    plt.xlabel("log10(LAMBDA_CHISQ)")
    plt.ylabel("counts")
    plt.show()
    """
if __name__ == '__main__':

    print "Hello world!"
    main_catalog_analysis()
