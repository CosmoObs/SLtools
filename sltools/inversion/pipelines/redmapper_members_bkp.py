#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package description
"""

import socket

import astropy.io.fits as pyfits

import matplotlib.pylab as plt

import read_catalogs as rc

import numpy as np

import functions as fc

def s82_ra( vec_ra):
    """
    docstring for s82_radec
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    #msk = (vec_ra > 180)
    for i in range(len(vec_ra)):
        if vec_ra[i] > 180:
            vec_ra[i] -= 360
#        print i, vec_ra[i]

    return vec_ra
def main_redmapper_members():
    """
    Main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    redm_data, redm_names = rc.read_catalog("cs82_redmapper_members")
    s82_ra(redm_data["RA"])

    sogr_data, sogr_names = rc.read_catalog("sogras_cris")
    s82_ra(sogr_data["RA"])

    red_data, red_names = rc.read_catalog("cs82_redmapper")

    ecgmm_data, ecgmm_names = rc.read_catalog("bcg_r200_old")
    s82_ra(ecgmm_data["RA"])

#    print redm_names
#    print red_names
    ids = np.unique(redm_data["MEM_MATCH_ID"])

    selec_reg = []
    for i in range(len(sogr_data["RA"])):
        selec_reg.append( fc.select_region(redm_data["RA"], redm_data["DEC"], \
                                     sogr_data["RA"][i], sogr_data["DEC"][i], \
                                     3) )
#        print i, len(selec_reg[i])

    """
    for i in range(len(selec_reg)):
        if len(selec_reg[i]) != 0:
            file_reg = open(sogr_data["ID"][i] + ".reg", "w")
            for j in range(len(selec_reg[i])):
                print i, j
            file_reg.close()
    """
    """
    msk_id =  (redm_data["MEM_MATCH_ID"] == ids[7])
    plt.scatter(redm_data["RA"][msk_id], redm_data["DEC"][msk_id], \
                s = np.exp(redm_data["MODEL_MAG"][msk_id] - \
                redm_data["MODEL_MAG"][msk_id].min()), alpha = 0.7)
    plt.plot(red_data["RA"][7], red_data["DEC"][7], "r", markersize = 10, \
             marker = "H", alpha = 0.5)
    plt.show()
    """
    
    for i_ids in np.unique(redm_data["MEM_MATCH_ID"]):
        msk = ( redm_data["MEM_MATCH_ID"] == i_ids )
        plt.plot(redm_data["RA"][msk], redm_data["DEC"][msk], ".")
    
    plt.plot(ecgmm_data["RA"], ecgmm_data["DEC"], 'mv', markersize = 7, \
             label = "ecgmm clusters")

    plt.plot(sogr_data["RA"], sogr_data["DEC"], "b^", \
            markersize = 7, label = "SOGRAS")

    for i in range(len(sogr_data["RA"])):
        circ = fc.make_circunference(2.0/60.0, sogr_data["RA"][i], \
                                     sogr_data["DEC"][i])
        if sogr_data["ID"][i][:1] == "S":
            plt.plot(circ[0], circ[1], 'b-')
        else:
             plt.plot(circ[0], circ[1], 'r-')
    plt.legend()
    plt.xlim(-51,61)
    plt.ylim(-1.25, 1.25)
    plt.show()
    
    print sogr_data["ID"]

if __name__ == '__main__':

    print "Hello world!"
    main_redmapper_members()
