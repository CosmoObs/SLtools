#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
Package to select objects within the CS82 footprint
To use other catalogs you need to change lines 30-32 in orther to make the
program read the catalog you want
"""

import coordinates as coo
import read_catalogs as rcs
import numpy as np
import pyfits
import functions as fuc

import matplotlib.pylab as plt



def main_cs82_footprint_cut():
    """
    Main function
    Input
     -
    Output
     -
    """

    gmb_data, gmb_names, gmb_hdu = \
                                rcs.read_catalog("ecgmm_old", return_hdu = True)
    gmb_data.sort(order = "NGALS")
    fuc.s82_ra(gmb_data["RA"])
    print "gmb_names", gmb_names, "\n"

    msk = np.zeros(len(gmb_data["RA"]), dtype = bool)
    cs82_lim = coo.cs82_limits()[0]
    for i in range(len(cs82_lim)):
        fuc.s82_ra(cs82_lim[i][0])
        msk = msk | ( gmb_data["RA"] > cs82_lim[i][0][0] ) & \
                    ( gmb_data["RA"] < cs82_lim[i][0][1] ) & \
                    ( gmb_data["DEC"] > cs82_lim[i][1][0]) & \
                    ( gmb_data["DEC"] < cs82_lim[i][1][1])


    gmb_cs82_dic = {}
    for i in gmb_names:
        gmb_cs82_dic[i] = gmb_data[i][msk]

    gmb_cs82_columns = []
    for i in range(len(gmb_hdu[1].columns.formats)):
        gmb_cs82_columns.append( pyfits.Column(name = gmb_names[i], \
                                 format = gmb_hdu[1].columns.formats[i], \
                                 array = gmb_cs82_dic[gmb_names[i]]) )
    tbhdu_s82 = pyfits.new_table(gmb_cs82_columns)
    tbhdu_s82.writeto("test.fits")

    for i in range(len(cs82_lim)):
        plt.plot([cs82_lim[i][0][0], cs82_lim[i][0][1], cs82_lim[i][0][1], \
                  cs82_lim[i][0][0], cs82_lim[i][0][0]], \
                 [cs82_lim[i][1][0], cs82_lim[i][1][0], cs82_lim[i][1][1], \
                  cs82_lim[i][1][1], cs82_lim[i][1][0]], \
                 'b-')
    plt.plot(gmb_data["RA"][msk], gmb_data["DEC"][msk], ".")
    plt.show()

    print len(gmb_data["RA"][msk])

if __name__ == '__main__':

    print "Hello world!"
    main_cs82_footprint_cut()
