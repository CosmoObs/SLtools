#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package description
"""

import read_catalogs as rc

import sdss2image as si

import coordinates as cor

import matplotlib.pylab as plt

import position_match as pm

import pyfits

import numpy as np

import time

import os

import sltools.coordinate.wcs_conversion as sltools_wcs_conv

def main_redmapper():
    """
    Main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    mem_data, mem_names = rc.read_catalog("cs82_redmapper_members")
    red_data, red_names = rc.read_catalog("cs82_redmapper")


    print mem_names
    print red_names

    print len(red_data["RA"])
    ids = np.unique(mem_data["MEM_MATCH_ID"])

    for i in range(1, 10):
        id_tmp = red_data["MEM_MATCH_ID"][i]
        msk = (mem_data["MEM_MATCH_ID"] == id_tmp)
        plt.plot(mem_data["RA"][msk], mem_data["DEC"][msk], ".")

        msk = (red_data["MEM_MATCH_ID"] == id_tmp)
        plt.plot(red_data["RA"][msk], red_data["DEC"][msk], 'o', alpha = 0.5)
        plt.show()


    """
    print len(ids)
    for i in ids:
        msk = (mem_data["MEM_MATCH_ID"] == i)
        plt.plot(mem_data["RA"][msk], mem_data["DEC"][msk], ".")
    plt.show()
    """

def make_redmapper_reg():
    """
    Function to make region files from redmapper cs82 cluster ctalog
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    mem_data, mem_names = rc.read_catalog("cs82_redmapper_members")
    red_data, red_names = rc.read_catalog("cs82_redmapper")


    print mem_names
    print red_names
    print len(red_data["RA"])
    ids = np.unique(mem_data["MEM_MATCH_ID"])

    str_reg = ""

    print names

    for i in range(len(ids)):
        print i, str(ids[i]) + ":"
        str_reg = "fk5\n"

        f_reg = open( str(ids[i]) + ".reg", "w" )

        msk_mem = (mem_data["MEM_MATCH_ID"] == ids[i])
        for j in range(len(mem_data["RA"][msk_mem])):
            str_reg += "circle( %f, %f, %.1f\")\n" % \
                    (mem_data["RA"][msk_mem][j], mem_data["DEC"][msk_mem][j], 5)

        #print str_reg
        f_reg.write(str_reg)
        f_reg.close()
        #print "\n"

def download_sdss_coadd():
    """
    Download coadimages of redmapper clusters
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    mem_data, mem_names = rc.read_catalog("cs82_redmapper_members")
    red_data, red_names = rc.read_catalog("cs82_redmapper")

    for i in range(len(red_data["RA"])):
        ra_red = red_data["RA"][i]
        dec_red = red_data["DEC"][i]
        print i, ra_red, dec_red
        time.sleep(2)
        names = si.get_sdss_tiles_within(ra_red, dec_red, "g", win = 1, \
                                         stripe82 = True)
        print names
        time.sleep(1)
        names = si.get_sdss_tiles_within(ra_red, dec_red, "r", win = 1, \
                                         stripe82 = True)
        print names
        time.sleep(1)
        names = si.get_sdss_tiles_within(ra_red, dec_red, "i", win = 1, \
                                         stripe82 = True)
        print names


def red_sdss_ds9(cat_position):
    """
    Function to open ds9 with colored image and reg files with cluster members
    for the cs82_redmapper clusters.
    It uses comand line, use the package ds9_control.py instead of this!
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    red_data, red_names = rc.read_catalog("cs82_redmapper")
    red_id = red_data["MEM_MATCH_ID"][cat_position]
    names = si.get_sdss_tiles_within(red_data["RA"][cat_position], \
                   red_data["DEC"][cat_position], "g", win = 1, stripe82 = True)
    print names
    sdss_dir = os.environ["SDSS_IMAGE_DIR"]
    imi = []
    img = []
    imr = []
    for i in range(len(names)):
        image_def = names[i][:len(sdss_dir)+12] + "i" + \
                    names[i][len(sdss_dir)+13:]
        imi.append(image_def)
        image_def = names[i][:len(sdss_dir)+12] + "r" + \
                    names[i][len(sdss_dir)+13:]
        imr.append(image_def)
        image_def = names[i][:len(sdss_dir)+12] + "g" + \
                    names[i][len(sdss_dir)+13:]
        img.append(image_def)

    str_cmd = "ds9 -rgb " + imi[0] + " -scale zscale -green " + imr[0] + \
              " -scale zscale -blue " + img[0] + " -scale zscale -region " + \
              str(red_id) + ".reg"


    print imi
    print imr
    print img
    print str_cmd
    os.system(str_cmd)


def red_cluster_size():
    """
    Computes the size of redmapper clusters.
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    mem_data, mem_names = rc.read_catalog("cs82_redmapper_members")
    red_data, red_names = rc.read_catalog("cs82_redmapper")

    ids = np.unique(mem_data["MEM_MATCH_ID"])

    dist = []

    for i in range(len(ids)):
    #    print i, str(ids[i]) + ":"
        msk_mem = (mem_data["MEM_MATCH_ID"] == ids[i])


        max_res = pm.max_diff(red_data["RA"][i], red_data["DEC"][i], \
                      mem_data["RA"][msk_mem], mem_data["DEC"][msk_mem])

        dist.append(max_res)

        if max_res[0] > 7:
            print i, str(ids[i]) + ":"
            print red_data["RA"][i], red_data["DEC"][i]
            print max_res
        """
        names = si.get_sdss_tiles_within(red_data["RA"][i], \
            red_data["DEC"][i], "g", win = 6, stripe82 = True)
        time.sleep(1)
        names = si.get_sdss_tiles_within(red_data["RA"][i], \
            red_data["DEC"][i], "r", win = 6, stripe82 = True)
        time.sleep(1)
        names = si.get_sdss_tiles_within(red_data["RA"][i], \
            red_data["DEC"][i], "i", win = 6, stripe82 = True)
        time.sleep(1)
        """

    print max(dist)

    """
    plt.plot(red_data["RA"][963], red_data["DEC"][963], ".")
    plt.plot(mem_data["RA"][msk_mem], mem_data["DEC"][msk_mem], 'o', alpha = 0.5
)
    plt.plot(mem_data["RA"][msk_mem][0], mem_data["DEC"][msk_mem][0], 'o', alpha = 0.5)
    plt.show()
    """

#            print mem_data["RA"][msk_mem][j], mem_data["DEC"][msk_mem][j]
#            dist.append()

def find_cs82_image():
    """
    Function to find in which cs82 image the redmapper clusters are and create
    the command to make a cut_out using imcp.py from sltools
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    print "find_cs82_im"
    red_data, red_names = rc.read_catalog("cs82_redmapper")

    limits, imgs = cor.cs82_limits()
    cs82_path = os.environ["CS82_IMAGE_DIR"]


    for i in range(len(limits)):
        #print imgs[i]
        """
        if limits[i][0][1] > 300 and limits[i][0][0] < 300:
            print imgs[i]
            print limits[i]
        """


        str_out = ""
        for j in range(len(red_data["RA"])):
            
            if red_data["RA"][j] < limits[i][0][1] and \
               red_data["RA"][j] > limits[i][0][0] and \
               red_data["DEC"][j] < limits[i][1][1] and \
               red_data["DEC"][j] > limits[i][1][0]:
                hdu = pyfits.open(imgs[i])
                hdr = hdu[0].header
                pix_xy = sltools_wcs_conv.radec2xy(hdr, red_data["RA"][j], \
                                                   red_data["DEC"][j])
            
                str_out = "imcp.py " + imgs[i] + " " 
                str_out += str(pix_xy[0]) + " " + str(pix_xy[1]) + " "
                str_out += "-s 600,600 -J -o "
                str_out +=  str(red_data["MEM_MATCH_ID"][j]) + "_"
                str_out += imgs[i][len(cs82_path) + 1 : len(cs82_path) + 10]
                if pix_xy[0] > 0 and  pix_xy[0] < 21000 and \
                   pix_xy[1] > 0 and  pix_xy[1] < 21000:
                   print str_out
#                print red_data["MEM_MATCH_ID"][j], red_data["RA"][j], \
#                      red_data["DEC"][j]

def find_sdss_image():
    """
    Function to make create the sdss image for the redmapper clusters
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    red_data, red_names = rc.read_catalog("cs82_redmapper")

    str_swarp = ""
    str_imcp = ""
#    i_tmp = 897

    for i in range(31, len(red_data["RA"])):
        time.sleep(1)
        print i
        names = si.get_sdss_tiles_within(red_data["RA"][i], \
                              red_data["DEC"][i], "i", win = 5, stripe82 = True)

        sdss_out_im =  str(red_data["MEM_MATCH_ID"][i]) + "_sdss_i.fits "

        str_swarp = "swarp -IMAGEOUT_NAME " + sdss_out_im
        for j in range(len(names)):
            str_swarp += names[j] + " "

        os.system(str_swarp)
        hdu = pyfits.open(sdss_out_im)
        hdr = hdu[0].header
        pix_xy = sltools_wcs_conv.radec2xy(hdr, red_data["DEC"][i], \
                                            red_data["RA"][i])

        str_imcp = "imcp.py " + sdss_out_im + " " 
        str_imcp += str(pix_xy[0]) + " " + str(pix_xy[1])
        str_imcp += " -s 600,600 -J -o " + sdss_out_im[:-6]
        os.system(str_imcp)

if __name__ == '__main__':

    print "Hello world!"
    #main_redmapper()
#    make_redmapper_reg()
#    for i in range(11, 20):
#        red_sdss_ds9(i)
#    red_cluster_size()
#    find_cs82_image()
    find_sdss_image()
