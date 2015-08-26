#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package description
"""
import getpass
import os
import ds9
import read_catalogs as rcs
import glob

def main_ds9_control():
    """
    Main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """

    red_data = rcs.read_catalog("cs82_redmapper")

    sdss_vec = []
    reg_vec = []
    cs82_vec = []

    for i in range(len(red_data[0]["MEM_MATCH_ID"])):
        #print red_data[0]["MEM_MATCH_ID"][i]

        obj_id = red_data[0]["MEM_MATCH_ID"][i]

        imi = "/data2/redmapper_stuff/sdss_cuts/"
        imi += str(obj_id)
        imi += "_sdss_i0.fits"

        img = "/data2/redmapper_stuff/sdss_cuts/"
        img += str(obj_id)
        img += "_sdss_g0.fits"

        imr = "/data2/redmapper_stuff/sdss_cuts/"
        imr += str(obj_id)
        imr += "_sdss_r0.fits"

        cs82_im = glob.glob("/data2/redmapper_stuff/cs82_cuts/" + \
                          str(obj_id) + "_*.fits")
        reg_file = glob.glob("/data2/redmapper_stuff/members_reg/" + \
                          str(obj_id) + ".reg")

        sdss_vec.append([imi, imr, img])
        cs82_vec.append(cs82_im)
        reg_vec.append(reg_file)
        #print reg_vec[i], cs82_vec[i]

    d1 = ds9.ds9()

    imcs82 = "/data2/redmapper_stuff/cs82_cuts/16_S82m13p_y0.fits"

    reg = "/data2/redmapper_stuff/members_reg/16.reg"

    var = ""

    count = -1

    user_reg = []

    while var != "x" and count < len(sdss_vec)-1:
        if var == "":
            user_reg = []
            count += 1
            print "Loading image %i of %i" % (count+1, len(sdss_vec))
            print "Image ID = %i" % red_data[0]["MEM_MATCH_ID"][count]
            d1.set("frame delete all")
            load_rgb(sdss_vec[count][0], sdss_vec[count][1], \
                     sdss_vec[count][2], d1)


            for j in range(len(cs82_vec[count])):
                d1.set("frame new")
                d1.set("zoom 1.0")
                d1.set("file " + cs82_vec[count][j])
                d1.set("zoom 1.0")
            d1.set("frame next")
#            ds9_match(d1)
        elif var == "reg":
            print reg_vec[count][0]
            d1.set("regions " + reg_vec[count][0])
        elif var == "regr":
            d1.set("regions delete all")
        elif var == "regs":
            reg = d1.get("regions -system wcs").split("\n")
            user_reg = get_reg_arc(reg)
            save_reg_arc(user_reg, red_data[0]["MEM_MATCH_ID"][count])
        elif var == "color":
            d1.set("cmap value 1 0.5")
        elif var == "m":
            d1.set("match frame wcs")
        elif var[:4] == "goto":
            print int(var[5:])
            count = int(var[5:]) - 2
        else:
            print "Command", var, "does not exists"
            print var[:3]

        var = raw_input("Command: ")
#    d1.set("exit")

def save_reg_arc(user_reg, im_id):
    #print getpass.getuser()
    count = 0
    fname = str(im_id) + "_arc_" + getpass.getuser() + "_" + str(count) + ".reg"
    #print fname
    #print glob.glob(fname)
    while len(glob.glob(fname)) > 0:
        count += 1
        fname = str(im_id) + "_arc_" + getpass.getuser() + "_" + str(count) + \
                ".reg"
        #print fname
    reg_out = open(fname, "w")
    reg_out.write("fk5\n")
    for i in user_reg:
        reg_out.write(i + "\n")
    reg_out.close()
def get_reg_arc(reg_in):
    reg_out = []
    for regs_i in range(len( reg_in )):
        if "color=blue" in reg_in[regs_i] or \
            "color=red" in reg_in[regs_i] or \
            "color=magenta" in reg_in[regs_i] or \
            "color=white" in reg_in[regs_i]:
            reg_out.append(reg_in[regs_i])
    return reg_out


def ds9_match(ds9_obj):
    ds9_obj.set("frame match wcs")

def load_rgb(image_i, image_r, image_g, ds9_obj):
    """
    Load a rgb image into ds9
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    ds9_obj.set("rgb")
    ds9_obj.set("rgb red")
    ds9_obj.set("file " + image_i)
    ds9_obj.set("scale zscale")

    ds9_obj.set("rgb green")
    ds9_obj.set("file " + image_r)
    ds9_obj.set("scale zscale")

    ds9_obj.set("rgb blue")
    ds9_obj.set("file " + image_g)
    ds9_obj.set("scale zscale")
    ds9_obj.set("rgb lock wcs yes")
    ds9_obj.set("rgb lock crop yes")
    ds9_obj.set("rgb lock slice yes")
    ds9_obj.set("rgb lock bin yes")
    ds9_obj.set("rgb lock scale yes")
    ds9_obj.set("rgb lock colorbar yes")
    ds9_obj.set("rgb lock smooth yes")


if __name__ == '__main__':

    print "Hello world!"
    main_ds9_control()
