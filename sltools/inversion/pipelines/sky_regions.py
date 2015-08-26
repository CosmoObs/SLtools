#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================

"""
FIXME - package description
"""
from __future__ import division

import pyfits
import math
from progress_bar import ProgressBar


import matplotlib.pylab as plt
import matplotlib.patches as patches
import numpy as np
import functions as fc

from matplotlib.path import Path
def plot_rectangle(v1, v2, v3, v4):
    """
    Plot shaded retangles
    #FIXME TBD
    Input
     - v1 float :
     - v2 float : 
     - v3 float : 
     - v4 float :
    Output
     -
    #FIXME - finish documentation
    """
    return 0
def main_sky_regions():
    """
    main function
    Input
     -
    Output
     -
    #FIXME - finish documentation
    """
    verts_stripe82 = [
        (-50.0, -1.25), # left, bottom
        (-50.0,  1.25), # left, top
        ( 59.0,  1.25), # right, top
        ( 59.0, -1.25), # right, bottom
        (0., 0.), # ignored
        ]

    verts_CS82 = [
        (-42.5, -1.00), # left, bottom
        (-42.5,  1.00), # left, top
        ( 45.0,  1.00), # right, top
        ( 45.0, -1.00), # right, bottom
        (0., 0.), # ignored
        ]

    dic_areas = {}
    dic_areas["Stripe 82"] = verts_stripe82
    dic_areas["CS82"] = verts_CS82

    codes = [Path.MOVETO,
             Path.LINETO,
             Path.LINETO,
             Path.LINETO,
             Path.CLOSEPOLY,
             ]

    fig = plt.figure()
    ax = fig.add_subplot(111)

    color = ["b", "r"]
    count = 0
    for i in dic_areas.keys():
        path = Path(dic_areas[i], codes)
        patch = patches.PathPatch(path, lw=2, facecolor = color[count], \
                                  alpha = 0.1, label = i)
        ax.add_patch(patch)
        count += 1

    plt.xlabel("RA")
    plt.ylabel("DEC")
    plt.legend()
    #plt.show()

def plot_obj():
    #---read sogras observed
    sogras_path = "/home/gcaminha/utils/catalogs/sogras_observed_cris.txt"
    sogras_pos = np.loadtxt( sogras_path, usecols = (1,2), unpack = True )
    sogras_id = np.loadtxt( sogras_path, usecols = [0],  dtype = "|S6")
    for i in range(len(sogras_pos[0])):
        if sogras_pos[0][i] > 180:
            sogras_pos[0][i] = sogras_pos[0][i] - 360
    #---end read sogras observed

    #---read readmapper CS82
    cs82_redmapper = pyfits.open( \
              "/home/gcaminha/utils/catalogs/cs82_redmapper_v3.21_catalog.fit" )
    redmapper_data = cs82_redmapper[1].data
    red_names = cs82_redmapper[1].columns.names
    #print "N-clusters redmapper_data =", len(redmapper_data["MEM_MATCH_ID"])
    #print red_names
    for i in range(len(redmapper_data["RA"])):
        if redmapper_data["RA"][i] > 180:
            redmapper_data["RA"][i] = redmapper_data["RA"][i] - 360
    #plot_hist_fits(red_names, redmapper_data, title = "cs82_redmapper")
    #---end read readmapper CS82

    #---read maxbcg_public_plus_lambda.fit
    max_plus = pyfits.open( \
                 "/home/gcaminha/utils/catalogs/maxbcg_public_plus_lambda.fit" )
    max_data = max_plus[1].data
    max_names = max_plus[1].columns.names
    #print "N-clusters max_data =", len(max_data["RA"])
    #print max_names18
    for i in range(len(max_data["RA"])):
        if max_data["RA"][i] > 180:
            max_data["RA"][i] = max_data["RA"][i] - 360
    #plot_hist_fits(max_names, max_data, title = "maxbcg")
    #---end read maxbcg_public_plus_lambda.fit


    #---read gmb_stripe82_clusters_toplevel_v1
    gmb = pyfits.open( \
        "/home/gcaminha/utils/catalogs/gmb_stripe82_clusters_toplevel_v1.fits" )
    gmb_data = gmb[1].data
    gmb_names = gmb[1].columns.names
    #print "N-clusters gmb_data =", len(gmb_data["RA"])
    #print gmb_names
    for i in range(len(gmb_data["ra"])):
        if gmb_data["ra"][i] > 180:
            gmb_data["ra"][i] = gmb_data["ra"][i] - 360
    #plot_hist_fits(gmb_names, gmb_data, title = "gmb_stripe82")
    #---end gmb_stripe82_clusters_toplevel_v1


    #---read stripe82_specObj-dr8.fits
    spc = pyfits.open( \
        "/home/gcaminha/utils/catalogs/stripe82_specObj-dr8.fits" )
    spc_data = spc[1].data
    spc_names = spc[1].columns.names
    print "N-cobjects spc_data =", len(spc_data["PLUG_RA"])
    #print spc_names
    #for i in range(len(spc_data["PLUG_RA"])):
    #    if gmb_data["PLUG_RA"][i] > 180:
    #        gmb_data["PLUG_RA"][i] = gmb_data["ra"][i] - 360
    #plot_hist_fits(gmb_names, gmb_data, title = "gmb_stripe82")
    #---end stripe82_specObj-dr8.fits

    # plot regions and objects #################################################
    
    main_sky_regions()
    plt.plot(sogras_pos[0], sogras_pos[1], "o", label = "sogras")
    #plt.plot(redmapper_data["RA"], redmapper_data["DEC"], ".", \
    #         label = "redmapper")

    #plt.plot(gmb_data["ra"], gmb_data["DEC"], ".", label = "gmb", alpha = .5)

    #plt.plot(max_data["RA"], max_data["DEC"], ".", label = "maxbcg")

    plt.plot(spc_data["PLUG_RA"], spc_data["PLUG_DEC"], ".", label = "spec")

    #plt.xlim(-52, 61)
    #plt.ylim(-1.3, 1.75)
    plt.legend(ncol = 2)
    plt.show()
    

    # match sogras #############################################################
    
    sogras_match = {}
    for i in sogras_id:
        sogras_match[i] = {}
    #print sogras_match
    #tolerance for the match in arcmin
    toll = 1
    separator = "|"

    str_out = "!SOGRAS_ID|!GMB|!MAX|!RED|"
    #print str_out
    for i in range(len(sogras_pos[0])):
        str_out = separator
        print i, "sogras id = ", sogras_match.keys()[i]
        out = min_diff_mask(sogras_pos[0][i], sogras_pos[1][i], \
                       gmb_data["ra"], gmb_data["dec"], tolerance = toll)

        str_out += sogras_match.keys()[i] + separator

        # Saving GMB data info #################################################
        if out[0] < toll:
            #print "GMB_DATA - %.2f" % (out[0]), "gmb-OK!"
            sogras_match[sogras_match.keys()[i]]["GMB"] = out[0]
            sogras_match[sogras_match.keys()[i]]["GMB_photz"] = \
                                                    gmb_data["redshift"][out[1]]
            str_out += "%.2f" % (out[0]) + separator
        else:
            #print "GMB_DATA - %.2f" % (out[0])
            sogras_match[sogras_match.keys()[i]]["GMB"] = -999
            sogras_match[sogras_match.keys()[i]]["GMB_photz"] = -999
            str_out += "---" + separator
        ########################################################################


        out = min_diff_mask(sogras_pos[0][i], sogras_pos[1][i], \
                       max_data["RA"], max_data["DEC"], tolerance = toll)

        # Saving MAX data info #################################################
        if out[0] < toll:
            #print "MAX_DATA - %.2f" % (out[0]), "max-OK!"
            sogras_match[sogras_match.keys()[i]]["MAX"] = out[0]
            sogras_match[sogras_match.keys()[i]]["MAX_photz"] = \
                                                           max_data["Z"][out[1]]
            str_out += "%.2f" % (out[0]) + separator
        else:
            #print "MAX_DATA - %.2f" % (out[0])
            sogras_match[sogras_match.keys()[i]]["MAX"] = -999
            sogras_match[sogras_match.keys()[i]]["MAX_photz"] = -999
            str_out += "---" + separator
        ########################################################################


        out = min_diff_mask(sogras_pos[0][i], sogras_pos[1][i], \
                       redmapper_data["RA"], redmapper_data["DEC"], \
                       tolerance = toll)

        # Saving RED data info #################################################
        if out[0] < toll:
            #print "RED_DATA - %.2f" % (out[0]), "red-OK!"
            sogras_match[sogras_match.keys()[i]]["RED"] = out[0]
            sogras_match[sogras_match.keys()[i]]["RED_photz"] = \
                                                     redmapper_data["Z"][out[1]]
            str_out += "%.2f" % (out[0]) + separator
        else:
            #print "RED_DATA - %.2f" % (out[0])
            sogras_match[sogras_match.keys()[i]]["RED"] = -999
            sogras_match[sogras_match.keys()[i]]["RED_photz"] = -999
            str_out += "---" + separator
        ########################################################################


        out = min_diff_mask(sogras_pos[0][i], sogras_pos[1][i], \
                       spc_data["PLUG_RA"], spc_data["PLUG_DEC"], \
                       tolerance = toll)

        # Saving SPC data info #################################################
        if out[0] < toll:
            sogras_match[sogras_match.keys()[i]]["SPC"] = out[0]
            sogras_match[sogras_match.keys()[i]]["SPEC_Z"] = \
                                                     spc_data["Z"][out[1]]
            str_out += "%.2f" % (out[0]) + separator
        else:
            #print "RED_DATA - %.2f" % (out[0])
            sogras_match[sogras_match.keys()[i]]["SPC"] = -999
            sogras_match[sogras_match.keys()[i]]["SPEC_Z"] = -999
            str_out += "---" + separator
        ########################################################################

        #print str_out

    dic_keys = sogras_match.keys()
    dic_keys.sort()
    for i in dic_keys:
        str_out = separator + i + separator
        #print i, sogras_match[i].keys()

        if sogras_match[i]["GMB"] == -999:
            str_out += "---" + separator
            str_out += "---" + separator
        else:
            str_out += "%.2f" % (sogras_match[i]["GMB"]) + separator
            str_out += "%.2f" % (sogras_match[i]["GMB_photz"]) + separator
        str_out += " " + separator

        if sogras_match[i]["MAX"] == -999:
            str_out += "---" + separator
            str_out += "---" + separator
        else:
            str_out += "%.2f" % (sogras_match[i]["MAX"]) + separator
            str_out += "%.2f" % (sogras_match[i]["MAX_photz"]) + separator
        str_out += " " + separator

        if sogras_match[i]["RED"] == -999:
            str_out += "---" + separator
            str_out += "---" + separator
        else:
            str_out += "%.2f" % (sogras_match[i]["RED"]) + separator
            str_out += "%.2f" % (sogras_match[i]["RED_photz"]) + separator
        str_out += " " + separator

        if sogras_match[i]["SPC"] == -999:
            str_out += "---" + separator
            str_out += "---" + separator
        else:
            str_out += "%.2f" % (sogras_match[i]["SPC"]) + separator
            str_out += "%.2f" % (sogras_match[i]["SPEC_Z"]) + separator

        print str_out
    

    # gmb-max match ############################################################
    """
    p = ProgressBar(len(gmb_data["RA"]))
    
    max_lambda = []
    max_n200 = []

    rich_gmb = []

    z_max = []
    z_gmb = []
    
    for i in range(len(gmb_data["RA"])):
        out = min_diff( gmb_data["ra"][i], gmb_data["dec"][i], \
                        max_data["ra"], max_data["dec"])
        if out[0] < 0.3:
            #print out
            #print redmapper_data["RA"][i], redmapper_data["DEC"][i]
            #print gmb_data["ra"][out[1]], gmb_data["dec"][out[1]]
            rich_gmb.append(gmb_data["ngal"][i])
            z_gmb.append(gmb_data["redshift"][i])

            max_lambda.append(max_data["LAMBDA"][out[1]])
            max_n200.append(max_data["N200"][out[1]])
            z_max.append(max_data["Z"][out[1]])

            p.update_time(i)
            print p

    print len(rich_gmb), " gmb-max maches"
    plt.plot([min(rich_gmb),max(rich_gmb)], [min(rich_gmb), max(rich_gmb)])
    plt.plot(rich_gmb, max_lambda, '.')
    plt.xlabel("ngal - gmb")
    plt.ylabel("LAMBDA - max")
    plt.show()

    plt.plot([min(rich_gmb),max(rich_gmb)], [min(rich_gmb), max(rich_gmb)])
    plt.plot(rich_gmb, max_n200, '.')
    plt.xlabel("ngal - gmb")
    plt.ylabel("N200 - max")
    plt.show()

    plt.plot([min(z_gmb),max(z_gmb)], [min(z_gmb), max(z_gmb)])
    plt.plot(z_gmb, z_max, '.')
    plt.xlabel("Z_PHOT - gmb")
    plt.ylabel("Z_PHOT - max")
    plt.show()
    """



    # red-gmb match ############################################################
    """
    p = ProgressBar(len(redmapper_data["RA"]))
    
    rich_red = []
    rich_gmb = []
    z_red = []
    z_gmb = []
    
    for i in range(len(redmapper_data["RA"])):
        out = min_diff( redmapper_data["RA"][i], redmapper_data["DEC"][i], \
                       gmb_data["ra"], gmb_data["dec"])
        if out[0] < 0.3:
            #print out
            #print redmapper_data["RA"][i], redmapper_data["DEC"][i]
            #print gmb_data["ra"][out[1]], gmb_data["dec"][out[1]]
            rich_red.append(redmapper_data["LAMBDA_CHISQ"][i])
            rich_gmb.append(gmb_data["ngal"][out[1]])
            z_red.append(redmapper_data["Z"][i])
            z_gmb.append(gmb_data["redshift"][out[1]])

            p.update_time(i)
            print p

    print len(rich_red), " red-gmb maches"
    plt.plot([min(rich_red),max(rich_red)], [min(rich_red), max(rich_red)])
    plt.plot(rich_red, rich_gmb, '.')
    plt.xlabel("LAMBDA_CHISQ - red")
    plt.ylabel("ngal - gmb")
    plt.show()


    plt.plot([min(z_red),max(z_red)], [min(z_red), max(z_red)])
    plt.plot(z_red, z_gmb, '.')
    plt.xlabel("Z_PHOT - red")
    plt.ylabel("Z_PHOT - gmb")
    plt.show()
    """

    #   red-max match ##########################################################
    """
    p = ProgressBar(len(redmapper_data["RA"]))

    max_lambda = []
    red_lambda = []

    max_n200 = []

    z_red = []
    z_max = []

    zspec_red = []
    zspec_max = []

    for i in range(len(redmapper_data["RA"])):
        out = min_diff( redmapper_data["RA"][i], redmapper_data["DEC"][i], \
                       max_data["ra"], max_data["dec"])
        if out[0] < 0.3:
            #print out
            red_lambda.append(redmapper_data["LAMBDA_CHISQ"][i])
            max_lambda.append(max_data["LAMBDA"][out[1]])

            max_n200.append(max_data["N200"][out[1]])

            z_red.append(redmapper_data["Z"][i])
            z_max.append(max_data["Z"][out[1]])

            if (redmapper_data["BCG_SPEC_Z"][i] != -1) and \
               ( max_data["Z_SPEC_BCG"][out[1]] != -1):
                zspec_red.append(redmapper_data["BCG_SPEC_Z"][i])
                zspec_max.append(max_data["Z_SPEC_BCG"][out[1]])

            p.update_time(i)
            print p

    print len(red_lambda), " red-max maches"

    plt.plot([min(red_lambda),max(red_lambda)], [min(red_lambda), \
             max(red_lambda)])
    plt.plot(red_lambda, max_lambda, '.')
    plt.xlabel("LAMBDA_CHISQ - red")
    plt.ylabel("LAMBDA - max")
    plt.show()

    plt.plot([min(red_lambda),max(red_lambda)], [min(red_lambda), \
             max(red_lambda)])
    plt.plot(red_lambda, max_n200, '.')
    plt.xlabel("LAMBDA_CHISQ - red")
    plt.ylabel("N200 - max")
    plt.show()

    plt.plot([min(z_red),max(z_red)], [min(z_red), max(z_red)])
    plt.plot(z_red, z_max, '.')
    plt.xlabel("Z_PHOT - red")
    plt.ylabel("Z_PHOT - max")
    plt.show()

    plt.plot([min(zspec_red),max(zspec_red)], [min(zspec_red), max(zspec_red)])
    plt.plot(zspec_red, zspec_max, '.')
    plt.xlabel("Z_SPEC - red")
    plt.ylabel("Z_SPEC - max")
    plt.show()
    """

def min_diff_mask(x_ref, y_ref, x_vec, y_vec, tolerance = 1.5):
    """
    Optimized functions to compute the mininum distance between a point and a array
    Input
     - x_ref  float : 
     - y_ref  float : 
     - x_vec  float : 
     - y_vec  float : 
     - toll   float : 
    Output
     -
    #FIXME - finish documentation
    """
    x_min = x_ref - tolerance
    x_max = x_ref + tolerance
    y_min = y_ref - tolerance
    y_max = y_ref + tolerance
    if (type(x_vec) != np.array):
        x_vec = np.array(x_vec)
    if (type(y_vec) != np.array):
        y_vec = np.array(y_vec)

    msk = (x_vec >= x_min) & (x_vec <= x_max) & \
          (y_vec >= y_min) & (y_vec <= y_max)

    x_vec_msk = x_vec[msk]
    y_vec_msk = y_vec[msk]
    #print len(x_vec_msk), len(y_vec_msk)

    if len(x_vec_msk) == 0:
        return 1000, -1

    diff_pos = []
    for i in range(len(x_vec_msk)):
        ra_diff = (x_ref - x_vec_msk[i]) * math.cos(y_ref - y_vec_msk[i])
        dist = fc.modulus( ra_diff, (y_ref - y_vec_msk[i]) )
        diff_pos.append( math.fabs(dist) )

    min_index = diff_pos.index( min(diff_pos) )

    min_index1 = np.where( x_vec == x_vec_msk[min_index] )
    #min_index2 = np.where( y_vec == y_vec_msk[min_index] )

    #print min_index1[0][0], min_index2[0][0]

    #print min(diff_pos)*60.0, min_index1[0][0]
    return min(diff_pos)*60.0, min_index1[0][0]

def min_diff(x_ref, y_ref, x_vec, y_vec):
    """
    functions to compute the mininum distance between a point and a array
    Input
     - x_ref  float : 
     - y_ref  float : 
     - x_vec  float : 
     - y_vec  float : 
    Output
     -
    #FIXME - finish documentation
    """
    diff_pos = []

    for i in range(len(x_vec)):
        ra_diff = (x_ref - x_vec[i]) * math.cos(y_ref - y_vec[i])
        dist = fc.modulus( ra_diff, (y_ref - y_vec[i]) )
        diff_pos.append( math.fabs(dist) )

#    print "----------"
#    print min(diff_pos) * 60.0, "arcmin"
    min_index = diff_pos.index( min(diff_pos) )
#    print x_ref, y_ref
#    print x_vec[min_index], y_vec[min_index]
#    print "----------"
    #print min(diff_pos)*60.0, min_index
    return min(diff_pos)*60.0, min_index

def plot_hist_fits(names, data, title = ""):
    for i in names:
        print data[i]
        print type(data[i][0])
        print i
        print len(data[i])
        #if (type(data[i][0]) != str) or (i != "R"):
        if (i != "redshift_code") and (i != "R") and (i != "Agc"):
            plt.hist(data[i])
            plt.xlabel(i)
            plt.title(title)
            plt.show()
if __name__ == '__main__':

    print "Hello world!"
    #main_sky_regions()
    plot_obj()
