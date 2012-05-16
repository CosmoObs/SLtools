#!/usr/bin/env python
# ========================================================
# Authors:
# Bruno Moraes - bruno.a.l.moraes@gmail.com - 04/Apr/2012
# Aldee Charbonnier - aldeebaran@gmail.com - 16/Apr/2012
# ========================================================

from __future__ import division
import sys
import numpy as np
import pyfits
import json
import math
import matplotlib.pyplot as plt
import argparse
from argparse import RawDescriptionHelpFormatter

from sltools.catalog import data_subsets as dts
from sltools.catalog import fits_data as fd
from sltools.catalog import ascii_data as ascd
from sltools.io import io_addons as ioadd
from sltools.plot import plot_templates as pltemp
from sltools.statistics import data_statistics as dstat
from sltools.catalog.cs82 import morpho_comparisons as mcomp

plot_formats = ['+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+', '+']
plot_colors = ['y', 'c', 'm', 'b', 'g', 'r', 'k', 'y', 'c', 'm', 'b', 'g']


def main(cat_A, cat_B, output_folder, plots, matching, multiple, radius, 
         coord_A_names, coord_B_names):

    ioadd.create_folder(output_folder)

    # Open FITS catalogs
    tbhdu_A = pyfits.open(cat_A, ignore_missing_end=True, memmap=True)[2]
    tbhdu_B = pyfits.open(cat_B, ignore_missing_end=True, memmap=True)[2]

    # Perform standard cuts on both catalogs
    cut_tbhdu_A = fd.sample_entries(tbhdu_A, CLASS_STAR=(0, 0.95),
                  MAG_AUTO=(0, 24), MAGERR_AUTO=(0, 0.2172),
                  FLAGS=(0, 3))
    data_A = cut_tbhdu_A.data

    cut_tbhdu_B = fd.sample_entries(tbhdu_B, CLASS_STAR=(0, 0.95),
                  MAG_AUTO=(0, 24), MAGERR_AUTO=(0, 0.2172),
                  FLAGS=(0, 3))
    data_B = cut_tbhdu_B.data

    # Perform matching when necessary, *for the moment only returns nearest 
    # neighbor*.

    if matching:

        data_A, data_B = mcomp.htm_matching(data_A, data_B, coord_A_names, 
                         coord_B_names, radius=radius, 
                         exclude_multiple=multiple, n_neigh=1,depth=10)

    else: pass

    

    # Import plotting parameters to a dictionary
    # Define fieldnames (get them from file later)
    fieldnames = ['FIELDNAME_A', 'UNIT_A_1', 'UNIT_A_2', 'FIELDNAME_B',
                  'UNIT_B_1', 'UNIT_B_2', 'FIELDNAME_G', 'UNIT_G_1',
                  'UNIT_G_2', 'CUTS_G', 'TITLE_I', 'XMIN_I', 'XMAX_I',
                  'YMIN_I', 'YMAX_I', 'XLABEL_I', 'YLABEL_I', 'TITLE_II',
                  'XMIN_II', 'XMAX_II', 'YMIN_II', 'YMAX_II', 'XLABEL_II',
                  'YLABEL_II', 'TITLE_III', 'XMIN_III', 'XMAX_III',
                  'YMIN_III', 'YMAX_III', 'XLABEL_III', 'YLABEL_III']

    # Create dictionary and get number of plots
    plot_dict = ascd.dict_from_csv(plots, fieldnames, header_lines=1,
                                   delimiter='\t')
    n = len(plot_dict['FIELDNAME_A'])
    print 'Number of Plots (I+II+III):', n

    # Loop over the different subsets of plots
    for i in range(n):
        # Affect plot variables
        field_A = plot_dict['FIELDNAME_A'][i]
        units_A = [float(plot_dict['UNIT_A_1'][i]),
                   float(plot_dict['UNIT_A_2'][i])]

        field_B = plot_dict['FIELDNAME_B'][i]
        units_B = [float(plot_dict['UNIT_B_1'][i]),
                   float(plot_dict['UNIT_B_2'][i])]

        field_G = plot_dict['FIELDNAME_G'][i]
        units_G = [float(plot_dict['UNIT_G_1'][i]),
                   float(plot_dict['UNIT_G_2'][i])]

        print '--- Plot', i + 1, ' ---'
        print field_A, '(cat A),', field_B, '(cat B) vs', field_G

        cuts_G = np.array(json.loads(plot_dict['CUTS_G'][i]))
        m = len(cuts_G)

        # Cuts for G are currently given with the preferred units (e.g. arcsec
        # instead of degrees). To apply these cuts, they must be rescaled.
        cuts_G_rescaled = (cuts_G - units_G[1]) / units_G[0]

        # Divide the catalogs in subsets defined by the cuts list
        cut_data_A, cut_data_B = dts.create_matched_subsets(
            data_A, data_B, field_G, cuts_G_rescaled)

        # Create data lists for plots
        data1 = []
        data2 = []
        dataG = []
        colors = []
        formats = []

        for j in range(m):
            colors.append(plot_colors[m-1-j])
            formats.append(plot_formats[j])

            data_rescaled = (units_A[0] * cut_data_A[m-1-j].field(field_A) +
                             units_A[1])
            data1.append(data_rescaled)
            del(data_rescaled)

            data_rescaled = (units_B[0] * cut_data_B[m-1-j].field(field_B) +
                             units_B[1])
            data2.append(data_rescaled)
            del(data_rescaled)

            data_rescaled = (units_G[0] * cut_data_A[m-1-j].field(field_G) +
                             units_G[1])
            dataG.append(data_rescaled)
            del(data_rescaled)

        # ========================================================
        # Create type I figure
        # ========================================================

        title_I = plot_dict['TITLE_I'][i]
        xlims_I = [float(plot_dict['XMIN_I'][i]),
                   float(plot_dict['XMAX_I'][i])]
        ylims_I = [float(plot_dict['YMIN_I'][i]),
                   float(plot_dict['YMAX_I'][i])]
        xlabel_I = plot_dict['XLABEL_I'][i]
        ylabel_I = plot_dict['YLABEL_I'][i]

        # Diagonal line on the figure, the ideal case
        min_temp = xlims_I[0]
        max_temp = xlims_I[1]
        xlims_id = [min_temp, max_temp]
        min_temp = ylims_I[0]
        max_temp = ylims_I[1]
        ylims_id = [min_temp, max_temp]
        del(min_temp, max_temp)

        # Management of reversed axes (e.g. for magnitude)
        if xlims_id[0] > xlims_id[1]:
            min_temp = xlims_id[0]
            max_temp = xlims_id[1]
            xlims_id[0] = max_temp
            xlims_id[1] = min_temp
            del(min_temp, max_temp)

        if ylims_id[0] > ylims_id[1]:
            min_temp = ylims_id[0]
            max_temp = ylims_id[1]
            ylims_id[0] = max_temp
            ylims_id[1] = min_temp
            del(min_temp, max_temp)

        ideal_x_I = np.arange(xlims_id[0], xlims_id[1], 0.01)
        ideal_y_I = np.arange(ylims_id[0], ylims_id[1], 0.01)

        # Legend of type I figure
        legend_loc_I = 'upper left'
        legend_I = []
        for j in range(m):
            legend_fig = (str(cuts_G[m-1-j][0]) + ' < ' + field_G +
                          ' < ' + str(cuts_G[m-1-j][1]))
            legend_I.append(legend_fig)
            del(legend_fig)

        # Plot and save type I figure
        pltemp.plot_data(data1, data2, formats, colors,
                         legend_I, legend_loc_I, xlabel_I, ylabel_I,
                         xlims_I, ylims_I, title_I)
        plt.plot(ideal_x_I, ideal_y_I, '--', color='k', label='')
        plt.xlim(xlims_I[0], xlims_I[1])
        plt.ylim(ylims_I[0], ylims_I[1])

        plt.savefig(output_folder + '/' + field_A + '_vs_' + field_G +
                    '_I.png')
        plt.clf()

        # ========================================================
        # Create type II figure
        # ========================================================

        title_II = plot_dict['TITLE_II'][i]
        xlims_II = [float(plot_dict['XMIN_II'][i]),
                    float(plot_dict['XMAX_II'][i])]
        ylims_II = [float(plot_dict['YMIN_II'][i]),
                    float(plot_dict['YMAX_II'][i])]
        xlabel_II = plot_dict['XLABEL_II'][i]
        ylabel_II = plot_dict['YLABEL_II'][i]

        # Horizontal line on the figure, the ideal case
        min_temp = xlims_II[0]
        max_temp = xlims_II[1]
        xlims_id = [min_temp, max_temp]
        del(min_temp, max_temp)

        # Management of reversed axes (e.g. for magnitude)
        if xlims_id[0] > xlims_id[1]:
            min_temp = xlims_id[0]
            max_temp = xlims_id[1]
            xlims_id[0] = max_temp
            xlims_id[1] = min_temp
            del(min_temp, max_temp)

        ideal_x_II = np.arange(xlims_id[0], xlims_id[1], 0.01)
        ideal_y_II = 0 * np.arange(xlims_id[0], xlims_id[1], 0.01)

        # Data for type II figure
        data1_II = []
        data2_II = []
        legend_II = []
        quantiles = np.array([0.05, 0.95])

        # Definitions of variables in preparation for type III figure
        median_diff = []
        error_diff_p = []   # positive error relative to the median
        error_diff_m = []   # negative error relative to the median

        for j in range(m):
            data_temp = (data1[j] + data2[j]) / 2
            data1_II.append(data_temp)
            del(data_temp)

            data_temp = data2[j] - data1[j]
            data2_II.append(data_temp)

            median = np.median(data_temp)
            error = [dstat.get_quantile_errors(data_temp, quantiles)[0, 1],
                     dstat.get_quantile_errors(data_temp, quantiles)[0, 0]]
            del(data_temp)

            legend_fig = (str(cuts_G[m-1-j][0]) + ' < ' + field_G + ' < ' +
                          str(cuts_G[m-1-j][1]) + '/ median(95%): ' +
                          str('%.2e (%.2e,%.2e)'
                          % (median, error[0], error[1])))
            legend_II.append(legend_fig)
            del(legend_fig)

            # for type III plot
            median_diff.append(median)
            error_diff_p.append(error[0] - median)
            error_diff_m.append(median - error[1])
            del(median, error)

        # Plot and save type II figure
        legend_loc_II = 'upper left'
        pltemp.plot_data(data1_II, data2_II, formats, colors,
                         legend_II, legend_loc_II, xlabel_II, ylabel_II,
                         xlims_II, ylims_II, title_II)
        plt.plot(ideal_x_II, ideal_y_II, '--', color='k', label='')
        plt.xlim(xlims_II[0], xlims_II[1])
        plt.ylim(ylims_II[0], ylims_II[1])

        plt.savefig(output_folder + '/' + field_A + '_vs_' + field_G +
                    '_II.png')
        plt.clf()

        # ========================================================
        # Create type III figure
        # ========================================================

        title_III = plot_dict['TITLE_III'][i]
        xlims_III = [float(plot_dict['XMIN_III'][i]),
                     float(plot_dict['XMAX_III'][i])]
        ylims_III = [float(plot_dict['YMIN_III'][i]),
                     float(plot_dict['YMAX_III'][i])]
        xlabel_III = plot_dict['XLABEL_III'][i]
        ylabel_III = plot_dict['YLABEL_III'][i]

        # Horizontal line on the figure, the ideal case
        min_temp = xlims_III[0]
        max_temp = xlims_III[1]
        xlims_id = [min_temp, max_temp]
        del(min_temp, max_temp)
        flag_inv = 0

        # Management of reversed axes (e.g. for magnitude)
        if xlims_id[0] > xlims_id[1]:
            min_temp = xlims_id[0]
            max_temp = xlims_id[1]
            xlims_id[0] = max_temp
            xlims_id[1] = min_temp
            del(min_temp, max_temp)
            flag_inv = 1

        ideal_x_III = np.arange(xlims_id[0], xlims_id[1], 0.01)
        ideal_y_III = 0 * np.arange(xlims_id[0], xlims_id[1], 0.01)

        # Data for type III figure
        vec_G = []
        data1_III = dataG
        data2_III = []
        legend_III = []

        # Definition of the extremal limits of intervals for the variable G:
        # keep the points into the frame
        G_min = cuts_G[0][0]
        if math.isnan(G_min):
            G_min = xlims_III[0]
        G_max = cuts_G[m - 1][1]
        if math.isnan(G_max):
            G_max = xlims_III[1]

        for j in range(m):
            data_temp = (data2[j] - data1[j])
            data2_III.append(data_temp)
            del(data_temp)

            G_mean = 0
            if j == m-1:
                G_mean = G_min + (cuts_G[m-1-j][1] - G_min) / 2
                if G_mean < xlims_III[0]:
                    G_mean = xlims_id[0] + (cuts_G[m-1-j][1] - xlims_id[0]) / 2
            elif j == 0:
                G_mean = cuts_G[m-1-j][0] + (G_max - cuts_G[m-1-j][0]) / 2
                if G_mean > xlims_III[1]:
                    G_mean = cuts_G[m-1-j][0] + (xlims_id[1] - cuts_G[m-1-j][0]) / 2
            else:
                G_mean = cuts_G[m-1-j][0] + (cuts_G[m-1-j][1] - cuts_G[m-1-j][0]) / 2
            vec_G.append(G_mean)

            legend_fig = (str(cuts_G[m-1-j][0]) + ' < ' + field_G + ' < ' +
                          str(cuts_G[m-1-j][1]) + str(' (%i objects)'
                          % (len(data2_III[j]))))
            legend_III.append(legend_fig)
            del(legend_fig)

        # Plot and save type III figure
        legend_loc_III = 'lower left'
        pltemp.plot_data(data1_III, data2_III, formats, colors,
                         legend_III, legend_loc_III, xlabel_III, ylabel_III,
                         xlims_III, ylims_III, title_III)
        plt.plot(ideal_x_III, ideal_y_III, '--', color='k')
        ##### Replace this piece with the "box" code from matplotlib website
        plt.errorbar(vec_G, median_diff,
                     yerr=[error_diff_m, error_diff_p], fmt='-o',
                     color='k', zorder=10, elinewidth=2, capsize=8)
        ####
        plt.xlim(xlims_III[0], xlims_III[1])
        plt.ylim(ylims_III[0], ylims_III[1])
        plt.savefig(output_folder + '/' + field_A + '_vs_' + field_G +
                    '_III.png')
        plt.clf()

    return 1


# ---
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='''\
        Pipeline to perform analyses comparisons between two catalogs
        (catalog A and catalog B) through a set of plots, which detailed
        properties are given as input in the 'plots.txt' file.

        The user may choose fields of the catalogs A and B he wants to compare
        (FIELDNAME_A and FIELDNAME_B), and a field of the catalog A 
        (FIELDNAME_G) which will be used to realise some cuts.

        The input file 'plots.txt' has to follow a specific organization, an
        example can be found on the following page: 
        http://twiki.linea.gov.br/bin/view/CS82/FittingRobustness, and as 
        well some examples of plots.

        -----
        NB: an input file may be called with the prefixe '@', containing one
        argument per line.

        '''
        , fromfile_prefix_chars='@')

    # WARNING: As they are, the compulsory arguments do *not* take default 
    # values when not given! How do we use "default" in these cases?

    parser.add_argument(dest='cat_A',
                        default='./',
                        help="Path of Catalog A.")
    parser.add_argument(dest='cat_B',
                        default='./',
                        help="Path of Catalog B.")
    parser.add_argument(dest='output_folder',
                        default='./output',
                        help="Output folder for the plots.")
    parser.add_argument(dest='plots',
                        default='plots.txt',
                        help="Parameters file for the required plots.")
    parser.add_argument('-m','--matching',
                        dest='matching', 
                        default=0,
                        help="Perform matching of the two catalogs if != 0")
    parser.add_argument('-x','--exclude_multiple',
                        dest='multiple',
                        default=0,
                        help="Exclude cases where multiple matching occurred")
    parser.add_argument('-r','--matching_radius',
                        dest='radius',
                        default='0.0002777777777777778', #1 arcsecond
                        help="Matching radius (in degrees)")
    parser.add_argument('--ra_A',
                        dest='ra_A',
                        default='ALPHA_J2000', # CS82 ra coordinate name
                        help="Right Ascension column name in catalog A")
    parser.add_argument('--dec_A',
                        dest='dec_A',
                        default='DELTA_J2000', # CS82 dec coordinate name
                        help="Declination column name in catalog A")
    parser.add_argument('--ra_B',
                        dest='ra_B',
                        default='ALPHA_J2000', # CS82 ra coordinate name
                        help="Right Ascension column name in catalog B")
    parser.add_argument('--dec_B',
                        dest='dec_B',
                        default='DELTA_J2000', # CS82 dec coordinate name
                        help="Declination column name in catalog B")


    opts = parser.parse_args()

    coord_A_names = [opts.ra_A,opts.dec_A]
    coord_B_names = [opts.ra_B,opts.dec_B]

    out = main(opts.cat_A, opts.cat_B, opts.output_folder, opts.plots, 
                   bool(int(opts.matching)), bool(int(opts.multiple)), 
                   float(opts.radius), coord_A_names, coord_B_names)

    if out == None:
        parser.print_help()

    sys.exit(0)





    

    


