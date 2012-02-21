#!/usr/bin/env python
# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to perform 'statistical' analisis of lens inversion 
"""

from lenstool import *

from functions import *

from plot_utils import *

from math import sqrt

from gravlenspipeline import GravlensPipeline as GP

import numpy as np

def inversion_several_sub(plt_show = False, main_lens = []):
    print '\ntesting\n'
    if main_lens == []:
        main_lens = [{'z_lens': 1.0, 'x_centre': 0.0, 'y_centre': 0.0, \
                      'profil': 1, 'angle_pos': 35.0, 'ellipticite': 0.3, \
                      'v_disp': 500.0}]

    sub_lens = generate_substructure ( z_lens = 1.0, n_subs = 3, \
                                       vdisp_limit = 150 )
    sub_lens.append (main_lens[0])

    compute_curves_lenstool(sub_lens) #nedded to define the region where the
                                      #source will be placed
    courves_file = 'ce.dat'
    x_ca, y_ca = np.loadtxt( courves_file, usecols=(3, 4), unpack=True )

    ca_limit = 1.01*max([modulus(x_ca[i], y_ca[i]) for i in range(len(x_ca))])

    xy_im_sub = []
    tries = 0       #defining the maximum number of tries to source position to 
    max_tries = 200 #avoid stuck the program here
    while len(xy_im_sub) != 4 and tries != max_tries:
        source_position = generate_sources( n_src = 1, radius = ca_limit)
        if tries == max_tries-1:
            source_position = [[0, 0]]
        xy_im_sub = generate_arcs_lenstool(source = source_position, \
                                          lenses_vec = sub_lens, z_source = 2.0)
        tries = tries + 1
    
    
    
    optimized_parameters = make_inversion_lenstool(xy_im_sub, z_source = 2.0)


    if plt_show:
        xy_lens = []
        for i in range(len(sub_lens)):
            xy_lens.append( [sub_lens[i]['x_centre'], sub_lens[i]['y_centre']] )
            
        best_lens = [{'z_lens': 1.0, 'x_centre': 0.0, 'y_centre': 0.0, \
                      'profil': 1, \
                      'angle_pos': optimized_parameters['angle_pos'], \
                      'ellipticite': optimized_parameters['ellipticite'], \
                      'v_disp': optimized_parameters['v_disp']}]

        plot_points (source = source_position, images = xy_im_sub, marker = '^')
        compute_curves_lenstool (best_lens)
        plot_curves_lenstool (label = 'simple lens', plt_show = False, \
                              marker = 'b.')

        str_out = 'ell = %(ellipticite).2f, angle = %(angle_pos).2f, ' + \
                  'v_disp =  %(v_disp).2f \n'
        plt.title(str_out % optimized_parameters)
        
        plot_points(images = xy_lens)
        
        compute_curves_lenstool(sub_lens)
        plot_curves_lenstool(label = 'subs. lens', plt_show = plt_show, \
                             marker = 'm.')
    return optimized_parameters, sub_lens, tries
################################################################################
def generate_hist_sub():
    ell, angle_pos, v_disp = np.loadtxt \
                    ('tables/sub_f8_1000.dat', usecols = {1,2,3}, unpack = True)
    v_disp_subs = np.loadtxt \
                    ('tables/sub_f8_1000.dat', usecols = {4,5,6}, unpack = True)

    chi2 = np.loadtxt \
                    ('tables/sub_f8_1000.dat', usecols = {7}, unpack = True)

                    
    v_disp_reduced = 500 - np.sqrt( v_disp**2 - (\
                    v_disp_subs[0]**2 + v_disp_subs[1]**2 + v_disp_subs[2]**2) )
    v_disp = 500 - v_disp

#-------------------------------------------------------------------------------
    ell1, angle_pos1, v_disp1 = np.loadtxt \
                    ('tables/sub_f1_1000.dat', usecols = {1,2,3}, unpack = True)
    v_disp_subs1 = np.loadtxt \
                    ('tables/sub_f1_1000.dat', usecols = {4,5,6}, unpack = True)

    chi21 = np.loadtxt \
                    ('tables/sub_f1_1000.dat', usecols = {7}, unpack = True)

                    
    v_disp_reduced1 = 500 - np.sqrt( v_disp1**2 - (\
                 v_disp_subs1[0]**2 + v_disp_subs1[1]**2 + v_disp_subs1[2]**2) )
    v_disp1 = 500 - v_disp1
#-------------------------------------------------------------------------------
    ellM1, angle_posM1, v_dispM1 = np.loadtxt \
                   ('tables/sub_fM1_1000.dat', usecols = {1,2,3}, unpack = True)
    v_disp_subsM1 = np.loadtxt \
                   ('tables/sub_fM1_1000.dat', usecols = {4,5,6}, unpack = True)

   
    chi2M1 = np.loadtxt \
                    ('tables/sub_fM1_1000.dat', usecols = {7}, unpack = True)

                    
    v_disp_reducedM1 = 500 - np.sqrt( v_dispM1**2 -( \
               _disp_subsM1[0]**2 + v_disp_subsM1[1]**2 + v_disp_subsM1[2]**2) )
    v_dispM1 =  500 - v_dispM1
#-------------------------------------------------------------------------------
    hist_legend = ['forme 8', 'forme 1', 'forme -1']
    #plot_histogram ([v_disp_reduced, v_disp_reduced1, v_disp_reducedM1], \
    #                xlabel = '500 - \sigma_v^{\\rm red}', legend = hist_legend)
    #plot_histogram ([v_disp, v_disp1, v_dispM1], xlabel = '500 - \sigma_v^{}', \
    #                legend = hist_legend)
    #plot_histogram ([0.3 - ell, 0.3 - ell1, 0.3 - ellM1], \
    #                xlabel = '0.3 - \\varepsilon', legend = hist_legend)
    #plot_histogram ([35 - angle_pos, 35 - angle_pos1, 35 - angle_posM1], \
    #                xlabel = '35 - \\theta', legend = hist_legend)

    #plot_scatter_hist(ell, v_disp)
    #plot_scatter_hist(chi2, v_disp, xlabel='\chi^2', ylabel='500-\sigma_v')
    plot_scatter_hist(chi2, v_disp_reduced, xlabel='\chi^2', \
                      ylabel='500-\sigma_v^{red}')
    plot_scatter_hist(chi2, 0.3 - ell, xlabel='\chi^2', \
                      ylabel = '0.3 - \\varepsilon')

################################################################################
def for_inv_sub():
    file_out = open('test.txt','w')

    str_out = '# id ell_opt angle_opt sigma_opt sigma_sub1 xs1 ys1 sigma_sub2'
    str_out = str_out + ' xs2 ys2 sigma_sub3 xs3 ys3 xi2 tries\n'
    file_out.write( str_out )
    for i in range(1000):
        file_out.write( str(i) )
        out_inv = inversion_several_sub(plt_show = False)

        str_out = ' %(ellipticite)f %(angle_pos)f %(v_disp)f ' % out_inv[0]
        file_out.write(str_out)

        str_out = '%(v_disp)f ' % out_inv[1][0]
        str_out = str_out + '%(x_centre)f %(y_centre)f ' % out_inv[1][0]
        
        str_out = str_out + '%(v_disp)f ' % out_inv[1][1]
        str_out = str_out + '%(x_centre)f %(y_centre)f ' % out_inv[1][1]
        
        str_out = str_out + '%(v_disp)f ' % out_inv[1][2]
        str_out = str_out + '%(x_centre)f %(y_centre)f ' % out_inv[1][2]
        
        file_out.write( str_out )

        str_out = str(extract_parameter('best.par', 'Chi2pos:')[0][0])
        str_out = str_out + ' ' + str(out_inv[2]) + '\n'
        file_out.write( str_out )

    print out_inv[1]
    file_out.close()

def hist_sigma(file_name, plot, legend):

    obj_gp = GP()
    data = np.loadtxt(file_name, unpack = True)
    xy1 = [(data[5][i], data[6][i]) for i in range(len(data[5]))]
    xy2 = [(data[8][i], data[9][i]) for i in range(len(data[8]))]
    xy3 = [(data[11][i], data[12][i]) for i in range(len(data[11]))]

    re_main = obj_gp.sis_convert(0.0, 500, 1.0, 2.0)

    sigma_opt = data[3]
    sigma1 = data[4]
    sigma2 = data[7]
    sigma3 = data[10]

    sigma_tmp = 0.0
    sigma_red_final = []
    
    for i in range(len(sigma_opt)):
        #print sigma_opt[i], sigma1[i], sigma2[i], sigma3[i]#xy1[i], xy2[i], xy3[i]
        re_main = obj_gp.sis_convert(0.0, sigma_opt[i], 1.0, 2.0)
        sigma_tmp = sigma_opt[i]**2
        
        mod1 = sqrt(xy1[i][0]**2 + xy1[i][1]**2)
        if mod1 < re_main:
            #print '1', mod1
            sigma_tmp = sigma_tmp - sigma1[i]**2

        
        mod2 = sqrt(xy2[i][0]**2 + xy2[i][1]**2)
        if mod2 < re_main:
            #print '2', mod2
            sigma_tmp = sigma_tmp - sigma2[i]**2


        mod3 = sqrt(xy3[i][0]**2 + xy3[i][1]**2)
        if mod3 < re_main:
            #print '3', mod3
            sigma_tmp = sigma_tmp - sigma3[i]**2

        
        sigma_red_final.append( 500 - sqrt(sigma_tmp) )

    plot_histogram ( sigma_red_final, xlabel = '500 - \sigma_v^{red}',  \
                     legend = legend, plot=plot)
    return sigma_red_final
    #print sigma_red_final

################################################################################
# MAIN #########################################################################
################################################################################
if __name__ == '__main__':
    #generate_hist_sub()
    #for_inv_sub()
    hist_sigma('tables/sub_f8_pos.dat', plot=True, legend='forme 8')
    hist_sigma('tables/sub_fM1_pos.dat', plot=True, legend='forme -1')
    hist_sigma('tables/sub_f1_pos.dat', plot=True, legend='forme 1')
        
