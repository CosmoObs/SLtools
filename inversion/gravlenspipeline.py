# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to automate some aplocations of gravlens
"""

#from sltools.gravlens.find_cc import find_cc


from math import sqrt
from cosmology import CosmologyLcdm 
from functions import extract_parameter, modulus, extract_second_identifiers

import numpy as np
import matplotlib.pyplot as plt

import random

class GravlensPipeline():
    """
    Class to automate some aplocations of gravlens
    """
    def __init__(self, z_lens = 1.0, z_source = 2.0):
        self.cosmology = CosmologyLcdm()

        self.__z_lens = z_lens
        self.__z_src  = z_source
        #self.__sigmav_lenstool = sigmav_lenstool
        self.__par_file = ''

        self.__cosmology_model = {'type' : 'LCDM', \
                                  'Om0' : 0.3, 'Ol0' : 0.7, 'Ok0' : 0.0}
################################################################################
    def show_param(self):
        """
        Print SIEP parameters 'ell' and 'b'.

        Input:
         - NONE
        Output:
         - NONE
        """

        print '__ell_pot = '#, self.__ell_gravlens
        print '__b_sis = '#, self.__b_sis
        print '__z_lens = ', self.__z_lens
################################################################################
    def read_gravlens_src(self, mag_min = 0):
        """
        Read the sources and images from a 'findimg' gravlens output file

        Input:
         - mag_min float : minimum image absolute magnification 
        Output:
         - NONE
        """
        image_xy = []
        source_xy = 0
        amplification = []

        file_in = 'src_out.dat'
        file_inf = open(file_in, 'r')
        for line in file_inf:
            line_inf = line.split()


            if '#' and 'source' in line_inf:
                source_xy = [float(line_inf[0]), float(line_inf[1])]

            elif not('#' in line_inf) and len(line_inf)>0 and \
                               np.absolute(float(line_inf[2])) > mag_min + 1E-8:
                                   
                image_xy.append([float(line_inf[0]), float(line_inf[1])])
                amplification.append(float(line_inf[2]))

        dic = {'source_xy' : source_xy, 'image_xy' : image_xy, \
               'amplification' : amplification}

        file_inf.close()
        return dic
################################################################################
    def plot_curves(self, courves_file = 'crit.dat', plt_show = False):
        """
        Function to plot caustics and critical curves only from a gravlens
        output file
        Input:
         - courves_file  str : file that contains the curves (gravlens output)
         - plt_show     bool : if true show the plot in a window
        Output:
         - NONE
        """
        x1_c, y1_c, u1_c, v1_c = \
                     np.loadtxt(courves_file, usecols = (0, 1, 2, 3), \
                                unpack=True)
        plt.figure(1 , figsize=(16, 8))
        plt.subplot(1, 2, 1).set_aspect(1)
        plt.plot(u1_c, v1_c, 'b-', linewidth=3, label='gravlens')
        axis_limit = max( 1.075*max(np.absolute(u1_c)), \
                          1.075*max(np.absolute(v1_c)) )
        plt.axis([-axis_limit, axis_limit, -axis_limit, axis_limit])

        plt.subplot(1, 2, 2).set_aspect(1)
        plt.plot(x1_c, y1_c, 'b-', linewidth=3, label='gravlens')
        axis_limit = max( 1.075*max(np.absolute(x1_c)), \
                          1.075*max(np.absolute(y1_c)) )
        plt.axis([-axis_limit, axis_limit, -axis_limit, axis_limit])
        plt.legend()
        if plt_show:
            plt.show()
            plt.close()
################################################################################
    def plot_gravlens(self, source_xy, image_xy, show = True):
        """
        Plot images, source, critical courves and caustics

        Input:
         - source_xy   [] : source position
         - image_xy    [] : images postiion
         - show      bool : flag to show the plot
        Output:
         - NONE
        """
        crit_file = 'crit.dat'

        x1_c, y1_c, u1_c, v1_c = \
            np.loadtxt(crit_file, usecols = (0, 1, 2, 3), unpack=True)
        
        #x1_s, y1_s, u1_s, v1_s, x2_s, y2_s, u2_s, v2_s = \
        #    np.loadtxt(src_file, usecols={0, 1, 2, 3, 4, 5, 6, 7}, unpack=True)

        x_img = [row[0] for row in image_xy]
        y_img = [row[1] for row in image_xy]
        x_src = source_xy[0]
        y_src = source_xy[1]

        plt.figure(1, figsize=(16, 7))

        plt.subplot(1, 2, 1).set_aspect(1)
        plt.plot(u1_c, v1_c, 'b.', linewidth=3, label='gravlens')
        plt.scatter(x_src, y_src, s=50, color='darkblue', marker='+')
        axis_limit = max( 1.075*max(np.absolute(u1_c)),  \
                          1.075*max(np.absolute(v1_c)),  \
                          1.075*max(np.absolute([x_src, y_src])) )
        plt.axis([-axis_limit, axis_limit, -axis_limit, axis_limit])
        plt.title('Source Plane')


        plt.subplot(1, 2, 2).set_aspect(1)
        plt.plot(x1_c, y1_c, 'b.', linewidth=3, label='gravlens')
        plt.scatter(x_img, y_img, s=50, color='darkblue', marker='+')
        axis_limit = max( 1.075*max(np.absolute(x1_c)),  \
                          1.075*max(np.absolute(y1_c)),  \
                          1.075*max(np.absolute(x_img)), \
                          1.075*max(np.absolute(y_img)) )
        plt.axis([-axis_limit, axis_limit, -axis_limit, axis_limit])
        plt.title('Lens Plane')
        plt.legend()

        if show:
            plt.show()
################################################################################
    def ell_convert(self, ell_lenstool):
        """
        Convert ellipticity from lenstool definition to gravlens definitio 

        Input:
         - ell_lenstool float : lenstool mass distribution ellipticity
        Output:
         - float              : gravlens ellipticity
        """
        return 1.0 - sqrt( (1.0 - ell_lenstool)/(1.0 + ell_lenstool) )
################################################################################
    def sis_convert(self, ell_pot, sigmav_lenstool, z_lens, z_src):
        """
        Convert the 'mass paramenter' from lenstool definition to gravlens
        definitiom

        Input:
         - ell_pot         float : lenstool potential ellipticity
         - sigmav_lenstool float : velocity dispertion defined at lenstool
         - z_lens          float : lens redshift
         - z_source        float : source redshift
        Output:
         - float                 : SIEP 'b' paramenter (see 'alphapot' at table
                                   3.2 in gravlens manual)
        """
        const = 4.0*3.14159265*206264.806246799834988155
        const = const/(299792.458**2)

        dist_ratio = self.cosmology.int_flat(z_lens, z_src)/ \
                     self.cosmology.int_flat(0.0, z_src)

        value1 = dist_ratio*(sigmav_lenstool)**2.0 * const

        return value1 * sqrt(1.0 - ell_pot)
################################################################################
    def kappas_convert(self, ell_pot, sigmav_lenstool, r0_lenstool, \
                       z_lens, z_src):
        """
        Convert the 'mass paramenter' kappa_s from lenstool definition to
        gravlens definition

        Input:
         - ell_pot         float : lenstool potential ellipticity
         - sigmav_lenstool float : velocity dispertion defined at lenstool
         - r0_lenstool     float : lenstool scale parameter (NFW rs)
         - z_lens          float : lens redshift
         - z_source        float : source redshift
        Output:
         - float                 : kappa_s in the gravlens definition
        """
        const = 1.5*3.14159265*206264.806246799834988155
        const = const/(299792.458**2)

        const = const * sigmav_lenstool**2 * (1.0 - ell_pot) / r0_lenstool 
        
        dist_ratio = self.cosmology.int_flat(z_lens, z_src)/ \
                     self.cosmology.int_flat(0.0, z_src)

        return const * dist_ratio
################################################################################
    def rs_convert(self, r0_lenstool, ell_lenstool):
        """
        Convert the 'scale paramenter' r_s from lenstool definition to
        gravlens definitiom 

        Input:
         - ell_pot     float : lenstool potential ellipticity
         - r0_lenstool float : lenstool scale parameter (NFW rs)
        Output:
         - float             : rs in the gravlens deffinition
        """
        return r0_lenstool / sqrt( 1.0 - ell_lenstool )
################################################################################
    def convert(self, par_file, gravlens_file, xy_src, generate_crit = True):
        """
        Convert a .par file to a gravlens input file. Only to generate point
        images in the SIEP model

        Input:
         - par_file       str: imput .par file
         - gravlens_file  str: output file to be used in gravlens
         - xy_src        [,] : source position
        Output:
         - NONE
        """
    #convert lenstool sintax to gravlens sintax
        self.__par_file = par_file
        profil = \
        extract_parameter(file_in = self.__par_file, identifier = 'profil')

        x_centre = \
        extract_parameter(file_in = self.__par_file, identifier = 'x_centre')

        y_centre = \
        extract_parameter(file_in = self.__par_file, identifier = 'y_centre')

        ellipticity = \
        extract_parameter(file_in = self.__par_file, identifier = 'ellipticite')

        angle_pos = \
        extract_parameter(file_in = self.__par_file, identifier = 'angle_pos')

        v_disp = \
        extract_parameter(file_in = self.__par_file, identifier = 'v_disp')

        z_lens = \
        extract_parameter(file_in = self.__par_file, identifier = 'z_lens')

        nplan = \
        extract_parameter(file_in = self.__par_file, identifier = 'nplan')

        nlens = \
        extract_parameter(file_in = self.__par_file, identifier = 'nlens')
        nlens = int(nlens[0][0])

        print nplan
        
        gravlens_profil = []
        gravlens_ellipticity = []
        gravlens_masspar = []


        file_opened = open(gravlens_file, 'w')

        file_opened.write('set omega 0.3\n')
        file_opened.write('set lambda 0.7\n')
        file_opened.write('set hval 0.7\n')
        file_opened.write('set zlens 1.0\n')
        file_opened.write('set zsrc 2.0\n')

        file_opened.write('set crittol 1.0e-10')

        file_opened.write('gridmode 1\n')
        file_opened.write('set gridlo1 0.0\n')
        file_opened.write('set gridhi1 10.0\n')

        file_opened.write('set ngrid1 200.0\n')
        file_opened.write('set ngrid2 200.0\n')

        file_opened.write('set maxlevel 3\n')

        str_out = 'startup ' + str(nlens) + ' 1\n'
        file_opened.write(str_out)

        for i in range(nlens):
            #convert potential name
            if profil[i][0] == '1':
                gravlens_profil.append('alphapot')

            #convert ellipticity
            gravlens_ellipticity.append(\
                                 self.ell_convert(float(ellipticity[i][0])/3.0))

            #convert mass parameter
            gravlens_masspar.append(\
                  self.sis_convert( float(ellipticity[i][0])/3.0 , \
                  float(v_disp[i][0]), float(z_lens[i][0]), float(nplan[0][1])))

            str_out = '    ' + str(gravlens_profil[i]) + ' ' + \
                      str(gravlens_masspar[i])
            str_out = str_out + ' ' + x_centre[i][0] + ' ' + y_centre[i][0]
            str_out = str_out + ' ' + str(gravlens_ellipticity[i]) + ' '
            str_out = str_out + str( float(angle_pos[i][0]) - 90.0 ) + ' '
            str_out = str_out + '0.0 0.0 0.0 0.0 1.0\n'

            file_opened.write(str_out)

        for i in range(nlens):
            file_opened.write('    0 0 0 0 0 0 0 0 0 0\n')

        if generate_crit:
            file_opened.write('plotcrit crit.dat\n')

        str_out = 'findimg ' + str(xy_src[0]) + ' ' + str(xy_src[1]) + \
                    ' src_out.dat\n'
        file_opened.write(str_out)

        file_opened.close()
################################################################################
    def convert2(self, par_file, gravlens_file = 'test.gravlens', \
                 generate_crit = True):
        all_identifiers = extract_second_identifiers(par_file, 'potential')

        nplan = \
        extract_parameter(par_file, identifier = 'nplan')
        z_src = float(nplan[0][1])

        nlens = \
        extract_parameter(file_in = par_file, identifier = 'nlens')
        nlens = int(nlens[0][0])

        #-----------------------------------------------------------------------
        file_opened = open(gravlens_file, 'w')

        file_opened.write('set omega 0.3\n')
        file_opened.write('set lambda 0.7\n')
        file_opened.write('set hval 0.7\n')
        file_opened.write('set zlens 1.0\n')
        file_opened.write('set zsrc 2.0\n')

        file_opened.write('set crittol 1.0e-10')

        file_opened.write('gridmode 1\n')
        file_opened.write('set gridlo1 0.0\n')
        file_opened.write('set gridhi1 10.0\n')

        file_opened.write('set ngrid1 200.0\n')
        file_opened.write('set ngrid2 200.0\n')

        file_opened.write('set maxlevel 3\n')

        str_out = 'startup ' + str(nlens) + ' 1\n'
        file_opened.write(str_out)
        #-----------------------------------------------------------------------



        
        for i in range(nlens):
            #print i
            #print all_identifiers[i]['profil']
            z_lens = float (all_identifiers[i]['z_lens'][0])

            if all_identifiers[i]['profil'][0] == '12': #NFW case
                r0_lenstool = float (all_identifiers[i]['core_radius'][0])
                ellipticite = float (all_identifiers[i]['ellipticite'][0])
                sigmav_lenstool = float (all_identifiers[i]['v_disp'][0])
                
                #converting the quantities
                rs_gl = self.rs_convert(r0_lenstool, ellipticite)
                ks_gl = self.kappas_convert(ellipticite, sigmav_lenstool, r0_lenstool, \
                                       z_lens, z_src)
                ell_gl = self.ell_convert(ellipticite)
                
                #making string with gravlens sintax
                str_out = '    nfwpot ' + str(ks_gl) + ' '
                str_out = str_out + str(all_identifiers[i]['x_centre'][0]) + ' '
                str_out = str_out + str(all_identifiers[i]['y_centre'][0]) + ' '
                str_out = str_out + str(ell_gl) + ' '
                str_out = str_out + str(\
                          float(all_identifiers[i]['angle_pos'][0]) - 90) + ' '
                str_out = str_out + '0.0 0.0 ' + str(rs_gl) + ' 0.0 0.0\n' 
                
                #print str_out
                file_opened.write(str_out)
                #print rs_gl, ks_gl, ell_gl

                
                
            if all_identifiers[i]['profil'][0] == '1': #SIS case
                ellipticite = float (all_identifiers[i]['ellipticite'][0])
                sigmav_lenstool = float (all_identifiers[i]['v_disp'][0])

                
                #converting the quantities
                ell_gl = self.ell_convert( ellipticite/3.0 )
                sis_par = self.sis_convert( ellipticite/3.0, sigmav_lenstool, \
                                            z_lens, z_src )
                
                #making string with gravlens sintax
                str_out = '    alphapot ' + str(sis_par) + ' '
                str_out = str_out + str(all_identifiers[i]['x_centre'][0]) + ' '
                str_out = str_out + str(all_identifiers[i]['y_centre'][0]) + ' '
                str_out = str_out + str(ell_gl) + ' '
                str_out = str_out + str(\
                          float(all_identifiers[i]['angle_pos'][0]) - 90) + ' '
                str_out = str_out + '0.0 0.0 0.0 0.0 1.0\n' 
                
                #print str_out
                file_opened.write(str_out)
                #print sis_par, ell_gl


        for i in range(nlens):
            file_opened.write('    0 0 0 0 0 0 0 0 0 0\n')

        if generate_crit:
            file_opened.write('plotcrit crit.dat\n')
            #print '\n'
            
        str_out = 'findimg ' + '1.347' + ' ' + '-0.093' + \
                    ' src_out.dat\n'
        file_opened.write(str_out)

        str_out = 'plotdef1 pot_gl.txt -5 5 200 -5 5 200\n'
        file_opened.write(str_out)

        str_out = 'plotkappa mass_gl.fits 3 -5 5 200 -5 5 200'
        file_opened.write(str_out)

        file_opened.close()
