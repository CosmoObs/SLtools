# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
"""
Package to deal with composite lens models
"""
from scipy import integrate

import numpy as np
import math

import matplotlib.pyplot as plt

import constants as cs
import cosmology

import functions as fu
import gravlenspipeline as gp
import pseudo_elliptical_models as pem
import radial_models as rm
import plot_utils as pu

class SisComp(object):
    def __init__(self, par_file, z_src_in):
        #print 'SisComp: Initializing SisComp'
        self.lens = fu.extract_second_identifiers(par_file, "potential")
        if len(self.lens) == 0:
            #print "SisComp: warning: using identifier potentiel instead" \
            #      "potential" 
            self.lens = fu.extract_second_identifiers(par_file, \
                                                              "potentiel")
        self.z_src = z_src_in
                     #float(fu.extract_parameter(par_file, 'nplan')[0][1])

        for j in range(len(self.lens)):
            for i in self.lens[j].keys():
                self.lens[j][i][0] = float( self.lens[j][i][0] )

        #print "    SisComp: Using file:", par_file
        #print "    that contains", len(self.lens), "lenses"
        self._read_par_file()
################################################################################
    def _read_par_file(self):
        sis_convert = gp.GravlensPipeline().sis_convert

        self.ell_sis = []
        self.v_sis = []
        self.z_lens = []
        self.angle_pos = []

        self.sis_par = []

        self.x_0 = []
        self.y_0 = []

        count = 0
        for i in self.lens:
            self.x_0.append( (i["x_centre"][0]) )
            self.y_0.append( (i["y_centre"][0]) )

            self.ell_sis.append( (i["ellipticite"][0]) )
            self.v_sis.append( (i["v_disp"][0]) )
            self.z_lens.append( (i["z_lens"][0]) )
            self.angle_pos.append( i["angle_pos"][0] )

            self.sis_par.append( sis_convert(0.0, self.v_sis[count], \
                            self.z_lens[count], self.z_src) )
            count = count + 1
################################################################################
    def kappa_composite(self, x_in, y_in):
        conv_tot = 0.0
        xy_rot = []
        psis = pem.kappa_sigma_test
        for i in range(len(self.ell_sis)):
            xy_rot.append( fu.rotate_vector( [x_in - self.x_0[i], \
                           y_in - self.y_0[i]], -self.angle_pos[i]) )
            x_eff = xy_rot[i][0]
            y_eff = xy_rot[i][1]
            conv_tot = conv_tot + psis( x_eff, y_eff, self.ell_sis[i], \
                                        self.sis_par[i])
        return conv_tot
################################################################################
    def plot_proj_kappa(self, plt_show = False):
        pu.plot_contour_func(self.kappa_composite, -1.1, 1.1, 75, -1.1, 1.1, 75)
        if plt_show:
            plt.show()
################################################################################
    def int_proj_kappa(self, radius):

        npts_r = 30
        npts_t = 30
        r_h = radius / float(npts_r)
        t_h = math.pi / float(npts_t)

        r_vec = np.linspace(0 , radius - r_h, npts_r)
        t_vec = np.linspace(0 , 2.0*math.pi - t_h, npts_t)

        x_coor = 0.0
        y_coor = 0.0

        out = 0.0
        for i_r in r_vec:
            for i_t in t_vec:
                x_coor = (i_r + r_h) * math.cos(i_t + t_h)
                y_coor = (i_r + r_h) * math.sin(i_t + t_h)
                out = out + r_h*t_h*self.kappa_composite(x_coor, y_coor) \
                      *(i_r + r_h)
        return out*2.0
################################################################################
    def proj2sigmav(self, int_kappa, z_lens, z_src):

        dos = cosmology.CosmologyLcdm().ang_dis(0.0, z_src)
        dls = cosmology.CosmologyLcdm().ang_dis(z_lens, z_src)

        sis_par_eff = math.sqrt(int_kappa/math.pi)

        out = sis_par_eff*cs.c2*dos/dls/4.0/math.pi

        #return math.pi*radius*self.sis_par[0]
        return math.sqrt(out/206264.806246799834988155)

################################################################################
################################################################################

def composing_sis_test(x_in, y_in, lens, z_src):
    print "func composing_sis not working properly - use composing_sis"
    x_eff = 0;
    y_eff = 0;

    conv_tot = 0.0

    psis = pem.PellModels( rm.SisLens() )

    sis_convert = gp.GravlensPipeline().sis_convert

    ell_sis = []
    v_sis = []
    z_lens = []
    angle_pos = []

    a_sis = []
    b_sis = []
    sis_par = []

    x_0 = []
    y_0 = []

    xy_rot = []

    count = 0
    for i in lens:
        x_0.append( (i["x_centre"][0]) )
        y_0.append( (i["y_centre"][0]) )

        ell_sis.append( (i["ellipticite"][0]) )
        v_sis.append( (i["v_disp"][0]) )
        z_lens.append( (i["z_lens"][0]) )
        angle_pos.append(i["angle_pos"][0])


        ell_sis[count] = 1.0 - \
                  math.sqrt( ( 1.0 - ell_sis[count] )/( 1.0 + ell_sis[count] ) )

        a_sis.append(1.0 - ell_sis[count])
        b_sis.append(1.0 + ell_sis[count])

        sis_par.append( sis_convert(ell_sis[count], v_sis[count], \
                        z_lens[count], z_src) )#/math.sqrt(1.0 - ell_sis[count]) )

        count = count + 1

    for i in range(len(ell_sis)):
        xy_rot.append( fu.rotate_vector( [x_in - x_0[i], y_in - y_0[i]], \
                                         - angle_pos[i] - 90.0) )
        x_eff = xy_rot[i][0]# - x_0[i]
        y_eff = xy_rot[i][1]# - y_0[i]
        conv_tot = conv_tot + psis.kappa_pell( x_eff, y_eff, a_sis[i], \
                                               b_sis[i], sis_par[i])
    return conv_tot
################################################################################
def composing_sis(x_in, y_in, lens, z_src):

    x_eff = 0;
    y_eff = 0;

    conv_tot = 0.0

    psis = pem.kappa_sigma_test

    sis_convert = gp.GravlensPipeline().sis_convert

    ell_sis = []
    v_sis = []
    z_lens = []
    angle_pos = []

    a_sis = []
    b_sis = []
    sis_par = []

    x_0 = []
    y_0 = []

    xy_rot = []

    count = 0
    for i in lens:
        x_0.append( (i["x_centre"][0]) )
        y_0.append( (i["y_centre"][0]) )

        ell_sis.append( (i["ellipticite"][0]) )
        v_sis.append( (i["v_disp"][0]) )
        z_lens.append( (i["z_lens"][0]) )
        angle_pos.append(i["angle_pos"][0])


        sis_par.append( sis_convert(0.0, v_sis[count], \
                        z_lens[count], z_src) )#/math.sqrt(1.0 - ell_sis[count]) )

        count = count + 1

    for i in range(len(ell_sis)):
        xy_rot.append( fu.rotate_vector( [x_in - x_0[i], y_in - y_0[i]], \
                                         - angle_pos[i]) )
        x_eff = xy_rot[i][0]# - x_0[i]
        y_eff = xy_rot[i][1]# - y_0[i]
        conv_tot = conv_tot + psis( x_eff, y_eff, ell_sis[i], sis_par[i])
    return conv_tot
################################################################################
def arg_mass_polar(r_in, t_in, lens, z_src):
    x_coor = r_in * math.cos(t_in)
    y_coor = r_in * math.sin(t_in)
    return composing_sis(x_coor, y_coor, lens, z_src)*r_in
################################################################################
def int_arg_mass_polar(t_in, lens, z_src, radius):
    out = integrate.quad(arg_mass_polar, a = 0.0, b = radius, \
                         args = (t_in, lens, z_src))
    return out[0]
################################################################################
def proj_mass_comp(lens, z_src, radius):
    out = integrate.quad(int_arg_mass_polar, a = 0.0, b = 2.0*math.pi, \
                         args = (lens, z_src, radius))
    return out
################################################################################
def proj_mass_comp2(lens, z_src, radius):
    npts_r = 100
    npts_t = 100
    r_h = radius / float(npts_r)
    t_h = math.pi / float(npts_t)

    r_vec = np.linspace(0 , radius - r_h, npts_r)
    t_vec = np.linspace(0 , 2.0*math.pi - t_h, npts_t)

    out = 0.0
    for i_r in r_vec:
        for i_t in t_vec:
            out = out + r_h*t_h*arg_mass_polar(i_r+ r_h, i_t + t_h, lens, z_src)

    return out*2.0
################################################################################
def proj_mass_comp3(lens, z_src, radius):
    npts_x = 100
    npts_y = 100
    x_h = radius / float(npts_x)
    y_h = radius / float(npts_y)
    x_vec = np.linspace(-radius + x_h/2.0, radius - x_h, npts_x)
    y_vec = np.linspace(-radius + y_h/2.0, radius - y_h, npts_y)

    out = 0.0

    vec_plt = [[],[]]

    for i_x in x_vec:
        for i_y in y_vec:
            if i_x**2.0 + i_y**2.0 < radius**2.0:
                out = out + x_h*y_h*composing_sis(i_x, i_y, lens, z_src)
    return out*4.0




