#!/usr/bin/env python

from plot_utils import plot_histogram, plot_scatter_hist, plot_mult_hist

from math import log, pi, sqrt

import math

from functions import extract_parameter, strip_list

from cosmology import CosmologyLcdm as Cosm

import numpy as np

import constants as cn

import os
################################################################################
def list_files(path = '.', file_extension = '.'):
    """
    Function to list all files with 'file_extension' in a 'path'

    The function selec the files names that ends with 'file_extension'

    Input:
     - path           str: path to directory where you whant to select the files
     - file_extension str: string with the extension to select the files
    Output:
     - list with file names
    """    
    dirList=os.listdir(path)
    files_out = []
    for fname in dirList:
        if fname[len(fname)-len(file_extension):] == file_extension:
            files_out.append(fname)
    return files_out

################################################################################
def plot_properties():
    statinfo_path = 'satellites'
    sub_files = list_files(path=statinfo_path, file_extension = '.sub');
    sub_files.sort()

    mass_mult = []
    dist_mult = []
    x_mult = []
    y_mult = []
    conc_mult = []

    for i in range(len(sub_files)):
        statinfo_path = 'satellites/' + sub_files[i]
        mass, dist, x, y, conc = np.loadtxt(statinfo_path, unpack=True, \
                                            usecols={0, 1, 2, 3, 4})
        mass_mult.append(mass)
        dist_mult.append(dist)
        x_mult.append(x)
        y_mult.append(y)
        conc_mult.append(conc)
        
    print len(sub_files)
    print len(x_mult)

    plot_histogram(conc_mult, plot=True, xlabel = 'conc')
    plot_histogram(dist_mult, plot=True, xlabel = 'dist')
    plot_histogram(mass_mult, plot=True, xlabel = 'mass')
    for i in range(len(x_mult)):
        plot_scatter_hist(x_mult[i], y_mult[i])

################################################################################
def v_disp(mass_vir, R_vir, conc_vir):
    """
    Function to convert MOKA outputs to velocity dispersion (lentool)

    Input:
     - mass_vir : halo virial mass [M_{sol}]
     - R_vir    : halo virial radius [Mega Parsec]
     - conc_vir : halo concentration parameter [R_vir / r_s]
    Output:
     - v_disp [km / s]
    """    
    f_conc = log(1.0 + conc_vir)/conc_vir - 1.0/(1.0 + conc_vir)
    const = 2.0/3.0 / pi / f_conc

    mass_solar = 1.98892E30    #[kilo-grams]
    megaparsec = 3.08568025E19 #[kilo-meters]
    grav_const = 6.67384e-11   #[meters^3 kilograms^-1 seconds^-2]

    #const = const*grav_const

    m_vir_sol = mass_vir * mass_solar
    R_vir_km = R_vir * megaparsec

    g_final = grav_const*mass_solar/1E9
    return sqrt(g_final*mass_vir/R_vir_km)*const

def rs(R_vir, conc_vir):
    return R_vir/conc_vir
################################################################################
def mass_vir(z_in, r_vir, cosmology):
    delta_vir = cosmology.virial_overdensity(z_in)
    omega_mz = cosmology.omega_m(z_in)
    omega_m0 = cosmology.omega_m(0)





    print delta_vir, omega_mz, omega_m0
    return 0.0
################################################################################
def subs():
    path_satinfo = 'satellites/satinfo.0.sub'
    np.loadtxt( path_satinfo, unpack=True, usecols=(1,2) )
    subhalo_mass, x_projec, y_projec, conc
    return 0.0







################################################################################
if __name__ == '__main__':
    Rvir0 = extract_parameter('moka.out', 'Rvir0 =')
    phy_size = extract_parameter('moka.out', 'physical size of the region:')
    axes_3d = extract_parameter('moka.out', 'initial (a, b, c) =')
    mass = extract_parameter('moka.out', 'with mass: ')
    conc, z_lens = np.loadtxt( 'info_haloes.dat', unpack=True, usecols=(1,2) )

    radian = 206264.806247 #arcsec
    cos_obj = Cosm()
    for i in range(len(Rvir0)):
        print '------------------------------------------------------------'
        print 'LENS ', i,':'
        print 'Rvir0 = ', Rvir0[i][0], "mpc"
        print 'Rvir0 = ', float(Rvir0[i][0])/cos_obj.ang_dis(0.0, z_lens[i])*radian, "arcsec"
        print 'phy_size = ', phy_size[i][0]

        axes_3d[i] = map( float, strip_list(axes_3d[i], ',)( ') )
        print 'axes = ', axes_3d[i]
        print 'mass = ', float(mass[i][0])
        print 'conc = ', conc[i]
        print 'z_lens = ', z_lens[i]
        print 'v_disp = ', v_disp(float(mass[i][0]), float(Rvir0[i][0]), conc[i])
        print 'r_s = ', rs(float(Rvir0[i][0]), conc[i])
        print '------------------------------------------------------------'
    #mass_vir(z_lens[5], float(Rvir0[5][0]), cos_obj)
    #print cos_obj.critical_density(0.0)
    #print cos_obj.mass_virial(1.10989, 0.25)
    print cos_obj.ang_dis(0, 0.35) * math.sin(1.0/cn.radians) * 1E3







    
