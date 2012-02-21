# ==================================
# Author: Gabriel Bartosch Caminha - gbcaminha@gmail.com
# ==================================
import pyfits

import numpy as np

from plot_utils import plot_scatter_hist

def create_region_file(ra, dec, file_out, text = []):
    """
    Generate a region file named 'file_out' from a list of 'ra' and 'dec'.
    Circles will be placed for each (RA,DEC) with coments from 'text'.

    Input:
     - ra        [, ...] : list of RA
     - dec       [, ...] : list of DEC
     - file_out      str : file name (.reg) for output
     - text      [, ...] : list with information to be added for each (RA,DEC) 
    Output:
     - NONE
    """
    #file_out = '/home/gcaminha/Downloads/obj.reg'
    file_opened = open(file_out, 'w')
    str_out = 'global color=green dashlist=8 3 width=1 font="helvetica 10 '
    str_out = str_out + 'normal roman" select=1 highlite=1 dash=0 fixed=0 '
    str_out = str_out + 'edit=1 move=1 delete=1 include=1 source=1\n'
    file_opened.write(str_out)

#    file_opened.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    file_opened.write('fk5\n')

    radius = 2.5

    if len(text) == 0:
        for i in range(len(ra)):
            str_out = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + \
                      str(radius) + '")\n'
            file_opened.write(str_out)
    else:
        for i in range(len(ra)):
            str_out = 'circle(' + str(ra[i]) + ',' + str(dec[i]) + ',' + \
                      str(radius) + '") # text={' + str(text[i]) + '}\n'
            file_opened.write(str_out)
    file_opened.close()

def objects_region():
    #Program to select objtects from a list in a sky region and make a .reg file
    #for all objects within the region.
    image_path = '/home/gcaminha/Downloads/mrgS20111116S0063_add.fits'
    im_hdulist = pyfits.open(image_path)
    im_header = im_hdulist[1].header

    table_path = '/home/gcaminha/Gabriel/BKP-2011-10-25/catalogues/galSpecInfo-dr8-S82.fits'
    tb_hdulist = pyfits.open(table_path)
    tb_data = tb_hdulist[1].data

    #plot_scatter_hist(tb_data.field('DEC'), tb_data.field('Z'))

    ra_max = 10.3422018
    ra_min = 10.2327202

    dec_min = -0.78519495
    dec_max = -0.67571335

    mask = tb_data.field('RA') > ra_min
    tb_data = tb_data[mask]
    mask = tb_data.field('RA') <  ra_max
    tb_data = tb_data[mask]

    mask = tb_data.field('DEC') > -0.78519495
    tb_data = tb_data[mask]
    mask = tb_data.field('DEC') < -0.67571335
    tb_data = tb_data[mask]

    #print tb_data['RA'], tb_data['DEC']
    print tb_data.names
    print tb_data['V_DISP']

    create_region_file(tb_data['RA'], tb_data['DEC'], \
                       file_out = '/home/gcaminha/Downloads/obj_spec_tmp.reg', \
                       text = tb_data['Z'])

    im_hdulist.close()
    tb_hdulist.close()


    
#MAIN PROGRAM

#file_path = '/home/gcaminha/Downloads/Gemini_196_3_gbcaminha_0.csv'
file_path = '/home/gcaminha/Downloads/spAll-v5_4_45.dat'

#ra, dec, g_mag, r_mag = np.loadtxt(file_path, usecols={1, 3, 7, 8}, unpack=True, delimiter=',')

ra, dec, z = np.loadtxt(file_path, usecols={5, 8, 9}, unpack=True, skiprows=2)

ra_max = 10.3422018
ra_min = 10.2327202

dec_min = -0.78519495
dec_max = -0.67571335

ra_final = []
dec_final = []
z_final = []

for i in range(len(z)):
    if ra[i] > ra_min and ra[i] < ra_max and dec[i] > dec_min and dec[i] < dec_max:
        print ra[i], dec[i], z[i]
        ra_final.append(ra[i])
        dec_final.append(dec[i])
        z_final.append(z[i])


#create_region_file(ra_final, dec_final, \
#                   file_out = '/home/gcaminha/Downloads/obj_spec_dr9.reg', \
#                   text = z_final)
#for i in range(len(ra)):
#    print ra[i], dec[i], g_mag[i], r_mag[i]

#objects_region()

















