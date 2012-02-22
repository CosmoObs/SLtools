from __future__ import division
from ast import literal_eval as lit
import sys
import os
import errno
import numpy as np
import math
import cmath
from scipy import optimize
import collections
import fit_mag as fmg
import glob
import pyfits
import matplotlib.pyplot as pl
import plot_templates as pltemp
from time import strftime


sys.path.append("/home/brunomor/Documents/Trabalho/Repositorio") ### Append Bruno
from sltools.catalog import fits_data as fd ### Import Bruno
from sltools.catalog import ascii_data as ascd ### Import Bruno
from sltools.catalog import table_matching as tbm
from sltools.coordinate import wcs_conversion as wcscnv

########################################
### SETTING FOLDERS, FILES, ETC   ######
########################################

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno == errno.EEXIST:
            pass
        else: raise


input_cs82 = sys.argv[1]

sdss_csv = sys.argv[2]

folder_path = os.path.dirname(input_cs82).replace('/AFTER_PSFEX','')

filename = os.path.basename(input_cs82)

mkdir_p(folder_path + '/plots_sdss/')

# Choose plot fields:

if filename.find('deV') > 0:

    fieldnames_sdss = ['deVRad_i','deVAB_i','deVPhi_i']
    fieldnames_cs82 = ['SPHEROID_REFF_WORLD','SPHEROID_ASPECT_WORLD','THETAMODEL_IMAGE']

elif filename.find('exp') > 0:

    fieldnames_sdss = ['expRad_i','expAB_i','expPhi_i']
    fieldnames_cs82 = ['DISK_SCALE_WORLD','DISK_ASPECT_WORLD','DISK_THETA_IMAGE']

else:

    print "There's something wrong with the input file"
    sys.exit(1)


######################################################################################
### IMPORTING THE DATA AND DOING THE NECESSARY MATCHING, BASIC CUTS AND REORDERING ###
######################################################################################

# ===================== Transformation from the SDSS skyserver CSV to a FITS catalog


sdss_csv = sys.argv[2]

fieldnames=['objID','ra','dec','psffwhm_i','deVRad_i','deVRadErr_i','deVAB_i','deVABErr_i','deVPhi_i','deVMag_i','deVMagErr_i','deVFlux_i','deVFluxIvar_i',
'expRad_i','expRadErr_i','expAB_i','expABErr_i','expPhi_i','expMag_i','expMagErr_i','expFlux_i','expFluxIvar_i','cModelMag_i','cModelMagErr_i',
'cModelFlux_i','cModelFluxIvar_i','modelMagErr_i','modelFlux_i','modelFluxIvar_i','run','rerun','camcol','field','type','modelMag_u','modelMag_g',
'modelMag_r','modelMag_i','modelMag_z']

sdss_dict = ascd.dict_from_csv(sdss_csv,fieldnames, header_lines=1, delimiter=',',dialect='excel')

sdss_tbhdu = fd.dict_to_tbHDU(sdss_dict)

sdss_data = sdss_tbhdu.data


# =============================== Getting CS82 Catalog Information

  
cs82_tbhdu = pyfits.open(input_cs82,ignore_missing_end=True,memmap=True)[2]

cs82_cut_hdu = fd.sample_entries(cs82_tbhdu,CLASS_STAR=(0,0.95), MAG_AUTO=(0,21), MAGERR_AUTO=(0,0.2172))

cs82_data = cs82_cut_hdu.data

print len(cs82_data)



# ============================= Doing a Nearest Neighbour Match and Getting Full Reordered cs82_data.

sdss_ra = sdss_data.field('ra').real
sdss_dec = sdss_data.field('dec').real
sdss_radec = zip(sdss_ra,sdss_dec)

cs82_ra = cs82_data.field('alpha_J2000')
cs82_dec = cs82_data.field('delta_J2000')
cs82_radec = zip(cs82_ra, cs82_dec)


near_neigh = np.array(tbm.nearest_neighbour(sdss_radec,cs82_radec))

radius = 0.000277 # degrees ~ 1 arc-second

true_match = np.array([ float(_d) < float(radius) for _d in near_neigh[:,1]])

sdss_data = sdss_data[true_match]
near_neigh = near_neigh[true_match]

mask = near_neigh[:,0].astype(int)

cs82_full = cs82_data[mask]

# I haven't made sure that the match was unique from the cs82 side, 
# i.e. two different sdss objects may be matching the same cs82 one.
# But in view of the cut-off radius chosen, this is a completely secondary 
# effect.

'''
sdss_ra = sdss_data.field('ra').real
sdss_dec = sdss_data.field('dec').real
sdss_radec = zip(sdss_ra,sdss_dec)

cs82_ra = cs82_full.field('alpha_J2000')
cs82_dec = cs82_full.field('delta_J2000')
cs82_radec = zip(cs82_ra, cs82_dec)

near_neigh3 = np.array(tbm.nearest_neighbour(cs82_radec,sdss_radec))


print near_neigh3[:,1].max()
'''
print len(cs82_full)

# Keeping only FLAGS = 0 for the analysis

mask_flags = (cs82_full.field('FLAGS') == 0)

cs82_full = cs82_full[mask_flags]
sdss_data = sdss_data[mask_flags]

print len(sdss_data)

'''
# Another ds9 mask

image_name = '/home/brunomor/Documents/Trabalho/S82p28m_y.V2.7A.swarp.cut.fits'

hdu_image = pyfits.open(image_name)

image_header = hdu_image[0].header

sdss_ra = sdss_data.field('ra').real
sdss_dec = sdss_data.field('dec').real

sdss_ra_list = sdss_ra.tolist()
sdss_dec_list = sdss_dec.tolist()

pixcrd = []
xlist = []
ylist = []

for i in range(len(sdss_ra)):
    pixcrd.append(wcscnv.radec2xy(image_header,sdss_ra_list[i],sdss_dec_list[i]))
    xlist.append(pixcrd[i][0])
    ylist.append(pixcrd[i][1])

ascd.write_ds9cat(xlist,ylist,30,'circle','yellow', outputfile='boss_new_matching.reg')


# Another ds9 mask


cs82_ra = cs82_full.field('alpha_J2000').real
cs82_dec = cs82_full.field('delta_J2000').real

cs82_ra_list = cs82_ra.tolist()
cs82_dec_list = cs82_dec.tolist()

pixcrd = []
xlist = []
ylist = []

for i in range(len(sdss_ra)):
    pixcrd.append(wcscnv.radec2xy(image_header,cs82_ra_list[i],cs82_dec_list[i]))
    xlist.append(pixcrd[i][0])
    ylist.append(pixcrd[i][1])

ascd.write_ds9cat(xlist,ylist,40,'circle','green',outputfile='cs82_new_matching.reg')
'''

# Selecting the fields to be used for the comparison

# bla bla

################################################
##########    SEVERAL CUTS & PLOTS    ##########
################################################

# Color cuts 

cs82_20_to_21 = cs82_full[(20 <= cs82_full.field('MAG_AUTO'))]

cs82_19_to_20 = cs82_full[(cs82_full.field('MAG_AUTO') < 20) & (19 <= cs82_full.field('MAG_AUTO'))]

cs82_19 = cs82_full[(cs82_full.field('MAG_AUTO') < 19)]


sdss_20_to_21 = sdss_data[(20 <= cs82_full.field('MAG_AUTO'))]

sdss_19_to_20 = sdss_data[(cs82_full.field('MAG_AUTO') < 20) & (19 <= cs82_full.field('MAG_AUTO'))]

sdss_19 = sdss_data[(cs82_full.field('MAG_AUTO') < 19)]


# Size cuts with CS82 FWHM

s_cut = 8

cs82_big_cs82_FWHM = cs82_19[(cs82_19.field('FWHM_IMAGE') > s_cut)]
cs82_small_cs82_FWHM = cs82_19[(cs82_19.field('FWHM_IMAGE') <= s_cut)]

sdss_big_cs82_FWHM = sdss_19[(cs82_19.field('FWHM_IMAGE') > s_cut)]
sdss_small_cs82_FWHM = sdss_19[(cs82_19.field('FWHM_IMAGE') <= s_cut)]



# Size cuts with SDSS PSF FWHM

psf_cut = 0.75

cs82_big_sloan_psf = cs82_full[(sdss_data.field('psffwhm_i').real > psf_cut)]
cs82_small_sloan_psf = cs82_full[(sdss_data.field('psffwhm_i').real <= psf_cut)]

sdss_big_sloan_psf = sdss_data[(sdss_data.field('psffwhm_i').real > psf_cut)]
sdss_small_sloan_psf = sdss_data[(sdss_data.field('psffwhm_i').real <= psf_cut)]


###########################################################################################################

# SECTION PLOTS - CAN I FORMULATE A SUFFICIENTLY GENERAL FUNCTION TO SIMPLIFY THIS PART?

plot_lims = [[[0,5],[0,1],[-90,90]],[[0,5],[0,1],[0,180]]]

plot_titles = ['Comparison 1 - Effective Radius', 'Comparison 2 - Aspect Ratio', 'Comparison 3 - Theta']

units = [[' (arcsec)',' (arcsec)'],['',''],[' (deg)',' (deg)']]


for i in range(len(fieldnames_cs82)):

    fit_x = np.arange(plot_lims[1][i][0],plot_lims[1][i][1],0.01)

    if i == 0:

        if fieldnames_cs82[i] == 'SPHEROID_REFF_WORLD':
            j=0
        else:
            j=1
        

        color_mag_data = [[cs82_20_to_21.field(fieldnames_cs82[i])*3600*(1.68**j),sdss_20_to_21.field(fieldnames_sdss[i]).real],[cs82_19_to_20.field(fieldnames_cs82[i])*3600*(1.68**j),sdss_19_to_20.field(fieldnames_sdss[i]).real],[cs82_19.field(fieldnames_cs82[i])*3600*(1.68**j),sdss_19.field(fieldnames_sdss[i]).real],[fit_x,fit_x]]


        cs82_size_data = [[cs82_big_cs82_FWHM.field(fieldnames_cs82[i])*3600*(1.68**j),sdss_big_cs82_FWHM.field(fieldnames_sdss[i]).real],[cs82_small_cs82_FWHM.field(fieldnames_cs82[i])*3600*(1.68**j),sdss_small_cs82_FWHM.field(fieldnames_sdss[i]).real],[fit_x,fit_x]]

        sdss_size_data = [[cs82_big_sloan_psf.field(fieldnames_cs82[i])*3600*(1.68**j),sdss_big_sloan_psf.field(fieldnames_sdss[i]).real],[cs82_small_sloan_psf.field(fieldnames_cs82[i])*3600*(1.68**j),sdss_small_sloan_psf.field(fieldnames_sdss[i]).real],[fit_x,fit_x]]

        print fieldnames_cs82[i]


    elif i == 2:

        color_mag_data = [[cs82_20_to_21.field(fieldnames_cs82[i]),sdss_20_to_21.field(fieldnames_sdss[i]).real],[cs82_19_to_20.field(fieldnames_cs82[i]),sdss_19_to_20.field(fieldnames_sdss[i]).real],[cs82_19.field(fieldnames_cs82[i]),sdss_19.field(fieldnames_sdss[i]).real],[fit_x-90,fit_x]]

    
        cs82_size_data = [[cs82_big_cs82_FWHM.field(fieldnames_cs82[i]),sdss_big_cs82_FWHM.field(fieldnames_sdss[i]).real],[cs82_small_cs82_FWHM.field(fieldnames_cs82[i]),sdss_small_cs82_FWHM.field(fieldnames_sdss[i]).real],[fit_x-90,fit_x]]


        sdss_size_data = [[cs82_big_sloan_psf.field(fieldnames_cs82[i]),sdss_big_sloan_psf.field(fieldnames_sdss[i]).real],[cs82_small_sloan_psf.field(fieldnames_cs82[i]),sdss_small_sloan_psf.field(fieldnames_sdss[i]).real],[fit_x-90,fit_x]]

        print fieldnames_cs82[i]


    else:

        color_mag_data = [[cs82_20_to_21.field(fieldnames_cs82[i]),sdss_20_to_21.field(fieldnames_sdss[i]).real],[cs82_19_to_20.field(fieldnames_cs82[i]),sdss_19_to_20.field(fieldnames_sdss[i]).real],[cs82_19.field(fieldnames_cs82[i]),sdss_19.field(fieldnames_sdss[i]).real],[fit_x,fit_x]]


        cs82_size_data = [[cs82_big_cs82_FWHM.field(fieldnames_cs82[i]),sdss_big_cs82_FWHM.field(fieldnames_sdss[i]).real],[cs82_small_cs82_FWHM.field(fieldnames_cs82[i]),sdss_small_cs82_FWHM.field(fieldnames_sdss[i]).real],[fit_x,fit_x]]


        sdss_size_data = [[cs82_big_sloan_psf.field(fieldnames_cs82[i]),sdss_big_sloan_psf.field(fieldnames_sdss[i]).real],[cs82_small_sloan_psf.field(fieldnames_cs82[i]),sdss_small_sloan_psf.field(fieldnames_sdss[i]).real],[fit_x,fit_x]]

        print fieldnames_cs82[i]

    
    pltemp.plot_data(color_mag_data, ['.','.','.','-'], ['b','g','r','k'], ['20 < MAG_AUTO < 21', '19 < MAG_AUTO < 20', 'MAG_AUTO < 19',''], 'CS82 - ' + fieldnames_cs82[i] + units[i][0], 'SDSS - ' + fieldnames_sdss[i] + units[i][1], plot_lims[0][i], plot_lims[1][i], plot_titles[i])
    mkdir_p(folder_path + '/plots_sdss/')
    pl.savefig(folder_path + '/plots_sdss/' + fieldnames_sdss[i] + '_mags.png')    
    pl.clf()

    pltemp.plot_data(cs82_size_data, ['.','.','-'], ['b','r','k'], ['CS82 FWHM_IMAGE > '+ str(s_cut*0.186)+'"','CS82 FWHM_IMAGE < '+ str(s_cut*0.186)+'"',''], 'CS82 - ' + fieldnames_cs82[i] + units[i][0], 'SDSS - ' + fieldnames_sdss[i] + units[i][1], plot_lims[0][i], plot_lims[1][i], plot_titles[i])
    pl.savefig(folder_path + '/plots_sdss/' + fieldnames_sdss[i] + '_cs82_fwhm.png')    
    pl.clf()

    pltemp.plot_data(sdss_size_data, ['.','.','-'], ['b','r','k'], ['SDSS psffwhm_i > '+ str(psf_cut)+'"','SDSS psffwhm_i < '+ str(psf_cut)+'"',''], 'CS82 - ' + fieldnames_cs82[i] + units[i][0], 'SDSS - ' + fieldnames_sdss[i] + units[i][1], plot_lims[0][i], plot_lims[1][i], plot_titles[i])
    mkdir_p(folder_path + '/plots_sdss/')
    pl.savefig(folder_path + '/plots_sdss/' + fieldnames_sdss[i] + '_sdss_psf.png')    
    pl.clf()


# Masters' Figure 11 Mag Color:

print j

rel_reff_21 = [cs82_20_to_21.field(fieldnames_cs82[0])*3600*1.68**j,(sdss_20_to_21.field(fieldnames_sdss[0]).real-cs82_20_to_21.field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_20_to_21.field(fieldnames_cs82[0])*3600*1.68**j)]

rel_reff_20 =[cs82_19_to_20.field(fieldnames_cs82[0])*3600*1.68**j,(sdss_19_to_20.field(fieldnames_sdss[0]).real-cs82_19_to_20.field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_19_to_20.field(fieldnames_cs82[0])*3600*1.68**j)]

rel_reff_19 =[cs82_19.field(fieldnames_cs82[0])*3600*1.68**j,(sdss_19.field(fieldnames_sdss[0]).real-cs82_19.field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_19.field(fieldnames_cs82[0])*3600*1.68**j)]

rel_reff = [rel_reff_21, rel_reff_20, rel_reff_19, [np.arange(0, 5, 0.01),0*np.arange(0, 5, 0.01)]]


pltemp.plot_data(rel_reff, ['.','.','.','--'], ['b','g','r','k'], ['20 < MAG_AUTO < 21 / median: ' + str('%.2f' % np.median(rel_reff_21[1])), '19 < MAG_AUTO < 20 / median: '+ str('%.2f' % np.median(rel_reff_20[1])), 'MAG_AUTO < 19 / median: '+ str('%.2f' % np.median(rel_reff_19[1])),''], 'CS82 - ' + fieldnames_cs82[0] + units[0][0], '(reff_sdss - reff_cs82)/reff_cs82', [0, 4], [-1, 5], 'Relative Effective Radius Difference SDSS x CS82')

pl.savefig(folder_path + '/plots_sdss/fig_11_'+fieldnames_sdss[0][0:3]+'.png')    
pl.clf()


# Masters' Figure 11 Size:

print j

rel_reff_big = [cs82_big_cs82_FWHM.field(fieldnames_cs82[0])*3600*1.68**j,(sdss_big_cs82_FWHM.field(fieldnames_sdss[0]).real-cs82_big_cs82_FWHM.field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_big_cs82_FWHM.field(fieldnames_cs82[0])*3600*1.68**j)]

rel_reff_small = [cs82_small_cs82_FWHM.field(fieldnames_cs82[0])*3600*1.68**j,(sdss_small_cs82_FWHM.field(fieldnames_sdss[0]).real-cs82_small_cs82_FWHM.field(fieldnames_cs82[0])*3600*1.68**j)/(cs82_small_cs82_FWHM.field(fieldnames_cs82[0])*3600*1.68**j)]


rel_reff_size = [rel_reff_big, rel_reff_small, [np.arange(0, 5, 0.01),0*np.arange(0, 5, 0.01)]]


pltemp.plot_data(rel_reff_size, ['.','.','--'], ['b','r','k'], ['CS82 FWHM_IMAGE > '+ str(s_cut*0.186)+'" / median: ' + str('%.2f' % np.median(rel_reff_big[1])), 'CS82 FWHM_IMAGE < '+ str(s_cut*0.186)+'" / median: '+ str('%.2f' % np.median(rel_reff_small[1])),''], 'CS82 - ' + fieldnames_cs82[0] + units[0][0], '(reff_sdss - reff_cs82)/reff_cs82', [0, 4], [-1, 5], 'Relative Effective Radius Difference SDSS x CS82')

pl.savefig(folder_path + '/plots_sdss/fig_11_size_'+fieldnames_sdss[0][0:3]+'.png')    
pl.clf()



# Half-light Radius Figures CS82:

print j

hlrad_21 = [cs82_20_to_21.field(fieldnames_cs82[0])*3600*1.68**j,cs82_20_to_21.field('FLUX_RADIUS')*0.186]

hlrad_20 =[cs82_19_to_20.field(fieldnames_cs82[0])*3600*1.68**j,cs82_19_to_20.field('FLUX_RADIUS')*0.186]

hlrad_19 =[cs82_19.field(fieldnames_cs82[0])*3600*1.68**j,cs82_19.field('FLUX_RADIUS')*0.186]

hlrad = [hlrad_21, hlrad_20, hlrad_19, [np.arange(0, 5, 0.01),np.arange(0, 5, 0.01)]]


pltemp.plot_data(hlrad, ['.','.','.','--'], ['b','g','r','k'], ['20 < MAG_AUTO < 21', '19 < MAG_AUTO < 20', 'MAG_AUTO < 19',''], 'CS82 - ' + fieldnames_cs82[0] + units[0][0], 'CS82 - FLUX_RADIUS [0.5]' + units[0][0], [0, 4], [0, 4], 'CS82 Half-Light Radius x Effective Radius')

pl.savefig(folder_path + '/plots_sdss/hlradius_'+fieldnames_sdss[0][0:3]+'.png')    
pl.clf()


# Half-light Radius Figures SDSS:

print j

hlrad_21 = [sdss_20_to_21.field(fieldnames_sdss[0]).real,cs82_20_to_21.field('FLUX_RADIUS')*0.186]

hlrad_20 =[sdss_19_to_20.field(fieldnames_sdss[0]).real,cs82_19_to_20.field('FLUX_RADIUS')*0.186]

hlrad_19 =[sdss_19.field(fieldnames_sdss[0]).real,cs82_19.field('FLUX_RADIUS')*0.186]

hlrad = [hlrad_21, hlrad_20, hlrad_19, [np.arange(0, 5, 0.01),np.arange(0, 5, 0.01)]]


pltemp.plot_data(hlrad, ['.','.','.','--'], ['b','g','r','k'], ['20 < MAG_AUTO < 21', '19 < MAG_AUTO < 20', 'MAG_AUTO < 19',''], 'CS82 - ' + fieldnames_cs82[0] + units[0][0], 'CS82 - FLUX_RADIUS [0.5]' + units[0][0], [0, 4], [0, 4], 'CS82 Half-Light Radius x SDSS Effective Radius')

pl.savefig(folder_path + '/plots_sdss/hlradius_sdss_'+fieldnames_sdss[0][0:3]+'.png')    
pl.clf()

# The following piece of code should be placed elsewhere, decide later where.

'''
# =============================== Creating the masking file for the full cs82 mag 21 S/N 5 CLASS_STAR < 0.95


# Transform from radec to xy:

image_name = '/home/brunomor/Documents/Trabalho/S82p28m_y.V2.7A.swarp.cut.fits'

hdu_image = pyfits.open(image_name)

image_header = hdu_image[0].header

cs82_ra_list = cs82_ra.tolist()
cs82_dec_list = cs82_dec.tolist()

pixcrd = []
xlist = []
ylist = []

for i in range(len(cs82_ra)):
    pixcrd.append(wcscnv.radec2xy(image_header,cs82_ra_list[i],cs82_dec_list[i]))
    xlist.append(pixcrd[i][0])
    ylist.append(pixcrd[i][1])

ascd.write_ds9cat(xlist,ylist,10,'circle','red', outputfile = 'cs82_full_cuts_objects.reg')


# =============================== Creating the Mask File for the full boss/sdss data


# Transform from radec to xy:

image_name = '/home/brunomor/Documents/Trabalho/S82p28m_y.V2.7A.swarp.cut.fits'

hdu_image = pyfits.open(image_name)

image_header = hdu_image[0].header

sdss_ra_list = sdss_ra.tolist()
sdss_dec_list = sdss_dec.tolist()

pixcrd = []
xlist = []
ylist = []

for i in range(len(sdss_ra)):
    pixcrd.append(wcscnv.radec2xy(image_header,sdss_ra_list[i],sdss_dec_list[i]))
    xlist.append(pixcrd[i][0])
    ylist.append(pixcrd[i][1])

ascd.write_ds9cat(xlist,ylist,30,'circle','yellow', outputfile = 'boss_full_objects.reg')
'''

print len(sdss_data)
