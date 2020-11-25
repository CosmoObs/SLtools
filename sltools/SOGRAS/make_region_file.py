import astropy.io.fits as pyfits
import glob
import os

#from sltools.coordinate.wcs_conversion import radec2xy_GBC
from sltools.coordinate.wcs_conversion import xy2radec


def s82_operations(catalog_path, image_directory, image_path):
    hdulist = pyfits.open(catalog_path)

    data = hdulist[1].data

    image_hdulist = pyfits.open(image_path)
    image_header = image_hdulist[0].header

    n_axis1 = image_header['NAXIS1']
    n_axis2 = image_header['NAXIS2']

    coord_radec_min = xy2radec(image_path, 0, 0)
    coord_radec_max = xy2radec(image_path, n_axis1, n_axis2)

    mask1 = data.field('RA') > min( [coord_radec_min[0], coord_radec_max[0]] )

    mask2 = data.field('RA') < max( [coord_radec_min[0], coord_radec_max[0]] )

    mask3 = data.field('DEC')> min( [coord_radec_min[1], coord_radec_max[1]] )
 
    mask4 = data.field('DEC')< max( [coord_radec_min[1], coord_radec_max[1]] )

    data_filtred = data[mask1 & mask2 & mask2 & mask3 & mask4]

    if len(data_filtred) != 0:
        reg_out_name = image_path[len(image_directory):] + '_dr8_phot.reg'
        print reg_out_name
        reg_out = open(reg_out_name,'w')

        print >> reg_out, '# Region file format: DS9 version 4.1'
        print >> reg_out, '# Filename:', image_path[len(image_directory):]
        print >> reg_out, '#file made by G. B. Caminha (gbcaminha@gmail.com) to identify the SDSS-DR8 positions in images'
        print >> reg_out, 'global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1'
        print >> reg_out, 'fk5'

        for i in range(len(data_filtred)):
            out_str = "circle(" + str(data_filtred.field('RA')[i]) + "," + str(data_filtred.field('DEC')[i]) + "," + str(1.5) + "\")"
            print >> reg_out, out_str


#Main Program



IMAGE_DIRECTORY = '/home/gbcaminha/SOGRAS/2008b_data/2008b_reduced/'
CATALOG_DIRECTORY = '/home/gbcaminha/SDSS/catalogs/SOGRAS_fields/'

IMAGE_FOLDERS = os.listdir(IMAGE_DIRECTORY)
IMAGE_FOLDERS.sort()

CALATOG_FILES = glob.glob( os.path.join(CATALOG_DIRECTORY, '*.fit'))
CALATOG_FILES.sort()

IMAGE_FILES = []

for i in range(len(IMAGE_FOLDERS)):
    IMAGE_FILES.append( glob.glob( os.path.join(IMAGE_DIRECTORY + IMAGE_FOLDERS[i], '*r.fits'))[0] )
    IMAGE_FOLDERS[i] = IMAGE_DIRECTORY + IMAGE_FOLDERS[i] + '/'





for i in range(len(IMAGE_FOLDERS)):
    s82_operations(CALATOG_FILES[i], IMAGE_FOLDERS[i], IMAGE_FILES[i])









