
import astropy.io.fits as pyfits

import sltools.pipelines.table_matching as p
from sltools.coordinate.wcs_conversion import radec2xy_GBC
import numpy as np

#                   
#                

cluster_id = '8103'

#sex_catalog_path = '/home/gbcaminha/SOGRAS/2008b_data/2008b_reduced/10755/10755_final_phot.cat'
sex_catalog_path = '/home/gbcaminha/SOGRAS/2008b_data/2008b_reduced/' + cluster_id + '/' + cluster_id + '_final_phot.cat'

obj_id, x_sex, y_sex, ra_sex, dec_sex = np.loadtxt(sex_catalog_path, usecols=(0,1,2,3,4), unpack=True)

xy_sex = tuple(zip(x_sex, y_sex))

radec_sex = tuple(zip(ra_sex, dec_sex))

#sdss_catalog_path = '/home/gbcaminha/SDSS/catalogs/SOGRAS_fields/SOGRAS_' + cluster_id + '_gbcaminha.fit'
sdss_catalog_path = '/home/gbcaminha/SDSS/galSpecInfo-dr8-S82.fits'
hdulist = pyfits.open(sdss_catalog_path)

#sogras_image_path = '/home/gbcaminha/SOGRAS/2008b_data/2008b_reduced/10755/10755i.fits'
#image_hdulist = pyfits.open(sogras_image_path)

#print hdulist[1].columns

data = hdulist[1].data

ra_sdss = data.field('RA')
dec_sdss = data.field('DEC')
#phot_z = data.field('z')
#phot_z_err = data.field('zErr')

#for the DR8 spectro
phot_z = data.field('Z')
phot_z_err = data.field('Z_ERR')
#print len(ra_sdss)
#xy_sdss = []

#for i in range(len(ra_sdss)):
#    xy_sdss.append(radec2xy_GBC(image_hdulist,ra_sdss[i],dec_sdss[i]))


#xy_sdss = tuple(xy_sdss)
radec_sdss = tuple(zip(ra_sdss, dec_sdss))

#match = p.nearest_neighbour(xy_sdss,xy_sex)
match = p.nearest_neighbour(radec_sdss,radec_sex)

max_dist = 1.5*0.000277778 #1'' in degrees


phot_z_match = []
phot_z_err_match = []

good_match = []

for i in range(len(match)):
    if match[i][1] < max_dist:
        good_match.append((match[i][0],phot_z[i],phot_z_err[i]))
        #print xy_sdss[i], xy_sex[match[i][0]], match[i]

#print good_match
good_match.sort()

#for i in range(len(good_match)): print good_match[i]
#print len(good_match)
#print len(radec_sex)

count = 0

for i in range(len(radec_sex)):
    if len(good_match)>0 and i == good_match[count][0] and count < len(good_match)-1:
        phot_z_match.append(good_match[count][1])
        phot_z_err_match.append(good_match[count][2])
        count = count+1
        #print i, count
    else:
        phot_z_match.append(99.000)
        phot_z_err_match.append(99.000)

#for i in range(len(phot_z_match)):
#    print i+1, phot_z_match[i], phot_z_err_match[i]



file_sex = open(sex_catalog_path, 'r')

count = 0
file_all = []
line_photz = ''

for line in file_sex:
    if line.split()[0] == '#':
        file_all.append(line[:len(line)-1])
    else:
        line_photz = line[:len(line)-1] + ' ' + str(phot_z_match[count]) + '   ' + str(phot_z_err_match[count])
        file_all.append(line_photz)
        count = count + 1

file_out_name = cluster_id + '_final_spec_valued.cat'
file_out = open(file_out_name, 'w')
for i in range(len(file_all)):
    print >> file_out, file_all[i]
    if file_all[i].split()[0] == '#' and file_all[i+1].split()[0] != '#':
        print >> file_out,'#  20 dr8_phot_z       spectrocopic redshift from DR8                  [redshift]'
        print >> file_out,'#  21 dr8_phot_z_err   spectrocopic redshift error from DR8            [redshift]'





