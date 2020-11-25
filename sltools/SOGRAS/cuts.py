import os
import astropy.io.fits as pyfits
import glob #useful to fory all files in a directory

from sltools.coordinate.wcs_conversion import radec2xy

def find_obj(path, coordinates_RA, coordinates_DEC, obj_id):

   file_list = [] #list to recive the fits files names that are in path
   for file2 in glob.glob( os.path.join(path, '*.fits') ):
     file_list.append(file2)
   del file2
   file_list.sort()


   for i in range(0,len(file_list)):
      print "Working in file:",file_list[i][len(path)+1:]

      im_hdulist = pyfits.open(file_list[i])
      im_header = im_hdulist[0].header

      n_axis1 = im_header['NAXIS1']
      n_axis2 = im_header['NAXIS2']

      #print n_axis1, n_axis1

      for j in range(0,len(coordinates_RA)):
         coord_XY = radec2xy(file_list[i],coordinates_RA[j],coordinates_DEC[j]) 
         if coord_XY[0]>0 and coord_XY[0]<n_axis1 and coord_XY[1]>0 and coord_XY[1]<n_axis2:
            print j, coord_XY
            comand = 'imcp.py ' + str(file_list[i]) + ' ' + str(coord_XY[0]) + ' ' + str(coord_XY[1]) + ' -s 900,900 -o ' +  '_'  + str(obj_id[n]) + "_" + str(file_list[i][len(path)+1:]) + '_' #comando para fazer o corte nas imagens



#print radec2xy(infile,0,0) 








#infile = '/home/gbcaminha/CS82/images/S82p46p_y.V2.3A.swarp.cut.fits'
path = '/home/gbcaminha/CS82/images'

geach_file = open('geach_deg2.txt', 'r')
coordinates_RA = []
coordinates_DEC = []
obj_id = []


for x in range(0,60):
  line = geach_file.readline()
  split = line.split()

  obj_id.append(int(split[0]))
  coordinates_RA.append(float(split[1]))
  coordinates_DEC.append(float(split[2]))


  #print x, obj_id[x], coordinates_RA[x], coordinates_DEC[x]

find_obj(path, coordinates_RA, coordinates_DEC, obj_id)








