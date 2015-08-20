#!/usr/bin/env python
# ==================================
# Authors:
# Gabriel Caminha - gbcaminha@gmail.com
# ==================================

"""Package to make cuts in images given a coordinate catalogue"""

##@package catalog2cutout
#
#
# This package make cutouts in a fits image given a fits calaog with a list of 
# RA and DEC
#

import pyfits
import sys
import string

from sltools.image.header_funcs import get_ra_dec_limits
from sltools.image.imcp import cutout

from sltools.coordinate.wcs_conversion import radec2xy, xy2radec

from optparse import OptionParser

def catalog2cutout_pipe(im_path, catalog_path, dim):
    """
    make cuts in images given a coordinate catalogue


    Input:
     - im_path : str
        Path for a fits image

     - catalog_path : str
        Path for a fits catalog

     - dim : str
        Output stamps shape

    
    ---
    """
    im_lim = get_ra_dec_limits(im_path)
    im_lim_xy = ( radec2xy(im_path, im_lim[0][0], im_lim[1][0]),
                  radec2xy(im_path, im_lim[0][1], im_lim[1][1]) )

    im_x_min = min(im_lim_xy[0][0], im_lim_xy[1][0])
    im_x_max = max(im_lim_xy[0][0], im_lim_xy[1][0])

    im_y_min = min(im_lim_xy[0][1], im_lim_xy[1][1])
    im_y_max = max(im_lim_xy[0][1], im_lim_xy[1][1])

    hdulist_cat = pyfits.open(catalog_path)

    hdulist_im = pyfits.open(im_path, memmap=True)
    im_header = hdulist_im[0].header
    im_data = hdulist_im[0].data

    #print hdulist_im.filename()



    cat_coord = zip( hdulist_cat[1].data.field("RA"),
                     hdulist_cat[1].data.field("DEC") )

    cat_coord_xy = []
    coord_ok = []
    for i in range(len(cat_coord)):
        cat_coord_xy.append(radec2xy( im_path, \
                                      cat_coord[i][0], cat_coord[i][1]) )

        test_lim_x = cat_coord_xy[i][0] > im_x_min \
                     and cat_coord_xy[i][0] < im_x_max
        test_lim_y = cat_coord_xy[i][1] > im_y_min \
                     and cat_coord_xy[i][1] < im_y_max
        if  test_lim_x and test_lim_y :
            coord_ok.append(cat_coord_xy[i])


    dx, dy = string.split(dim,sep=',')
    imgcut_out, hdr_out = cutout( im_data, im_header, coord_unit='pixel', \
                                    xo=coord_ok[0][0], yo=coord_ok[0][1], \
                        size_unit='pixel', x_size=dx, y_size=dy, mask=None)

    for i in range(len(coord_ok)):
        outfits = "out_catalog2cutout_"+str(i)+".fits"
        pyfits.writeto( outfits, imgcut_out, hdr_out )

    return 0




if __name__ == "__main__" :

    USAGE = "\n  %prog como se a chamada do programa FIXME  [options]"
    PARSER = OptionParser(usage=USAGE)

    PARSER.add_option('--imPATH', dest='im_path', default='',
        help = "FIXME")
    PARSER.add_option('--catPATH', dest='cat_path', default='',
        help = "FIXME")
    PARSER.add_option('-s',
                dest='dim', default='250,250',
                help="Output stamps shape default is 250,250")

    (OPTS, ARGS) = PARSER.parse_args()

    if not( OPTS.im_path or OPTS.cat_path ):
        PARSER.print_help()
        sys.exit(0)


    catalog2cutout_pipe(OPTS.im_path, OPTS.cat_path, OPTS.dim)


