""" Module to deal with FITS (image) header entries """

import sys;
import math as m;

def read_pixelscale( hdr ):
    """ Read out header dimpix (pixel_scale) from CD_* keys """

    try:
        pixel_scale = hdr['PIXSCALE'];
        return pixel_scale;
    except:
        pass;
        
    try:
        CD1_1 = float(hdr['CD1_1']);
        CD1_2 = float(hdr['CD1_2']);
        CD2_1 = float(hdr['CD2_1']);
        CD2_2 = float(hdr['CD2_2']);
    except:
        CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1;

    try:
        CUNIT1 = str(hdr['CUNIT1']);
        CUNIT2 = str(hdr['CUNIT2']);
    except:
        print >> sys.stderr, "Warning: Unable to find 'CUNIT[1|2]' parameters in header instance for image coordinates unit.";
        print >> sys.stderr, "         Degrees ('deg') is being used as default unit value.";
        CUNIT1 = 'deg';
        CUNIT2 = 'deg';

    conv_factor = 1.;

    # Convertion from degrees to arcsec units..
    #
    if ( CUNIT1 == 'deg' and CUNIT2 == 'deg' ):
        conv_factor = 3600.0;

    # Scale is given, using CD1_1, CD1_2, CD2_1, CD2_2 header keys:
    #
    pixel_scale = m.sqrt((CD1_1**2 + CD1_2**2 + CD2_1**2 + CD2_2**2)/2.) * conv_factor;

    return (pixel_scale);

# ---

def update_coordinates(hdr, x_ini, y_ini):
    """Update header information regarding world coordinate system"""


    LTV1 = None;
    LTV2 = None;

    NAXIS1 = int(hdr['NAXIS1']);
    NAXIS2 = int(hdr['NAXIS2']);

    try:
        CRPIX1 = float(hdr['CRPIX1']);
        CRPIX2 = float(hdr['CRPIX2']);
    except:
        CRPIX1 = CRPIX2 = 1.0;

    try:
        CD1_1 = float(hdr['CD1_1']);
        CD1_2 = float(hdr['CD1_2']);
        CD2_1 = float(hdr['CD2_1']);
        CD2_2 = float(hdr['CD2_2']);
    except:
        CD1_1 = CD1_2 = CD2_1 = CD2_2 = 1.0;
        
    # Check the edges and update header keyworks for new image..
    #
    #	if ( x_ini >= 0  and  y_ini >= 0 ):
    #		LTV1 = -1*x_ini;
    #		LTV2 = -1*y_ini;
    #		CRPIX1 = CRPIX1 + LTV1;
    #		CRPIX2 = CRPIX2 + LTV2;
    #	if ( x_ini < 0  and  y_ini >= 0 ):
    #		LTV2 = -1*y_ini;
    #		CRPIX2 = CRPIX2 + LTV2;
    #	if ( y_ini < 0  and  x_ini >= 0 ):
    #		LTV1 = -1*x_ini;
    #		CRPIX1 = CRPIX1 + LTV1;
    LTV1 = -1*x_ini;
    LTV2 = -1*y_ini;
    CRPIX1 = CRPIX1 + LTV1;
    CRPIX2 = CRPIX2 + LTV2;


    # Add some keys to header..
    #
    WCSDIM = 2
    CDELT1  =  CD1_1;
    CDELT2  =  CD2_2;
    LTM1_1 = 1.0
    LTM2_2 = 1.0
    WAT0_001= 'system=image'
    WAT1_001= 'wtype=tan axtype=ra'
    WAT2_001= 'wtype=tan axtype=dec'


    # Header update..
    #
    if (LTV1 != None) :
        hdr.update('LTV1',LTV1)
        hdr.update('CRPIX1',CRPIX1)
    if (LTV2 != None) :
        hdr.update('LTV2',LTV2)
        hdr.update('CRPIX2',CRPIX2)
    hdr.update('WCSDIM',WCSDIM)
    hdr.update('CDELT1',CDELT1)
    hdr.update('CDELT2',CDELT2)
    hdr.update('LTM1_1',LTM1_1)
    hdr.update('LTM2_2',LTM2_2)
    hdr.update('WAT0_001',WAT0_001)
    hdr.update('WAT1_001',WAT1_001)
    hdr.update('WAT2_001',WAT2_001)


    return (hdr);
