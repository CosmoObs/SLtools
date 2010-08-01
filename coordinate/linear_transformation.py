""" linear_transformation """

##@package linear_transformation

import sys;

"""Linear transformation of coordinates"""

##@package transform

# =================================================
def linear_trasnformation( ra, dec, RAshift=0, DECfactor=1 ):
    """Shift and reorder of RA/DEC coordinates

    * Since mock catalogues and DC images objects don't agree
    in RA/DEC position, coordinates values have to be re-centered.

    shift_radec( ra.float, dec.float [,...])

    Function shift given RA and DEC values as follow:
    RA_out = RA_in + RAshift
    DEC_out = DEC_in * DECfactor

    RA and DEC values are constrained to 0 < RA < 360 and -90 < DEC < 90

    Input:
     - ra        : scalar or tuple (RA_min,RA_max) of Right Ascension
     - dec       : scalar or tuple (DEC_min,DEC_max) of Declination
     - RAshift   : value to translate given RA value(s)
     - DECfactor : value to multiply given DEC value(s)

    Output:
     - (ra,dec) : scalars or tuples of coordinates

    """
	
    try:
        if len(ra) > 1 and len(dec) > 1 :
            ra_min = float(ra[0]) + float(RAshift);
            ra_max = float(ra[1]) + float(RAshift);
            dec_min = float(dec[0]) * float(DECfactor);
            dec_max = float(dec[1]) * float(DECfactor);

            if ( dec_max < dec_min ) :
                _temp = dec_min;
                dec_min = dec_max;
                dec_max = _temp;
            if ( ra_max < ra_min ) :
                _temp = ra_min;
                ra_min = ra_max;
                ra_max = _temp;

            if ra_min > 360 and ra_max > 360 :
                ra_min -= 360;
                ra_max -= 360;
            if ra_min < 0 and ra_max < 0 :
                ra_min += 360;
                ra_max += 360;

            ra = (ra_min,ra_max);
            dec = (dec_min,dec_max);

    except TypeError:
        ra = float(ra) + float(RAshift);
	if ra > 360 :
	    ra -= 360;
	elif ra < 0 :
	    ra += 360;
	else:
	    pass;

        dec = float(dec) * float(DECfactor);


    return ra,dec;

# ---
