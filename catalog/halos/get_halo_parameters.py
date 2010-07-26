
def get_halo_parameters( hdulist, haloid, RAshift=0, DECfactor=1 ):
    import numpy;

    # OK. So, now we are able to read some parameters...
    tbdata = hdulist[1].data;          # First file extension

    ind = numpy.where( tbdata.field('HALOID') == haloid )[0][0];
    halo = tbdata[ind];
    
    M200 = float(halo.field('M200'));
    zL = float(halo.field('Z'));
    RA = float(halo.field('RA')) + RAshift;
    DEC = float(halo.field('DEC')) * DECfactor;

    # Return asked Halo parameters as a 4 length list
    return (M200,zL,RA,DEC)


