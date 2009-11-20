

##@package get_halo_parameters [formerly GetHaloParameters]
# Get from Catalog M200, z, RA, DEC for given HaloID
#
#@param HDUlist : pyfits catalog
#@param HaloID : a Halo Identification number
#@return [Mass200, zL, RA,DEC] : a list with four parameters
def get_halo_parameters( hdulist, haloids, RAshift=0, DECfactor=1 ):
    import numpy;

    # OK. So, now we are able to read some parameters...
    tbdata = hdulist[1].data          # First file extension

    M200 = [];
    zL = [];
    RA = [];
    DEC = [];
    for _id in haloids:
        # Read given Halo_ID index in table
        ind = numpy.where( tbdata.field('HALOID') == _id )[0][0]
        halo = tbdata[ind]

         # Nice...now we have our Halo selected. Lets get what we came for... ;-)
        M200.append( halo.field('M200') );                     # units: Msun / h
        zL.append( halo.field('Z') );                                # redshift
        RA.append( halo.field('RA') - RAshift );                # Necessary conversion for stanford(Mock) catalog.
        DEC.append( halo.field('DEC') * DECfactor );       # .. same for DEC. Both in degrees :p

    # Return asked Halo parameters as a 4 length list
    return (M200,zL,RA,DEC)

# -----


