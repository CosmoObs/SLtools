
##@package get_header_parameter
#
# Read parameter value from image header.
# 
#@param image_file, param : image FITS filename and parameter, (string,string) 
#@return parvalue : header parameter value (string)

# ---
import pyfits;
# ---
def get_header_parameter( image_file, *parargs ):

    _header = pyfits.getheader( image_file );

    # Get image parameters in header..
    #
    param_list = [];
    for _param in parargs :
        param_list.append( _header[_param] ); 

    return (param_list);
# -
