

##@package get_header_parameter [formerly GetHeaderParameter] 
# Reads parameter value from image header.

# 
#@param image_file, param : image FITS filename and parameter, (string,string) 
#@return parvalue : header parameter value (string)

def get_header_parameter( image_file, param ):
    import pyfits; 

    # Select image parameters in header:
    _header = pyfits.getheader(image_file);

    parvalue = _header[param]; 

    return (parvalue);


# ----- 


