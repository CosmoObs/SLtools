

##@package get_header_parameter [formerly GetHeaderParameter] 
# Reads parameter value from image header.

# 
#@param image_file, param : image FITS filename and parameter, (string,string) 
#@return parvalue : header parameter value (string)

def get_header_parameter( image_file, *parargs ):
    import pyfits; 

    # Select image parameters in header:
    _header = pyfits.getheader(image_file);

    param_list = [];
    for _param in parargs :
        param_list.append( _header[_param] ); 

    return (param_list);


# ----- 


