

##@package get_catalog_data 
# Get lists of parameters from given pyfits catalog 
# 
#@param HDUlist: pyfits catalog 
#@param comma separated list of strings with parameters name 
#@return dictionary with the lists of parameters 
def get_catalog_data( hdulist, *args ): 
   

    if ( len(args) == 0 ): 
        return (None); 

    tbdata = hdulist[1].data; 

    dic = {}; 
    for arg in args: 
        dic[arg] = tbdata.field(arg); 

    return (dic); 

# ----- 


