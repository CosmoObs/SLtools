def get_catalog_data( hdulist, *args ): 
   

    if ( len(args) == 0 ): 
        return (None); 

    tbdata = hdulist[1].data; 

    dic = {}; 
    for arg in args: 
        dic[arg] = tbdata.field(arg); 

    return (dic); 

# ----- 


