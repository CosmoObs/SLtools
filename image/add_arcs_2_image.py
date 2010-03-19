
##@package add_arcs_2_image [formerly Add_Arcs2Image]

# Addition of Arc images to sky image on respective Halos.

# 
#@param hdulist,image_file,list of tuples with haloIDs and arc filenames

#@return image_array, image_header

def add_arcs_2_image( image_array, arc_file, Xhalo, Yhalo ):
    import string;
    import pyfits;
    import logging

    logging.debug(' Xhalo,Yhalo = %d, %d' % (Xhalo,Yhalo) )

    # Read Fits file to array
    arcimg_array, arcimg_header = pyfits.getdata( arc_file, header=True );

    # Definicao da janela de pixels a ser utilizada na soma do arco `a tile
    DY2 = int(len(arcimg_array)/2);
    if ( len(arcimg_array)%2 ) :
        minDY = Yhalo-DY2-1;
    else:
        minDY = Yhalo-DY2;
    maxDY = Yhalo+DY2;

    DX2 = int(len(arcimg_array[0])/2);
    if ( len(arcimg_array[0])%2 ) :
        minDX = Xhalo-DX2-1;
    else:
        minDX = Xhalo-DX2;
    maxDX = Xhalo+DX2;


    # Finalmente, a adicao das imagens:
    image_array[minDY:maxDY,minDX:maxDX] = image_array[minDY:maxDY,minDX:maxDX] + arcimg_array;

    return (image_array); #,image_header);

