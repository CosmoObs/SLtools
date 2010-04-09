
##@package add_arcs_2_image
#
# Given the Halo central position (X,Y) inside
# an image (DC tile), add arc poststamp to image.
#
#@param image_array, arc_image_array, X_halo, Y_halo
#@return image_array, image_header

# ---
def add_arcs_2_image( image_array, arc_image_array, X_halo, Y_halo ):
    import string;
    import logging

    logging.debug(' X_halo,Y_halo = %d, %d' % (X_halo,Y_halo) )


    # Definicao da janela de pixels a ser utilizada na soma do arco `a tile
    #
    DY2 = int(len(arc_image_array)/2);
    if ( len(arc_image_array)%2 ) :
        minDY = Y_halo-DY2-1;
    else:
        minDY = Y_halo-DY2;
    maxDY = Y_halo+DY2;

    DX2 = int(len(arc_image_array[0])/2);
    if ( len(arc_image_array[0])%2 ) :
        minDX = X_halo-DX2-1;
    else:
        minDX = X_halo-DX2;
    maxDX = X_halo+DX2;


    # Finalmente, a adicao das imagens:
    #
    image_array[minDY:maxDY,minDX:maxDX] = image_array[minDY:maxDY,minDX:maxDX] + arc_image_array;


    return (image_array);

# -
