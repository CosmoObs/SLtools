
##@package add_arcs_2_image [formerly Add_Arcs2Image]

# Addition of Arc images to sky image on respective Halos.

# 
#@param hdulist,image_file,list of tuples with haloIDs and arc filenames

#@return image_array, image_header

def add_arcs_2_image( hdulist, image_file, halos_and_arcs ):
    
    image_array, image_header = pyfits.getdata( image_file, header=True );


    lista = [];
    tupla = ('Halo','Arc_conv_file','RA(deg)','DEC(deg)','X(px)','Y(px)');
    lista.append(tupla);


    for pair in HalosandArcs_list:

        # Read Halo ID and corresponding arc(s) image file
        HaloID = pair[0];
        Arc_filename = pair[1];

        HaloParameters = GetHaloParameters( hdulist,HaloID );
        ra = HaloParameters[2];
        dec = HaloParameters[3];
        _out = commands.getoutput('%s %s %s %s' % (sky2xy,Image_filename,ra,dec) );
        _out = string.split(_out);

        Xhalo = int(string.atof(_out[4]));
        Yhalo = int(string.atof(_out[5]));

        if ( string.find(Arc_filename,'.fits') == -1 ) :
            continue;

        # Read Fits file to array
        arcimg_array, arcimg_header = pyfits.getdata( arc_file, header=True );

        # Definicao da janela de pixels a ser utilizada na soma do arco `a tile
        DY2 = int(len(ArcImg_array)/2);
        if ( len(ArcImg_array)%2 ) :
            minDY = Yhalo-DY2-1;
        else:
            minDY = Yhalo-DY2;
        maxDY = Yhalo+DY2;

        DX2 = int(len(ArcImg_array[0])/2);
        if ( len(ArcImg_array[0])%2 ) :
            minDX = Xhalo-DX2-1;
        else:
            minDX = Xhalo-DX2;
        maxDX = Xhalo+DX2;


        # Finalmente, a adicao das imagens:
        image_array[minDY:maxDY,minDX:maxDX] = image_array[minDY:maxDY,minDX:maxDX] + arcimg_array;


    return (image_array,image_header);

