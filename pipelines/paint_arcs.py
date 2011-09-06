#!/usr/bin/python
# =====================================================
# Authors:
# Cristina Furlanetto - furlanetto.cristina@gmail.com
# =====================================================

def run(params,seeing,background,suffix_name=''):

    hdu = paint_arcs(params)
    
    arc_array = hdu.data
    arc_header = hdu.header
    
    pyfits.writeto("arc"+suffix_name+".fits",arc_array,arc_header,clobber=True)
    
    dim_x = arc_header["NAXIS1"]
    dim_y = arc_header["NAXIS2"]
    
    # Convolving the arc
    seeing = convert_arcsec_pix(seeing,pix_scale)	
    sigma = seeing / (2.*math.sqrt(2.*math.log(2.)))
    arc_header.update('APSF',seeing,comment='Arc psf (arcsec)')	
    
    arc_conv_array = convolve.gauss_convolution_fft(arc_array,3.,sigma)
    pyfits.writeto("arc_conv"+suffix_name+".fits",arc_conv_array,arc_header,clobber=True)
        
    # Adding Poisson noise
    arc_conv_noise_array = add_noise.add_poisson_noise(arc_conv_array)
    pyfits.writeto("arc_conv_noise"+suffix_name+".fits",arc_conv_noise_array,arc_header,clobber=True)	
    
    # Creating background image with noise
    bkg_array = background * numpy.ones((dim_y,dim_x)) 
    
    bkg_noise_array = add_noise.add_poisson_noise(bkg_array)
    
    pyfits.writeto("background"+suffix_name+".fits",bkg_noise_array,clobber=True)
    
    # Adding arc to the background image
    arc_conv_noise_bkg = arc_conv_noise_array + bkg_noise_array
    arc_header.update('BKG',background,comment='Background mean value')	
    
    pyfits.writeto("arc_conv_noise_bkg"+suffix_name+".fits",arc_conv_noise_bkg,arc_header,clobber=True)	

    

#params = [7.,1.,0.6,30.,230.,20.,100.,1000.,31.04]
#params = [7.,1.,0.6,30.,230.,20.,100.,1000.,31.04]
#params = [10.,4.,0.6,100.,300.,20.,100.,1000.,31.04]






#run(params,0.6,500.)

