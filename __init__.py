"""
Strong Lensing Tools

SLtools is developed with the main goal of supporting astronomical daily-common 
data processing tools as well as products from our group research.

Documentation is provided by python docstrings *and* PDF (also HTML) pages 
located in library's "doc/" directory. File pages are generated with Doxygen
() plus doxypy() filter.

SLtools makes massive use of Numpy() and PyFits() packages. Take the chance to
read their documentation also if you didn't do that yet.

Since the library has multi-purposed packages, to have access to each of them 
one has to do explicit import on it.

The following examples show how one can import a package (e.g, image):

    import sltools.image
    
    from sltools import image
    
    import sltools.image as slimage

Use the help command to see the available modules or details for this package:

    help(sltools.image)

    help(image)
    
    help(slimage)

The following command makes (image) modules available in the current namespace:

    from sltools.image import *


SLtools' contents are breafly presented below:

    - arcellipse :
    
    - catalog : FITS and ASCII data tables handlers
        - halos :
        - sources :
    
    - coordinate : Sky/Image position/data convertion routines
    
    - geometry : general mathematical/geometrical data processing
    
    - gravlens :
    
    - image : FITS image processing routines
    
    - io : Input/Output handlers used within packages
    
    - lens :
    
    - pipelines : larger pieces of codes using more packages from library
        - finders : arc finder procedures

    - plot :
    
    - string : strings handling routines
    
---
"""
