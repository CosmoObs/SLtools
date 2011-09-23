#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import sys
import os
import glob
import numpy as np
import pylab as pl 
import pyfits as pf
from sltools.image import sextractor

pathfile=sys.argv[1]

def get_paths(pathfile):
    path_dict=sextractor.read_config(pathfile)
    return path_dict

def get_configs(path_dict):
    image_name= path_dict['image_path']+'S82p28m_2000.fits'

    header=pf.getheader(image_name)

    first_sex_config = sextractor.read_config(path_dict['sex_path']+'default_CS82_Maria.sex')

    header_info = {"GAIN":header['GAIN'], "SEEING_FWHM":header['SEEING'], "MAG_ZEROPOINT":header['MAGZP']}

    first_sex_config.update(header_info)

    cs82_sex_config = {"-c":path_dict['sex_path']+'default_CS82_Maria.sex', "DETECT_MINAREA": '5', "DETECT_THRESH": '1.5', "ANALYSIS_THRESH": '1.5', "DEBLEND_NTHRESH":'64' , "DEBLEND_MINCONT": '0.0005', "CATALOG_NAME": path_dict['out_catalog_path']+'cs82_2000.fit',"CATALOG_TYPE": 'FITS_LDAC', "CHECKIMAGE_NAME": path_dict['out_segs_path']+'cs82_2000_seg.fits',"CHECKIMAGE_TYPE": 'SEGMENTATION', "STARNNW_NAME": path_dict['nnw_path']+'default.nnw', "FILTER_NAME": path_dict['conv_path']+'default.conv',"WEIGHT_IMAGE":''} #, "WEIGHT_IMAGE":'',  , OBJECTS  path_dict['out_obj_path']+'cs82_2000_obj.fits'


    first_sex_config.update(cs82_sex_config)
    return first_sex_config  


def get_params(path_dict):
    filename=path_dict['param_path']+'default_CS82_Maria_SE_v2.13.2.param'
    param_list= sextractor.read_config(filename).keys()
    return param_list

path_dict = get_paths(pathfile)
image = path_dict['image_path']+'S82p28m_2000.fits'
first_config_dict = get_configs(path_dict)
first_param_list=get_params(path_dict)


sextractor.run(image, params=first_param_list, args=first_config_dict, preset='', temp_dir='', quiet=False)




