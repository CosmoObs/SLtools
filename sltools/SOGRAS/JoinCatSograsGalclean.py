#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import sys
import os
import glob

# path where are the cluster's files (galclean files)
my_path='/home/maria/Area_de_Trabalho/SOAR/morfologia_sogras/'  # for morphology
# my_path='/home/maria/Area_de_Trabalho/SOAR/Photoz_Sogras/'   # for photoz

files=os.listdir(my_path)

# separating between low and high z
low_z=['16177i_sub_back_final_morphology.txt', '15173i_sub_back_final_morphology.txt', '13733i_sub_back_final_morphology.txt', '6731i_sub_back_final_morphology.txt']
high_z=['6214i_sub_back_final_morphology.txt', '15106i_sub_back_final_morphology.txt',  '10755i_sub_back_final_morphology.txt', '13733i_sub_back_final_morphology.txt', '10920i_sub_back_final_morphology.txt', '28721i_sub_back_final_morphology.txt', '13945i_sub_back_final_morphology.txt', '12573i_sub_back_final_morphology.txt', '15260i_sub_back_final_morphology.txt', '5330i_sub_back_final_morphology.txt', '8103i_sub_back_final_morphology.txt', '13380i_sub_back_final_morphology.txt', '15127i_sub_back_final_morphology.txt', '16628i_sub_back_final_morphology.txt']

'''
# for files with photo-z:
low_z=['16177_final_phot_valued.cat', '15173_final_phot_valued.cat', '13733_final_phot_valued.cat', '6731_final_phot_valued.cat']
high_z=['6214_final_phot_valued.cat', '15106_final_phot_valued.cat',  '10755_final_phot_valued.cat', '13733_final_phot_valued.cat', '10920_final_phot_valued.cat', '28721_final_phot_valued.cat', '13945_final_phot_valued.cat', '12573_final_phot_valued.cat', '15260_final_phot_valued.cat', '5330_final_phot_valued.cat', '8103_final_phot_valued.cat', '13380_final_phot_valued.cat', '15127_final_phot_valued.cat', '16628_final_phot_valued.cat']

'''

# reading the data
def open_file(filename, output):
    fdata=open(my_path+filename, 'r')
    
    for line in fdata.readlines(): 
        if "#" in line:
            pass
        else:
            f=open(my_path+output,'a')
            f.write(line)
            f.close()
    fdata.close()

# outputting low and high z files
for i in xrange(len(low_z)):
    open_file(low_z[i], 'Sogras_Galclean_Low_z.txt')

for i in xrange(len(high_z)):
    open_file(high_z[i],'Sogras_Galclean_High_z.txt')
