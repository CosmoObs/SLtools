#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import sys
import os
import glob

# path where are the cluster's files
my_path='/home/maria/Area_de_Trabalho/SOAR/clusters_sogras/'

files=os.listdir(my_path)

# getting .cat files -- photometry 
clusters= [os.listdir(my_path+files[i]) for i in xrange(len(files))]
clusters_catalog =[]
for i in xrange(len(clusters)):
    for j in xrange(len(clusters[0])):
        if '.cat' in clusters[i][j]:
            clusters_catalog.append(clusters[i][j])

# separating between low and high z
low_z=['16177_final_phot.cat', '15173_final_phot.cat', '13733_final_phot.cat', '6731_final_phot.cat']
high_z=['6214_final_phot.cat', '15106_final_phot.cat',  '10755_final_phot.cat', '13733_final_phot.cat', '10920_final_phot.cat', '28721_final_phot.cat', '13945_final_phot.cat', '12573_final_phot.cat', '15260_final_phot.cat', '5330_final_phot.cat', '8103_final_phot.cat', '13380_final_phot.cat', '15127_final_phot.cat', '16628_final_phot.cat']

def open_file(filename,filedir,output):
    fdata=open(my_path+filedir+filename, 'r')
    
    for line in fdata.readlines(): 
        if "#" in line:
            print line
            pass
        else:
            f=open(my_path.replace('clusters_sogras','join_cat_sogras')+output,'a')
            f.write(line)
            f.close()
    fdata.close()

# outputing in a specific directory
for i in xrange(len(low_z)):
    open_file(low_z[i], low_z[i].replace('_final_phot.cat','/'),'SograsLow_z.txt')

for i in xrange(len(high_z)):
    open_file(high_z[i], high_z[i].replace('_final_phot.cat','/'),'SograsHigh_z.txt')


