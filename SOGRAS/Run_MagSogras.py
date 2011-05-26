#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import sys
import os
import glob

my_path='/home/maria/Area_de_Trabalho/SOAR/clusters_sogras/'
files=os.listdir(my_path)

'''
# for all 18 clusters files
clusters= [os.listdir(my_path+files[i]) for i in xrange(len(files))]
clusters_catalog =[]
for i in xrange(len(clusters)):
    for j in xrange(len(clusters[0])):
        if '.cat' in clusters[i][j]:
            clusters_catalog.append(clusters[i][j])
#print clusters_catalog

for i in range(len(files)):
    os.system('python MagSogras.py '+my_path+files[i]+'/'+clusters_catalog[i]+' 1')
'''

# for files in low and high z
my_path='/home/maria/Area_de_Trabalho/SOAR/'
os.system('python MagSogras.py '+my_path+'join_cat_sogras/SograsLow_z.txt'+' 0.5')


