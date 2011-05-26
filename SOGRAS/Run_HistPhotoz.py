#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import sys
import os
import glob

my_path='/home/maria/Area_de_Trabalho/SOAR/Photoz_Sogras/'
files=os.listdir(my_path)

# --- combinations for individual and low and high z files

zcatalogs = []
for i in xrange(len(files)):
    if '.cat' in files[i]:
            zcatalogs.append(files[i])

for i in xrange(len(zcatalogs)):
    os.system('python HistPhotoz.py '+zcatalogs[i])    #+my_path+'/'+zcatalogs[i])

'''
my_path='/home/maria/Area_de_Trabalho/SOAR/'
os.system('python HistPhotoz.py '+my_path+'/SograsLow_z.txt'+' 0.5')

'''
