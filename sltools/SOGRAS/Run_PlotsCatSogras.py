#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import sys
import os
import glob
import PlotsCatSogras as PCS
import re
import numpy as np

# path where are my files
my_path='/home/maria/Area_de_Trabalho/SOAR/Galclean_Analisys_Final/'

files=os.listdir(my_path)

cat_ord, cat_no_ord = [],[]
for i in xrange(len(files)):
    if re.search('_z.txt', files[i]):
        cat_no_ord.append(files[i])                        # list of ordenated cat
    elif re.search('_z_ordenado.txt', files[i]):
        cat_ord.append(files[i])                     # list of no odernated cat

name_of_result_dir = ['Final_Results_Galclean_1plus2sersic', 'Final_Results_Galclean_1sersic', 'Final_Results_Galclean_2sersic']

#for i in name_of_result_dir:
#    os.system('mkdir '+i)

if sys.argv[1]=='1':
    save_path=my_path+name_of_result_dir[1]+'/'
    #filename= cat_ord[0]
elif sys.argv[1]=='2':
    save_path=my_path+name_of_result_dir[2]+'/'
    #f=cat_ord[0].replace('.txt'.'')
elif sys.argv[1]=='0':
    save_path=my_path+name_of_result_dir[0]+'/'

############### mudar pra low e highz 
f= 'Low_z'   #cat_ord[0].replace('.txt','')
#f= 'High_z'

sersic_1, ra1, dec1, tot_mag1, R_e1, exp1, ax_rat1, pos_ang1, disk_box1, sersic_2, ra2, dec2, tot_mag2, R_e2, exp2, ax_rat2, pos_ang2, disk_box2 = np.loadtxt( cat_ord[0], usecols=(0, 1, 2, 3, 4, 5, 6 ,7 ,8 ,9 ,10 ,11 ,12 ,13 ,14 ,15 ,16 ,17), unpack=True)
################

def filter_1(s1, s2):
    filter_1=(s1==1)*(s2==1)   # if sersic _1 and _2 == 1   for 2 sersic analisys
    return filter_1

def filter_2(s1, s2, m1, m2):
    filter_2=(s1==1)*(s2==0)#*(m2!=-99.999)   # if sersic _1==1 and _2 == 0 and tot_mag2==-99.00   for 1 sersic analisys of mag e individual plots
    return filter_2

# filtering sersics
k1=filter_1(sersic_1, sersic_2)
k2=filter_2(sersic_1, sersic_2, tot_mag1, tot_mag2)

# join the 'columns'
m = np.append(tot_mag1[k1], tot_mag2[k1])
tot_mag = np.append(m, tot_mag1[k2])

d = np.append(disk_box1[k1], disk_box2[k1])
disk_box = np.append(d, disk_box1[k2])

r = np.append(R_e1[k1], R_e2[k1])
R_e = np.append(r, R_e1[k2])

e = np.append(ax_rat1[k1], ax_rat2[k1])
ax_rat = np.append(e, ax_rat1[k2])

p = np.append(pos_ang1[k1], pos_ang2[k1])
pos_ang = np.append(p, pos_ang1[k2])

n = np.append(exp1[k1], exp2[k1])
exp = np.append(n, exp1[k2])

### ---------------- combinations to plot low z, high z, join, separated, order, no order, etc ---------
# 1 -> maior ns , 2-> menos ns
'''
PCS.hist_PA(pos_ang1[k1], pos_ang2[k1], tot_mag1[k1], tot_mag2[k1],save_path, f)    
PCS.hist_R(R_e1[k1], R_e2[k1], tot_mag1[k1], tot_mag2[k1],save_path, f)
PCS.hist_ns(exp1[k1], exp2[k1], tot_mag1[k1], tot_mag2[k1],save_path, f)
PCS.hist_ax_ratio(ax_rat1[k1], ax_rat2[k1], tot_mag1[k1], tot_mag2[k1],save_path, f)
PCS.hist_disk_boss(disk_box1[k1], disk_box2[k1], tot_mag1[k1], tot_mag2[k1],save_path, f)

PCS.scatter_ns(tot_mag1[k1], exp1[k1], exp2[k1], save_path, f)
PCS.hist_disk_bulge(tot_mag1[k1], tot_mag2[k1],save_path, f)
PCS.hist_mag_soma(tot_mag1[k1+k2], tot_mag2[k1+k2],save_path, f)
PCS.hist_mag(tot_mag1[k1], tot_mag2[k1],save_path, f)
'''

''' # 
PCS.hist_PA(pos_ang, pos_ang, tot_mag, tot_mag,save_path, f)    
PCS.hist_R(R_e1, R_e, tot_mag, tot_mag,save_path, f)
PCS.hist_ns(exp, exp, tot_mag, tot_mag,save_path, f)
PCS.hist_ax_ratio(ax_rat, ax_rat, tot_mag, tot_mag,save_path, f)
PCS.hist_disk_boss(disk_box, disk_box, tot_mag, tot_mag,save_path, f)
PCS.hist_mag(tot_mag, tot_mag,save_path, f)
'''

'''  # 2 sersics
PCS.scatter_ns(tot_mag1[k1], exp1[k1], exp2[k1], exp1[k2], exp2[k2], save_path, f)
PCS.hist_disk_bulge(tot_mag1[k1], tot_mag2[k1],save_path, f)
PCS.hist_mag_soma(tot_mag1[k1+k2], tot_mag2[k1+k2],save_path, f)
'''
PCS.hist_mag_soma_tot_mag(tot_mag1[k1], tot_mag2[k1], tot_mag1[k2], tot_mag, save_path, f)

