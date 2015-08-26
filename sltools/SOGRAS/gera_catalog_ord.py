#!/usr/bin/env python
# ====================================================================
# Authors:
# Maria Elidaiana da S. Pereira - mariaeli@cbpf.br
# ====================================================================

import sys
import numpy as np

filename = sys.argv[1]
sersic_1, ra1, dec1, tot_mag1, R_e1, exp1, ax_rat1, pos_ang1, disk_box1, sersic_2, ra2, dec2, tot_mag2, R_e2, exp2, ax_rat2, pos_ang2, disk_box2 = np.loadtxt(filename, usecols=(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18), unpack=True)

#----- oredamento do cs82-------------
#idx1 = (tot_mag2!=-99.999)*(exp1 > exp2)  # ns2 menor que ns1
#idx2 = (tot_mag2!=-99.999)*(exp1 <= exp2) # ns2 maior ou igual a ns1

#----- ordenamento do sogras-----------
idx1 = (exp1 > exp2)  # ns2 menor que ns1
idx2 = (exp1 <= exp2) # ns2 maior ou igual a ns1

ns_ordenado = np.array(zip(append(sersic_1[idx1],sersic_2[idx2]), append(ra1[idx1],ra2[idx2]), append(dec1[idx1], dec2[idx2]), append(tot_mag1[idx1],tot_mag2[idx2]),append(R_e1[idx1],R_e2[idx2]),append(exp1[idx1], exp2[idx2]),append(ax_rat1[idx1], ax_rat2[idx2]),append(pos_ang1[idx1],pos_ang2[idx2]), append(disk_box1[idx1], disk_box2[idx2]), append(sersic_2[idx1],sersic_1[idx2]), append(ra2[idx1],ra1[idx2]), append(dec2[idx1], dec1[idx2]), append(tot_mag2[idx1],tot_mag1[idx2]),append(R_e2[idx1],R_e1[idx2]),append(exp2[idx1], exp1[idx2]),append(ax_rat2[idx1], ax_rat1[idx2]),append(pos_ang2[idx1],pos_ang1[idx2]), append(disk_box2[idx1], disk_box1[idx2])))
         
np.savetxt(str(argv[1]).replace('.txt','_')+'ordenado.txt', ns_ordenado, fmt='%i %10.6f %10.6f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %i %10.6f %10.6f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f')


