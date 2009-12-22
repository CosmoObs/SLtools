#!/usr/bin/env python
import math
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt


import get_nfw_concentration 

n=50000

c=range(n)

mass=1e14
redshift=0
	

get_nfw_concentration.read_config("get_nfw_concentration.conf")

for i in range(n):
	c[i] = get_nfw_concentration.c(mass, redshift, 0)

# statistics

c_mean=np.mean(c)

# the histogram of the data

n, bins, patches = plt.hist(c, 50, normed=1, alpha=0.75)

plt.xlabel('c')
plt.ylabel('P(c,M,z)')
plt.title(r'$\mathrm{P.D.F \,of \, c},\, M_{200}=' + str(mass) +'M_{\odot}/h,\, z=' + str(redshift) + ',\, <c>='+ str(c_mean) +',\ \sigma_{log_c}=0.1$')
plt.axis([0, 12, 0, 0.5])
plt.grid(True)

plt.show()

