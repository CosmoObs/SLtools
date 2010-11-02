# Code to get E and B mode components for each galaxy 
# by ME, BALM, MSSG
# Start: Oct 2010

from __future__ import division
from sys import argv
from math import sqrt, atan2, cos, sin
from numpy import *
from numpy import trapz #-----Method chosen for the numerical integration-----
import time #----For calculating the running time


#--------------------------------------FIRST PART: DEFINITION OF CONSTANTS AND VARIABLES------------------------------------------

#------Timing the program--------

t0 = time.time()


#-----------Cluster Information------------------- 


zl = 0.33 #----Cluster redshift, hardcode for now


#The coordinates of the center of the cluster are not used since the input file already gives the relative positions of the galaxies. For completeness:
#-------------------------------------
# For the DES Simln: x_c = 12005.0, y_c = 14518.0


#-----------Cosmological Information---------------


H0 = 2.300951306*10E-18 # In 1/sec
gpctocm = 3.08568025*10E27 #---Conversion factor from Gpc to centimeters
OmegaM0 = 0.27 # Omega_m0
OmegaDE0 = 0.73 # Omega_Lambda0
c = 2.99792458*10E10 #----In centimeters/sec


#-------------Geometrical variables for the positions of the galaxies--------------

x_rel, y_rel, r, phi = [], [], [], []  # Define empty vectors

#------------Redshift of the sources and beta factor---------------------

betas, zs, indices = [], [], []

#------------Shape variables---------------------------

gamma1, gamma2, Eorig, Borig, E, B = [], [], [], [], [], []

#---------------------SECOND PART: EXTRACTING INPUT FILE INFORMATION AND CALCULATION OF THE MODES 'Eorig' AND 'Borig'-----------------------------

#---------------------Incorporation of the input file information-----------------

X = open(argv[1]).readlines() # Open the file


for i in range(len(X)):
	if float(X[i].split()[4]) > 0.35:#-----This line establishes the redshift criterium for selecting the galaxies---
	 x_rel.append(float(X[i].split()[0]))
	 y_rel.append(float(X[i].split()[1]))
	 gamma1.append(float(X[i].split()[2]))
	 gamma2.append(float(X[i].split()[3]))
	 zs.append(float(X[i].split()[4]))

	else:pass  # i.e. skip this object

print len(zs) #------ To check the number of galaxies in the redshift range selected above------


#---------------------Calculation of 'Eorig' and 'Borig'---------------------------


for i in range(len(zs)):
  r.append(sqrt(x_rel[i]**2 + y_rel[i]**2)) # Distance from center of the cluster
  phi.append(atan2(y_rel[i],x_rel[i]))      # Angle wrt +x axis
  Eorig.append(-(gamma1[i]*cos(2*phi[i])) - (gamma2[i]*sin(2*phi[i]))) # E signal
  Borig.append(-(gamma1[i]*sin(2*phi[i])) + (gamma2[i]*cos(2*phi[i]))) # B signal

#--------------------------------------THIRD PART: CALCULATION OF THE DISTANCE FACTOR 'beta'----------------------------------------

#------------------Definition of the Hubble function-------------------

def h(z):
 return sqrt(OmegaM0*(1+z)**3+OmegaDE0)

#------------------Numerical integration of the comoving distance ----------------------------------------

step=10E-4 #-----Integration step------

def integrand(z):  # Define function
  discrete=[]      # Empty vetor
  for i in xrange(1, int(z/step)+1):  # Step all the way to max z
    discrete.append(c/(H0*gpctocm*h(i*step)))
  return discrete

def Dm(z): # Define the basic measure of distance 
  return  trapz(integrand(z), dx=step) # Do the integ'n


#-----------Calculation of the beta factor (with definitions of the angular diameter distances inside)-------------------------------

def beta(z):
 Da=(1./(1. + z))*Dm(z)          # General def'n of the angular diameter distance
 Da12=(Dm(z) - Dm(zl))/(1 + z)   # Dls
 Dls=Da12
 Ds=Da                           # Ds
 beta_value=Dls/Ds               # beta
 return beta_value



#--------------------------------------FOURTH PART: GETTING 'E' AND 'B'----------------------------------------

#----------Getting the beta factor for galaxies behind the cluster---------------------------

for i in range(len(zs)): 
	betas.append(beta(zs[i]))  # Get beta for each galaxy

#----------Calculating and Printing 'E' and 'B'--------------------------------------------


for i in range(len(zs)):
	E.append(Eorig[i]/betas[i])  # Get the E and B modes by dividing through by beta factor
	B.append(Borig[i]/betas[i])
 
#---------Columns in the output file are, in order: relative x position, relative y position, gamma1, gamma2, source redshift, Eorig, Borig, E = Eorig/beta and B = Borig/beta-------------------------

file=open('des_clustersim_reduced_out.txt','w')  # Print to file
for i in range(len(zs)):
	file.write("%-6s" % str(x_rel[i])+'     '+str(y_rel[i])+'     '+str(gamma1[i])+'     '+str(gamma2[i])+'     '+str(zs[i])+'     '+str(Eorig[i])+'     '+str(Borig[i])+'     '+str(E[i])+'     '+str(B[i])+'\n')
file.close()


#------Printing total time it took to get through this ---------

print time.time() - t0, "seconds"
