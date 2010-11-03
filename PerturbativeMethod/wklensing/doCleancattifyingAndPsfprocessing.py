#! /usr/bin/env python

# Start: Nov 2, 2010
# By MSSG



from os import *
from glob import *

# Get list of all files in current directory 
total_file_list=listdir(getcwd())

# Make Cleancat
print "___________________________________________________Make Clean Cats"
mastercatlist= glob('*.master.cat') 

cleancatname=[file_i.replace("master.cat","clean.cat") for file_i in mastercatlist]


for file_i in cleancatname:
  cleanfile_i = str(file_i)
  masterfile_i = str(file_i.replace("clean.cat","master.cat"))
  if cleanfile_i in total_file_list:     # If cleanfile already exists    
    print 'File '+cleanfile_i+' already exists, skipping '
    print " Num objs in " + masterfile_i + " = "
    #    system('lc -c < '+ masterfile_i)
    print " Num objs in " + cleanfile_i + " = "
    #    system('lc -c < '+ cleanfile_i)
    print
  else:
    system('lc -i "%flags 1 <" < '+ masterfile_i+' > '+ cleanfile_i)
    print file_i


# PSFCORRECT LOOP*****************************************************
print "_____________________________________________________PSFCORRECT"
objslist=[file_i for file_i in glob('*.clean.cat') if file_i.find('objs')!= -1] # list of .clean.cat files with word objs in its name 
objslist.sort()                                             # put the list in alphabetic order 
psflist=[file_i for file_i in glob('*.clean.cat') if file_i.find('psfs') != -1]       # list of files with word psfs in its name
psflist.sort()                                              # put the list in alphabetic order

print objslist
print
print psflist

for i in xrange(len(objslist)):                             # start loop
  gammafile_i = str(objslist[i].replace("clean.cat","clean.psfcorr.4.gamma.cat"))
  if gammafile_i in total_file_list:
    print 'File '+gammafile_i+' already exists, skipping '
  else:
    system('psfcorrect.pl -clean -p 4 -q 4 '+str(psflist[i])+' '+str(objslist[i]))
    print objslist[i]

for i in xrange(len(objslist)):                             # start loop
  reddfile_i = str(objslist[i].replace("clean.cat","redd"))
  gammafile_i = str(objslist[i].replace("clean.cat","clean.psfcorr.4.gamma.cat"))
  if reddfile_i in total_file_list:
    print 'File '+reddfile_i+' already exists, skipping '
  else:
    system('lc x es gamma smag fwhm < ' + gammafile_i + ' >' + reddfile_i)



