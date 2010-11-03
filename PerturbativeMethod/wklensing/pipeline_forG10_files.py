#! /usr/bin/env python

# Start: Nov 1, 2010
# By ME, modifications by MSSG

# To do full imcat processing of Great10 files in one go

from os import *
from glob import *

###### SEXTRACTOR LOOP*****************************************************
print "____________________________________________________SEXTRACTION"
total_file_list=listdir(getcwd())                        # Get list of all files in current directory 
fitlist= glob('*.fits')                                  # Get all files with extension .fits in the current directory 
sexlist=[filename.replace(".fits",".sexcat") for filename in fitlist]  # Create list of files with extension .sexcat

### Loop over all files
for file_i in sexlist:  # Now file_i is a list of all the sexcat filenames
  fitsname_i = str(file_i.replace('.sexcat','.fits')) # Get back the fits file name for the given sexcat file
  if file_i in total_file_list:     # If SE'd file already exists    
    print 'File '+str(file_i)+' already exists, skipping SExtraction of '+ fitsname_i +'$...'
  else:
    system('sex -c STEP2.sex '+ fitsname_i )       # System is a call to the Linux terminal
    system('mv tempcat '+str(file_i))
    print 'File '+ fitsname_i +' SExtracted to '+str(file_i)+' -- darn tootin ! '
    
# SEXCAT TO MARUSA LOOP***********************************************
print "_______________________________________________SEXCAT TO MARUSA"
sexcatlist= glob('*.sexcat')                                     # Get all files with extension .sexcat in the current dir
catlist=[file_i.replace(".sexcat",".cat") for file_i in sexcatlist]   # Create list of files with extension .cat

for file_i in catlist: # Now file_i is a list of all the cat filenames

  sexname_i = str(file_i.replace('.cat','.sexcat'))
  catname_i =str(file_i)

  if file_i in total_file_list:      # If cat file already exists    
    print 'File '+ catname_i +' already exists, skipping SExtraction of '+ sexname_i + '$...'
  else:
    system('sex2imcat_marusa.pl -z 25 -o '+ catname_i + ' ' + sexname_i )
    print 'File '+ sexname_i +' made into catfile '+ catname_i + ' - yessum ! '

# SHAPEST LOOP********************************************************
print "________________________________________________________SHAPEST"
catlist= glob('*.cat')
mastercatlist=[file_i.replace(".cat",".master.cat") for file_i in catlist]

for file_i in mastercatlist:
  catname_i = str(file_i.replace('.master.cat','.cat'))
  master_i =str(file_i)
  if file_i in total_file_list:      # If mastercat file already exists    
    print 'Mastercat file exists, skipping shapest for '+ catname_i
  else:
    system('shapest.pl '+ catname_i ) 
    print 'File '+ master_i +' created from ' + catname_i + " -- yahoozy b's and g's!! " 

# system('rm *got*')  # Clean up 'got' files after processing

