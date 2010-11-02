#! /usr/bin/env python

# Start: Nov 1, 2010
# By ME, modifications by MSSG

# To do full imcat processing of Great10 files in one go
# Best run with:  nohup python pipeline_forG10_files.modified.py &

from os import *
from glob import *

###### SEXTRACTOR LOOP*****************************************************
print "____________________________________________________SEXTRACTION"
list_all=listdir(getcwd())                               # list of all files in current directory 
fitlist= glob('*.fits')                                  # get all files with extension .fits in the current directory 
sexname=[i.replace(".fits",".sexcat") for i in fitlist]  # create the files with extention .sexcat

for i in sexname:
  if i in list_all:
    print 'File '+str(i)+' already exists, skipping SExtraction of '+str(i.replace('.sexcat','.fits'))+'$...'
  else:
    system('sex -c STEP2.sex '+str(i.replace('.sexcat','.fits')))       #system is a function of os module to type terminal commands 
    system('mv tempcat '+str(i))
    print 'File '+str(i.replace('.sexcat','.fits'))+' SExtracted to '+str(i)+'!'
    
# SEXCAT TO MARUSA LOOP***********************************************
print "_______________________________________________SEXCAT TO MARUSA"
sexcatlist= glob('*.sexcat')                                     # getting in the current directory all files with extention .sexcat
catname=[i.replace(".sexcat",".cat") for i in sexcatlist]        # create the files with extention .cat

for i in catname:
  if i in list_all:
    print 'File '+str(i)+' already exists, skipping SExtraction of '+str(i.replace('.cat','.sexcat'))+'$...'
  else:
    system('sex2imcat_marusa.pl -z 25 -o '+str(i)+' '+str(i.replace('.cat','.sexcat'))) 
    print 'File '+str(i.replace('.cat','.sexcat'))+' SExtracted to '+str(i)+'!'

# SHAPEST LOOP********************************************************
print "________________________________________________________SHAPEST"
list_all_master=listdir('../mastercats')  # go up one directory and list all files in /mastercats (I'm guessing the /mastercat is a directory above the current, right)
#list_all_master=listdir('mastercats')  # go to /mastercats in current directory and list all files 
catlist= glob('*.cat')
mastercatname=[i.replace(".cat",".master.cat") for i in catlist]

for i in mastercatname:
  if i in list_all_master:
    print 'Mastercat file exists, skipping shapest for '+str(i.replace('.cat','.sexcat'))
  else:
    system('shapest.pl '+str(i.replace('.master.cat','.cat'))) 
    print 'File '+str(i)+' created from '+str('.master.cat','.cat')+'!'
system('rm *got*')


############### If you need to interrupt the processes, use
############### move_all_done_files.py at this step, kill all related
############### imcat/python procs, and rerun this code
 
# REMOVE FLAGGED OBJECTS FOR SHEAR************************************
print "___________________________________________________REMOVE FLAGS"
mastercatlist= glob('*.master.cat') 
cleancatname=[i.replace(".master.cat","clean.cat") for i in mastercatlist]

for i in cleancatname:
  system(r'lc -i "%flags 1 <" < '+str(i.replace("clean.cat",".master.cat"))+' > '+str(i))
  print i

# PSFCORRECT LOOP*****************************************************
print "_____________________________________________________PSFCORRECT"
cleancatlist=[i for i in glob('*.clean.cat') if i.find('objs')] # list of .clean.cat files with word objs in its name 
cleancatlist.sort()                                             # put the list in alphabetic order 
psffilelist=[i for i in list_all if i.find('psfs') != -1]       # list of files with word psfs in its name
psffilelist.sort()                                              # put the list in alphabetic order

for i in xrange(len(cleancatlist)):                             # start loop
	system('psfcorrect.pl -clean -p 4 -q 4 '+str(psffilelist[i])+' '+str(cleancatlist[i]))
	print cleancatlist[i]
