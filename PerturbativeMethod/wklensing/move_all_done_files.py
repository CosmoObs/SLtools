#! /usr/bin/env python

# Start: Nov 2, 2010
# By MSSG

# To move all done files into donefiles dir, if you accidentally interrupt the processes,

from os import *
from glob import *

mastercatlist=glob('*.master.cat') # Get all master cats -- this defines the done files

for file_i in mastercatlist:  
  fitsname_i = str(file_i.replace(".master.cat",".fits")) # Get back the fits file name for the given mastercat file
  sexname_i = str(file_i.replace(".master.cat",".sexcat")) # Get back the sexcat file name for the given mastercat file
  catname_i = str(file_i.replace(".master.cat",".cat")) # Get back the cat file name for the given mastercat file
  
# Now move them all to the donefiles subdir
  system("mv "+ file_i+ " donefiles ")  
  system("mv "+ fitsname_i+ " donefiles ")  
  system("mv "+ sexname_i+ " donefiles ")
  system("mv "+ catname_i+ " donefiles ")  
