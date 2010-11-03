# Oct 27, 2010

Tmp dir for wk lensing code by M.E., BALM, and MSSG

We will change location later

############# General WL utilities (for DES, ClusterStep, LBT etc.)

--- shearprofile.py - Extracting E and B mode signal for each galaxy
 at a given z, and at relative transverse distance from a cluster
 center

--- shr.pl - legacy perl code that is now superceded by shearprofile.py 

############# G10

--- pipeline_forG10_files.py - to run Imcat all the way up and through
    shapest

--- move_all_done_files.py - if interrupted, this moves all fully proc'd files to donefiles subdir.


--- doCleancattifyingAndPsfprocessing.py - makes cleancats, runs psfcorrect, and then makes redd files --> final step


