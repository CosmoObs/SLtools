#!/usr/local/bin/perl -w
# ======================================================================
#+
$usage = "

NAME
        shearprofile.pl

PURPOSE
        Read in an object catalogue, offset it and compute tangential shear.
        Also compute tangential shear after rotating the galaxies by
        45 degrees (the 'B modes'). Bin the galaxies according to some scheme
        and output the shear profiles as 3 column text (x,y,e).


USAGE
        shearprofile.pl [flags] [options] \$project.cat \$other.cat ...

FLAGS
        -errors   Include error propagation
        -u        Print this message
        -h        Ignore catalogue header (marked with \#)

INPUTS
        \$project.cat   Object catalogue (multi-column text)
        \$other.cat     Another object catalogue (multi-column text)

OPTIONAL INPUTS
        -oE     Efile   Write output to \"file\" (def=\$project.Eprofile.txt)
        -oB     Bfile   Write output to \"file\" (def=\$project.Bprofile.txt)
        -cx       f     X position of cluster center
        -cy       f     Y position of cluster center
        -rmin     f     Minimum radius [0.0]
        -rmax     f     Maximum radius [10.0]
        -nbin     i     No. of radius bins [20.0]
        -xcol     i     Use ith column for x positions (def=0)
        -ycol     i     Use ith column for y positions (def=1)
        -e1col    i     Use ith column for e1 (def=-4)
        -e1errcol i     Use ith column for e1err (def=-3)
        -e2col    i     Use ith column for e2 (def=-2)
        -e2errcol i     Use ith column for e2err (def=-1)
        -nucol    i     Use ith column for nu 
        -radcol   i     Use ith column for size
        -lc             Output in lc catalog format
        -wtcutoff f     Use this for the weight def'n cutoff 
        -clus     s     Which cluster to do the redshift correction for

OUTPUTS
        \$project.Eprofile.txt     or explicitly named output catalogue file
        \$project.Bprofile.txt

OPTIONAL OUTPUTS

COMMENTS

EXAMPLES

BUGS
 - no option to reduce verbosity

REVISION HISTORY:
  2006-07-26 Started Marshall (KIPAC)
  2007-05-30 Added cluster center offset, changed to arbitrary coords Applegate (KIPAC)
  2009-04-06 Added integral table feature (jpd)

\n";
#-
# ======================================================================

# $\="\n";

use Getopt::Long;
GetOptions("oE=s", \$Eoutfile,
           "oB=s", \$Boutfile,
	   "cx=f", \$cx,
	   "cy=f", \$cy,
           "rmin=f", \$rmin,
           "rmax=f", \$rmax,
           "nbin=i", \$nbin,
           "xcol=i", \$xcol,
           "ycol=i", \$ycol,
           "e1col=i", \$e1col,
           "e2col=i", \$e2col,
           "e1errcol=i", \$e1errcol,
           "e2errcol=i", \$e2errcol,
 	   "zcol=f", \$zcol,
           "nucol=i", \$nucol,
           "radcol=i", \$radcol,
           "errors", \$errors,
           "h", \$header,
           "u", \$help,
	   "lc", \$lc,
	   "verb=i", \$verb,
	   "wtcutoff=f", \$wtcutoff,
	   "clus=s", \$clus
           );


(defined($help)) and die "$usage\n";
$num=@ARGV;
($num == 0) and die "$usage\n";

$sensible = 0;
(defined($Eoutfile)) or ($sensible = 1);
($num>1) and $sensible=1;

print " Eoutfile = $Eoutfile \n";
print " Boutfile = $Boutfile \n";

# print $num;
print "\n";
print "\n";
print @ARGV;
print "\n";
print "\n";

print " Hit me wid it!";
 <STDIN>;


# Default verbosity
(defined($verb)) or ($verb = 0);
$v = $verb ;

# Default cluster
(defined($clus)) or ($clus = 'des');

print "Cluster = $clus";


# Default column numbers:
(defined($zcol)) or ($zcol = 0); # Default, no z correc
(defined($xcol)) or ($xcol = 0);
(defined($ycol)) or ($ycol = 1);
(defined($e1col)) or ($e1col = -4);
(defined($e1errcol)) or ($e1errcol = -3);
(defined($e2col)) or ($e2col = -2);
(defined($e2errcol)) or ($e2errcol = -1);


print " e1col, e2col = $e1col, $e2col \n";
# <stdin>;

(defined($nucol)) or ($nucol = -1);
(defined($radcol)) or ($radcol = -1);

# Default binning parameters:
(defined($cx)) or ($cx = 0);
(defined($cy)) or ($cy = 0);
(defined($rmin)) or ($rmin = 0.0);
#$rmin = $rmin*60.0; 
(defined($rmax)) or ($rmax = 10.0);
#$rmax = $rmax*60.0; 
(defined($nbin)) or ($nbin = 20);
$dr = ($rmax - $rmin) / $nbin;

print " \n $dr \n ";

(defined($wtcutoff)) or ($wtcutoff = 40.0);



################### Some initzns

$numgals =0;
$totwt=0;

$countstep = 500;
 (defined($verb)) or ($verb = 0);

# Loop over catalogues:

while (defined($file = shift)){

# Set up binning grid:

    $root = $file;
    print " root = $root  \n";

$wtless =0;
    
    for ($i = 0; $i < $nbin; $i++){
      $rbin[$i] = $rmin + ($dr * ($i+0.5) ) ; # Center of i'th bin with Lower edge  $rmin + ($dr * $i) ) of i'th bin -- Upper edge is $rmin + $dr * ($i+1)
      $Ebin[$i] = 0.0;
      $Ebinwtd[$i] = 0.0;
      $Ebinerr[$i] = 0.0;
      $Bbin[$i] = 0.0;
      $Bbinwtd[$i] = 0.0;
      $Bbinerr[$i] = 0.0;
      $Eorigbin[$i] = 0.0;
      $Eorigbinerr[$i] = 0.0;
      $Borigbin[$i] = 0.0;
      $Borigbinerr[$i] = 0.0;
      $area[$i] = 0.0;
      $wtsum[$i] = 0.0;
  }

    print "  \n rbin[i]   =  @rbin \n"; 
#    <stdin>;

    open (IN, $file) or die "$file: $!";

    print  "\n Reading data from $file \n";

# Sort out sensible filename:
    if ($sensible) {
       $root = $file;
       $root =~ s/(.*)\..*/$1/;
       $Eoutfile = $root.".Eprofile.txt";
       $Boutfile = $root.".Bprofile.txt";
    }

       $wtoutfile = $root.".weights.txt";

#    open (WTOUT, " > $wtoutfile") or die "can't fork: $!";  # to write to weights file

    print " Eoutfile = $Eoutfile \n";
    print " Boutfile = $Boutfile, and wait: \n";

#    <stdin>;
    
    $numdensityfile = "numdensityProfile.txt";

# Count objects:

    $count = 0;
    $headcount = 0;

    (defined($header)) and (<IN>);

    print " Starting into file \n" ;


#    $wtcutoff = 40;
    $wtmin = 5;

############################ Step through lines of catalogue,
    while (<IN>){
      $headcount ++;
      next if ($headcount < 7);
      chomp;
# dealing with header lines,
      if (/^\#/ or /^$/){
        next;
      } else {
        @line = split;
      }
# and working out appropriate object information:



# First get all the quantities:
      $id = $line[0];
      $x = $line[$xcol];
      $y = $line[$ycol];
      $e1 = $line[$e1col];
      $e2 = $line[$e2col];
      $z = $line[$zcol];
      if ($zcol == 0) { $z = 0};  # Default, no z correc
      if ($zcol < 0 ) { $z = -$zcol}; 

      if ($nucol != -1)  {     $nu = $line[$nucol]} else {$nu=100000};

      if ($radcol != -1)  {      $rad = $line[$radcol]} else {$rad = 1} ;


      if ($v > 0)  {     print " ************************** Obj num:  $id \n";
			 print  "x   y    $x   $y  \n";
			 print  "e1  e2   $e1  $e2 \n";
		     };
      
# First compute radius and hence bin number:
      $xrel = $x - $cx;
      $yrel = $y - $cy;
      $r = sqrt($xrel*$xrel + $yrel*$yrel);
      $i = int(($r - $rmin) / $dr); # Bin num
      $phi = atan2($yrel,$xrel);
      $cos2phi = cos(2.0*$phi);
      $sin2phi = sin(2.0*$phi);




      if ($r < $rmin) {
#	  print " r = $r \n ";
	  next ;  # Skip stuff below edge
      }

      if ($r < $rmin)  {  	  print   "r, i'th bin,  phi    $r $i  $phi \n"; 
				  print " E[i] =   $E \n" ; 
				  print " Eorig[i] =   $Eorig \n" ; 			
#				  <stdin> ; 
			      } ;
      


      if (abs($e1) < 1e-6 && abs($e2) <1e-6) {
#	  print " $e1, $e2  \n ";
	  next ;  # Skip stuff below edge
      }


      if ($e1 ==0 && $e2 ==0) {  	  print   "r, i'th bin,  phi    $r $i  $phi \n"; 
				  print " E[i] =   $E \n" ; 
				  print " Eorig[i] =   $Eorig \n" ; 			
#				  <stdin> ; 
			      } ;
      


# Now eval shrfac in distances.mw for the midpoint redshifts of each range (i.e. include Dls/Ds factor)



# jpd

       


      if ($v > 0)  { 
	  print "\n dr = $dr \n";


	  print   "r, i'th bin,  phi    $r $i  $phi \n";
	  print  "cos2phi, sin2phi :   $cos2phi, $sin2phi \n"
	  };





# Calculate tangential shear and add it to the (weighted) sum:     

#      print "i, nucol, nu, E, E / (1/nu + 1),  Ebin[i], Ebinunwtd[i] = $i, $nucol, $nu, $E, $E / (1/$nu + 1),  $Ebin[$i], $Ebinunwtd[$i] ";
#     <stdin>;
      
      $Eorig = -($e1*$cos2phi+$e2*$sin2phi)  ;

#  <STDIN>;

 
        $nbinn[$i]++;

#      if ($nucol != -1) { # if nucol def'd
	
      if ($nu < $wtcutoff) {
	  $w = ($nu - $wtmin) / ($wtcutoff-$wtmin);
	  $wtless++;
#	  print " wtless = $wtless \n";
     }
      if ($nu >= $wtcutoff) {$w = 1;}    ;


      if ($i<$nbin) {
	  $numgals++;
	  $totwt+=$w;
      };

#      if ($i < $nbin) {      print WTOUT "$i  $nu   $w \n"}; # weights

	 
#     print " $wtcutoff, $nu, $w  \n";
# <stdin>;

	  $Eorigbin[$i] += $Eorig ; 	      
	  $Ebin[$i] += $E ; 	      
	  $Ebinwtd[$i] += $w*$E;
	  $wtsum[$i] += $w;

	  $Eorigbinerr[$i] += $Eorig*$Eorig; 
	  $Ebinerr[$i] += $E*$E; 

    
      if ($v > 0)  {     
	  print  "e1, cos2phi, e2 , sin2phi,  E :  $e1, $cos2phi, $e2, $sin2phi, $E\n";	  
	  print  " \n"
	  };

# Now rotate ellipticity by 45 degrees:

      $b1 =  $e2;
      $b2 = -$e1;
      
      if  ($v >0) { 
	  print " e1, e2 = $e1, $e2 \n"  ;
	  print " b1, b2 = $b1, $b2  \n" ; 	
      }; 
      

      $Borig = -($b1*$cos2phi+$b2*$sin2phi);

	  print  "x, y,  e1, e2 ,  r, E , B :  $x,$y, $e1, $e2, $r,  $Eorig, $Borig\n";	  


      $Borigbin[$i] += $Borig; 
      $Bbin[$i] += $B; 
     
      $Bbinwtd[$i] += $w*$B;
      
      $Borigbinerr[$i] += $Borig*$Borig; 
      $Bbinerr[$i] += $B*$B; 

      if  ($v >0) { 
	  print " Ebin[i] =   $Ebin[$i] \n" ; 
	  print " Bbin[i] =   $Bbin[$i] \n" ; 

	  print " Eorigbin[i] =   $Eorigbin[$i] \n" ; 
	  print " Borigbin[i] =   $Borigbin[$i] \n" ; 

#	  <STDIN> 
	  } ;



      $count++;
#      if ($count/$countstep == int($count/$countstep)){ print "At event $count  count, i, Ebin, Ebinerr =  $count, $i,  $Ebin[$i], $Ebinerr[$i] \n";  } ;

      if ($count/$countstep == int($count/$countstep)){ print "At event $count  count, i, Eorigbin, Eorigbinerr =  $count, $i,  $Eorigbin[$i], $Eorigbinerr[$i] \n";  } ;

      if ($  v>0){
	  print " i, nbinn[i] = $i, $nbinn[$i] \n";
	  <stdin>;
      }

    }
    close(IN);   #end of looping over objs

	print "Done: $count objects processed\n";    
	print "Hit me with your bestest shot, fire away:  " ;

    if ($verb > 0)  { 
	<STDIN>;
    } 
    


# Now work out averages and errors and finish!
# Open output files:
    if (defined($lc)){
	print "eout = ", EOUT;
	print " eoutfile = $Eoutfile \n";
	print "  numdensityfile: =  $numdensityfile: \n ";
	
	open (EOUT, "| lc -C -n x -n y -n e > $Eoutfile") or die "can't fork: $!";
	open (BOUT, "| lc -C -n x -n y -n e > $Boutfile") or die "can't fork: $!";
#	open (NUMOUT, "| lc -C -n x -n num  -n area -n numdens  > $numdensityfile") or die "can't fork: $!";
    } else {
	open (EOUT, ">$Eoutfile") or die "Couldn't open $Eoutfile: $!";
	open (BOUT, ">$Boutfile") or die "Couldn't open $Boutfile: $!";

#	open (NUMOUT, ">$numdensityfile") or die "Couldn't open $numdensityfile: $!";
    }


# ##################### Get final vals in bins

    for ($i = 0; $i < $nbin; $i++){
      if (defined($errors)){
        if ($Ebinerr[$i] > 0.0){
          $Ebin[$i] = $Ebin[$i] / $Ebinerr[$i];
          $Ebinerr[$i] = sqrt(1.0 / $Ebinerr[$i]);
          $Bbin[$i] = $Bbin[$i] / $Bbinerr[$i];
          $Bbinerr[$i] = sqrt(1.0 / $Bbinerr[$i]);      
      }
    } else {

	print "	\n bin, nbinn[i]  =   $i,   $nbinn[$i] \n" ;

        if ($nbinn[$i] > 1){
	    
######### E
	    $Eorigbin[$i] = $Eorigbin[$i] / $nbinn[$i];
	    $Eorigbinerr[$i] = sqrt($Eorigbinerr[$i] / $nbinn[$i] - $Eorigbin[$i]*$Eorigbin[$i]);
	    $Eorigbinerr[$i] = $Eorigbinerr[$i] / sqrt($nbinn[$i] - 1.0);

	    $Ebin[$i] = $Ebin[$i] / $nbinn[$i];
	    $Ebinerr[$i] = sqrt($Ebinerr[$i] / $nbinn[$i] - $Ebin[$i]*$Ebin[$i]);
	    $Ebinerr[$i] = $Ebinerr[$i] / sqrt($nbinn[$i] - 1.0);

	    $Ebinwtd[$i] = $Ebinwtd[$i] / $wtsum[$i];


########### B	    
	    $Borigbin[$i] = $Borigbin[$i] / $nbinn[$i];
	    $Borigbinerr[$i] = sqrt($Borigbinerr[$i] / $nbinn[$i] - $Borigbin[$i]*$Borigbin[$i]);
	    $Borigbinerr[$i] = $Borigbinerr[$i] / sqrt($nbinn[$i] - 1.0);

	    $Bbin[$i] = $Bbin[$i] / $nbinn[$i];
	    $Bbinerr[$i] = sqrt($Bbinerr[$i] / $nbinn[$i] - $Bbin[$i]*$Bbin[$i]);
	    $Bbinerr[$i] = $Bbinerr[$i] / sqrt($nbinn[$i] - 1.0);

	    $Bbinwtd[$i] = $Bbinwtd[$i] / $wtsum[$i];

	    print "nbinn[i]-1 > 1 \n", $nbinn[$i]-1 > 1;
	    
	    print  "  nbinn[$i],  Ebin[$i],  Bbin[$i]  =   $nbinn[$i],  $Ebin[$i], $Bbin[$i] \n";
	    print  "  Ebinerr[$i],  Bbinerr[$i] =  $Ebinerr[$i], $Bbinerr[$i]  \n" ;
	    
	    $area[$i] =  (($i*$dr -$rmin + $dr )**2 - ($i*$dr - $rmin)**2)/10**5;
	    
#	  printf "\n  $i ,   $dr , $rmin ,  $dr,  %f, %f   ,  %f , %f \n", $i*$dr -$rmin, ($i*$dr -$rmin)**2,  , ($i*$dr -$rmin + $dr ), ($i*$dr -$rmin + $dr)**2 ;
	    
	    print "\n      Low edge =", $i*$dr - $rmin, "  High edge = ", $i*$dr -$rmin + $dr, "     Area/10^5 = ", $area[$i];
	    print "\n\n\n\n\n";
	    
#      $ndens = $nbinn[$i]/$area[$i] ;
	    
	} # If nbinn>1 loop
    } # else loop
      
#       print STDERR "i, rbin, Ebin, Ebinerr: $i, $rbin[$i] $Ebin[$i], $Ebinerr[$i]\n";
#   Print out x y e text into files:

      print EOUT "$rbin[$i]  $Eorigbin[$i]  $Eorigbinerr[$i]  $Ebin[$i]  $Ebinerr[$i]  $numgals   $Ebinwtd[$i]  $Ebinerr[$i]   $totwt \n";
      print BOUT "$rbin[$i]  $Borigbin[$i]  $Borigbinerr[$i]  $Bbin[$i]  $Bbinerr[$i] $numgals   $Bbinwtd[$i]  $Bbinerr[$i]   $totwt \n";


#     print "\n  $rbin[$i]  $nbinn[$i]  $area[$i] $nbinn[$i]/$area[$i]  $ndens \n";

#      print NUMOUT "$rbin[$i]  $nbinn[$i]  $area[$i] $ndens  \n";
  }   # End running over bins

      print "\n";


    close(EOUT) or die "Couldn't close $Eoutfile: $!";
    close(BOUT) or die "Couldn't close $Boutfile: $!";
    
#    close (WTOUT) ;

    print STDERR "E-mode profile in $Eoutfile\n";
    print STDERR "B-mode profile in $Boutfile\n";
    print STDERR "Num density profile in $numdensityfile\n";
    

}

# ======================================================================
