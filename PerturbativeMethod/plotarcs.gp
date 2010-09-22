# Start Sept 2010
# Script to plot the variations in arcs

set pointsize 0.5

 set terminal postscript color

set out 'arcs1.ps'
plot [-3:3][-3:3] 'arcdata/arcs1.txt' u 1:2
pause -1

set out 'arcs2.ps'
plot  [-3:3][-3:3]  'arcdata/arcs2.txt' u 1:2
pause -1

set out 'arcs3.ps'
plot  [-3:3][-3:3]  'arcdata/arcs3.txt' u 1:2
pause -1

set out 'arcs4.ps'
plot  [-3:3][-3:3]  'arcdata/arcs4.txt' u 1:2
pause -1

set out 'arcs5.ps'
plot  [-3:3][-3:3]  'arcdata/arcs5.txt' u 1:2
pause -1

set out 'arcs6.ps'
plot  [-3:3][-3:3]  'arcdata/arcs6.txt' u 1:2
pause -1

set out 'arcs7.ps'
plot  [-3:3][-3:3]  'arcdata/arcs7.txt' u 1:2
pause -1

set out 'arcs8.ps'
plot  [-3:3][-3:3]  'arcdata/arcs8.txt' u 1:2
pause -1

set out 'arcs9.ps'
plot  [-3:3][-3:3]  'arcdata/arcs9.txt' u 1:2
pause -1

set out 'arcs10.ps'
plot  [-3:3][-3:3]  'arcdata/arcs10.txt' u 1:2
pause -1

set out 'arcs11.ps'
plot  [-3:3][-3:3]  'arcdata/arcs11.txt' u 1:2
pause -1

set out 'arcs12.ps'
plot  [-3:3][-3:3]  'arcdata/arcs12.txt' u 1:2
pause -1

set out 'arcs13.ps'
plot  [-3:3][-3:3]  'arcdata/arcs13.txt' u 1:2
pause -1

set out 'arcs14.ps'
plot  [-3:3][-3:3]  'arcdata/arcs14.txt' u 1:2
pause -1

set out 'arcs15.ps'
plot  [-3:3][-3:3]  'arcdata/arcs15.txt' u 1:2
pause -1

set out 'arcs16.ps'
plot  [-3:3][-3:3]  'arcdata/arcs16.txt' u 1:2
pause -1

set out 'arcs17.ps'
plot  [-3:3][-3:3]  'arcdata/arcs17.txt' u 1:2
pause -1

set out 'arcs18.ps'
plot  [-3:3][-3:3]  'arcdata/arcs18.txt' u 1:2
pause -1

set out 'arcs19.ps'
plot  [-3:3][-3:3]  'arcdata/arcs19.txt' u 1:2
pause -1

set out 'arcs20.ps'
plot  [-3:3][-3:3]  'arcdata/arcs20.txt' u 1:2
pause -1


set terminal X11