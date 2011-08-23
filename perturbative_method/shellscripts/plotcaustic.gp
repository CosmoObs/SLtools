# Start Sept 2010
# Script to plot the variations in caustic
# MSSG

set pointsize 0.5

# set terminal postscript color

set out 'caustic1.ps'
plot [-3:3][-3:3] 'arcdata/caustic1.txt' u 1:2
pause -1

set out 'caustic2.ps'
plot  [-3:3][-3:3]  'arcdata/caustic2.txt' u 1:2
pause -1

set out 'caustic3.ps'
plot  [-3:3][-3:3]  'arcdata/caustic3.txt' u 1:2
pause -1

set out 'caustic4.ps'
plot  [-3:3][-3:3]  'arcdata/caustic4.txt' u 1:2
pause -1

set out 'caustic5.ps'
plot  [-3:3][-3:3]  'arcdata/caustic5.txt' u 1:2
pause -1

set out 'caustic6.ps'
plot  [-3:3][-3:3]  'arcdata/caustic6.txt' u 1:2
pause -1

set out 'caustic7.ps'
plot  [-3:3][-3:3]  'arcdata/caustic7.txt' u 1:2
pause -1

set out 'caustic8.ps'
plot  [-3:3][-3:3]  'arcdata/caustic8.txt' u 1:2
pause -1

set out 'caustic9.ps'
plot  [-3:3][-3:3]  'arcdata/caustic9.txt' u 1:2
pause -1

set out 'caustic10.ps'
plot  [-3:3][-3:3]  'arcdata/caustic10.txt' u 1:2
pause -1

set out 'caustic11.ps'
plot  [-3:3][-3:3]  'arcdata/caustic11.txt' u 1:2
pause -1

set out 'caustic12.ps'
plot  [-3:3][-3:3]  'arcdata/caustic12.txt' u 1:2
pause -1

set out 'caustic13.ps'
plot  [-3:3][-3:3]  'arcdata/caustic13.txt' u 1:2
pause -1

set out 'caustic14.ps'
plot  [-3:3][-3:3]  'arcdata/caustic14.txt' u 1:2
pause -1

set out 'caustic15.ps'
plot  [-3:3][-3:3]  'arcdata/caustic15.txt' u 1:2
pause -1

set out 'caustic16.ps'
plot  [-3:3][-3:3]  'arcdata/caustic16.txt' u 1:2
pause -1

set out 'caustic17.ps'
plot  [-3:3][-3:3]  'arcdata/caustic17.txt' u 1:2
pause -1

set out 'caustic18.ps'
plot  [-3:3][-3:3]  'arcdata/caustic18.txt' u 1:2
pause -1

set out 'caustic19.ps'
plot  [-3:3][-3:3]  'arcdata/caustic19.txt' u 1:2
pause -1

set out 'caustic20.ps'
plot  [-3:3][-3:3]  'arcdata/caustic20.txt' u 1:2
pause -1


set terminal X11