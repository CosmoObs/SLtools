from matplotlib.pyplot import plot, show
from numpy import loadtxt, where, unique1d, array
import time 
x, y, u, v = loadtxt('cc.dat', unpack=True)

delta_curve = 1 # maximum distance between consecutive points in gravlens output (~ gridhi1/ngrid1) definining a curve

start_point=list()
end_point=list()

i=0 # points
while i < (len(x)- 1):
	
    	start_point.append(i)

	while (i < (len(x) - 1)) and ((pow(x[i+1]-x[i],2)+pow(y[i+1]-y[i],2)) < delta_curve): 
		i=i+1

	i=i+1
	end_point.append(i)

print start_point
print end_point

print "# curves ", len(start_point)
# Image Plane

for j in range(len(start_point)):
	plot(x[start_point[j]:end_point[j]],y[start_point[j]:end_point[j]])

show()
# Source Plane 

for j in range(len(start_point)):
	plot(u[start_point[j]:end_point[j]],v[start_point[j]:end_point[j]])

show()


