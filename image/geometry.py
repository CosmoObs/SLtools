"""
The following functions execute some geometrical operations such as calculate the perpedicular bisector between two points.
"""

from math import sin, cos ,sqrt, fabs, atan,tan 
#from pylab import *

def Threepoints_circle(x,y):

    """
    Function that calculates circle parameters defined by three points.
    
    Input:
     - x          <list float> : three x coordinates from the points
     - y          <list float> : three y coordinates from the points
        
    Output:
     - <float> : center x coordinate
     - <float> : center y coordinate
     - <float> : circle radius
    """
    p1=[x[0],y[0]]
    p2=[x[1],y[1]]
    p3=[x[2],y[2]]
    theta1,b1=Perpendicular_Bisector(p1,p2)
    theta2,b2=Perpendicular_Bisector(p2,p3)
    xc,yc,Flag=Intercept_point(theta1,b1,theta2,b2)
    if Flag==1:
        return 0
    r=((p1[0]-xc)*(p1[0]-xc))+((p1[1]-yc)*(p1[1]-yc))
    r=sqrt(r)
    return xc,yc,r


def Intercept_point(theta,b_t,phi,b_p):

    """
    Function that calculates the interception point from two straight lines.
    
    Input:
     - theta      <float> : angle defined by the first line and the x axis
     - b_t        <float> : linear coeficient from the first line in the line equation y=ax+b
     - phi        <float> : angle defined by the second line and the x axis
     - b_p        <float> : linear coeficient from the second line in the line equation y=ax+b
   
    Output:
     - <float> : intercept x coordinate
     - <float> : intercept y coordinate
     - <int>   : a Flag. if Flag= 1 the lines are parallels, otherwise Flag=0
    """

    a_t=tan(theta)
    a_p=tan(phi)
    Flag=0
    if(fabs(a_t-a_p)<0.001):
        print "parallels found"
        Flag=1
        xi=0
        yi=0
    elif (theta==2*atan(1) and phi==0.):
        xi=b_t
        yi=b_p
    elif (theta==0. and phi==2*atan(1)):
        xi=b_p
        yi=b_t
    elif theta==2*atan(1):
        xi=b_t
        yi=a_p*xi+b_p
    elif phi==2*atan(1):
        xi=b_p
        yi=a_t*xi+b_t
    elif theta==0:
        yi=b_t
        xi=(yi-b_p)/a_p
    elif phi==0:   
        yi=b_p
        xi=(yi-b_t)/a_t
    else:
        xi=(b_p-b_t)/(a_t-a_p)
        yi=a_t*xi+b_t
    return xi,yi,Flag


def Perpendicular_Bisector(p1,p2): 
    """
    Function that calculates the perpendicular bisector line between two points.
    
    Input:
     - p1      <float> : angle defined by the first line and the x axis
     - p2      <float> : linear coeficient from the first line in the line equation y=ax+b
        
    Output:
     - <float> : angle defined by the perpendicular bisector line and the x axis
     - <float> : perpendicular bisector linear coeficient in the line equation y=ax+b
     
    """

	x1=p1[0]
	y1=p1[1]	
	x2=p2[0]
	y2=p2[1]	
	Dy=float(y2-y1)
	Dx=float(x2-x1)
	b1=0
	My=0
	Mx=0
	if (fabs(Dy)>0): 
		m=-(Dx)/(Dy)
		theta=atan(m)
		My=float(y2+y1)/2.
		Mx=float(x2+x1)/2.
		
		b1=My-m*Mx			
	elif Dy==0:
		m=0
		theta=2*atan(1)
		b1=float(x2+x1)/2.		
			
	return theta,b1	
