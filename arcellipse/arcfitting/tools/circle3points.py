from math import sin, cos ,sqrt, fabs, atan,tan 
#from pylab import *

def Threepoints_circle(x,y):
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
#  given two points it returns the perpendicular bissector between them y=mx+b1
# if the line is vertical it retuns a Flag=1 and the line equation x=b1
# if the line is horizontal it returns Flag=2 and the line equation y=b1
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
		#print " Coeficientes "+str(b1)+":"+str(m)
		#print "Mx e My "+str(Mx)+":"+str(My)
			
	elif Dy==0:
		m=0
		theta=2*atan(1)
		b1=float(x2+x1)/2. 
		
			
	return theta,b1	
