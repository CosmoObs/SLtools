"""
The following functions execute elementar geometrical calculations.
"""

from math import sin, cos ,sqrt, fabs, atan, tan 


def three_points_to_circle(x,y):

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


def get_intercept_point(theta,b_t,phi,b_p):

    """
    Function that calculates the interception point from two straight lines.
    
    Input:
     - theta      <float> : angle defined by angular coefficient from first line.
     - b_t        <float> : linear coeficient from the first line.
     - phi        <float> : angle defined by angular coefficient from second line.
     - b_p        <float> : linear coeficient from the second line.
   
    Output:
     - <float> : intercept x coordinate.
     - <float> : intercept y coordinate.
     - <int>   : a Flag. if Flag= 1 the lines are parallels, otherwise Flag=0,
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


def define_perpendicular_bisector(p1,p2): 
    """
    Function that calculates the perpendicular bisector line between two points p_i=[x_i,y_i].
    
    Input:
     - p1      <list> : first point.
     - p2      <list> : second point.
        
    Output:
    - <float> : angle defined by angular coefficient.
    - <float> : linear coefficient.
     
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

def two_points_to_line(p1,p2): 
    """
    Function that calculates the line defined by two points  p_i=[x_i,y_i].
    
    Input:
     - p1      <list> : first point.
     - p2      <list> : second point.
        
    Output:
     - <float> : angle defined by angular coefficient.
     - <float> : linear coefficient.
     
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
    if (fabs(Dx)>0): 
        a=(Dy)/(Dx)
        theta=atan(a)
        b1=y1-a*x1			
    elif Dx==0:
        a=0
        theta=2*atan(1)
        b1=x1		
            
    return theta,b1	





def get_distance_from_line_to_point(p,theta,c):
    """
    This Function calculates the distance from a line to a point.

    Input:
      - p              <list>  : a point p=[x,y]. 	
      - theta          <float> : angle defined by line angular coefficient.
      - c              <float> : linear coefficient.
    Output:
     - <float> : The distance from the line to the point p.

    """
    if(theta!=2*atan(1)) and (theta!=0.):
        a=tan(theta)
        m=-1./a
        b=p[1]-m*p[0]
        x=(c-b)/(m-a)
        y=a*x+c
        dx=p[0]-x
        dy=p[1]-y
        Distance=sqrt(dx*dx + dy*dy)
                
    elif (theta==2*atan(1)): # case of the perpendicular bisector is a vertical line a-->inf.
        dx=fabs(p[0]-c)
        Distance=dx

    elif theta==0.: # case of the perpendicular bisector is a horizontal line a-->0
        dy=fabs(p[1]-c)
        Distance=dy
    return Distance 




