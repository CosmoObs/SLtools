"""
The following functions execute elementary geometrical calculations.
"""

from math import sin, cos ,sqrt, fabs, atan, tan 


def three_points_to_circle(p1, p2, p3):

    """
    Function that calculates circle parameters defined by three points.
    
    Input:
	 - p1 [float,float] : pair (x,y) of coordinates of the 1st point
	 - p2 [float,float] : pair (x,y) of coordinates of the 2nd point
	 - p3 [float,float] : pair (x,y) of coordinates of the 3rd point
        
    Output:
     - float : center x coordinate
     - float : center y coordinate
     - float : circle radius
    """

    theta1,b1=define_perpendicular_bisector(p1,p2)
    theta2,b2=define_perpendicular_bisector(p2,p3)
    xc,yc,Flag=get_intercept_point(theta1,b1,theta2,b2)
    if Flag==1:
        return 0,0,0
    r=((p1[0]-xc)*(p1[0]-xc))+((p1[1]-yc)*(p1[1]-yc))
    r=sqrt(r)
    return xc,yc,r


def three_points_to_circle_2(p1, p2, p3):
	"""
	Determines the circunference that passes through 3 points.

	Given 3 points (x,y) on a plane, returns circunference radius and center. 
        The method implemented
	in this funcion is the one described in 
        <A HREF="http://en.wikipedia.org/wiki/Circumcircle"> 
	CircumcircleWikipedia</A> .

	Input:
	 - p1 [float,float] : pair (x,y) of coordinates of the 1st point
	 - p2 [float,float] : pair (x,y) of coordinates of the 2nd point
	 - p3 [float,float] : pair (x,y) of coordinates of the 3rd point

	Output: 
	 - float : x position of the circunference center
	 - float : y position of the circunference center
         - float : circunference radius


	"""

	x1, y1 = p1
	x2, y2 = p2
	x3, y3 = p3

	# definitions
	P1_P2 = ( (x1 - x2)**2 + (y1 - y2)**2 )**0.5
	P2_P3 = ( (x2 - x3)**2 + (y2 - y3)**2 )**0.5
	P3_P1 = ( (x3 - x1)**2 + (y3 - y1)**2 )**0.5
	P1_P2_vec_P2_P3 = abs( (x1 - x2)*(y2 - y3) - (x2 - x3)*(y1 - y2) )
	# ==================================================================
	# radius

	R = (P1_P2 * P2_P3 * P3_P1)/( 2*P1_P2_vec_P2_P3 )

	# definitions
	alpha = ( P2_P3**2 * ( ( (x1 - x2)* (x1 - x3) ) + ( (y1 - y2)* (y1 - y3) ) ) )/ \
(2*P1_P2_vec_P2_P3**2)
	beta  = ( P3_P1**2 * ( ( (x2 - x1)* (x2 - x3) ) + ( (y2 - y1)* (y2 - y3) ) ) )/ \
(2*P1_P2_vec_P2_P3**2)
	gamma = ( P1_P2**2 * ( ( (x3 - x1)* (x3 - x2) ) + ( (y3 - y1)* (y3 - y2) ) ) )/ \
(2*P1_P2_vec_P2_P3**2)
	# ==================================================================
	# center (x0, y0)

	x0 = alpha*x1 + beta*x2 + gamma*x3
	y0 = alpha*y1 + beta*y2 + gamma*y3


	return x0, y0, R


def get_intercept_point(theta,b_t,phi,b_p):

    """
    Function that calculates the interception point from two straight lines.
    
    Input:
     - theta      float : angle defined by angular coefficient from first line.
     - b_t        float : linear coeficient from the first line.
     - phi        float : angle defined by angular coefficient from second line.
     - b_p        float : linear coeficient from the second line.
   
    Output:
     - float : intercept x coordinate.
     - float : intercept y coordinate.
     - int   : a Flag. if Flag= 1 the lines are parallels, otherwise Flag=0,
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
    Function that calculates the perpendicular bisector line between two points 
    p_i=[x_i,y_i].
    
    Input:
     - p1      [float,float] : first point.
     - p2      [float,float] : second point.
        
    Output:
    - float : angle defined by angular coefficient.
    - float : linear coefficient.
     
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
     - p1      [float,float] : first point.
     - p2      [float,float] : second point.
        
    Output:
     - float : angle defined by angular coefficient.
     - float : linear coefficient.
     
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
      - p       [float,float]  : a point p=[x,y]. 	
      - theta          float : angle defined by line angular coefficient.
      - c              float : linear coefficient.
    Output:
     - float : The distance from the line to the point p.

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

def width_ellipse(l,area):
    """
    Function to calculate wiidth of an ellipse of given area and length l.

    Input:
     - l          float : the ellipse length.
     - area       float : the ellipse area.
     
     
    Output: 
     - float : the width.  
     
    """


	w=area/(atan(1)*l)
	return w

def length_from_connected_dots(points): 
    """
    Function to calculate the sum of distance from a previous point to the next on 
    a list of points.

    Input:
     - points  [[float,float],...] : list  of p_i points where p_i=(x_i,y_i).
     
     
    Output: 
     - float : the length.  
     
    """
    x=[]
    y=[]
    for i in range(0,len(points)):
        x.append(points[i][0])
        y.append(points[i][1])
    l=0
    err=0
    for j in range(0,len(X)-1):
        dx=x[j]-x[j+1]
        dy=x[j]-y[j+1]
        dx=dx*dx
        dy=dy*dy
        dl=sqrt(dx+dy)
        l=l+dl

    return float(l)	



