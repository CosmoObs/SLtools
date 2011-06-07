"""
Set of functions used in Mediatrix Decomposition.

"""



##@package mediatrix_decomposition 
#
# This module includes all functions developed to perform mediatrix decomposition method.
#
#



from get_extrema_2loops import get_extrema_2loops
from numpy import where
from math import sin, cos ,sqrt, fabs, atan, tan 
from pylab import subplot, Rectangle, Arrow, xlim, ylim, ylabel, xlabel, title, savefig, Circle
from pyfits import getdata
from geometry import define_perpendicular_bisector, get_distance_from_line_to_point, two_points_to_line, three_points_to_circle

def find_keydots (p1,p2,image_pixels,image,keydots,area, method="medium",alpha=1,near_distance=(sqrt(2)/2)):
    """
    Function to calculate the keydots Points in Mediatrix Decomposition.
    
    Input:
     - p1      <array> : coordinates (x,y) of the first extreme point.
     - p2      <array> : coordinates (x,y) of the second extreme point.
     - image_pixels   <list> : list of points coordidates fitting the object.
     - image   <array> : the image matrix.
     - keydots  <array> : list with the two p_i extreme points and p_i=[p_i_x,p_i_y].
     - area  <array> : the object area.
     - method   <string> : possible values are 'medium' or 'brightest'.
     - alpha      <float> : the factor alpha=l_i/w.
     - near_distance      <float> : the distance to consider a point near to the perpendicular bisector.
     
    Output:
     - <array> : list  with the keydots points.  
     
    """
    if (p1 in keydots) and (p2 in keydots):
        index1=keydots.index(p1)
        index2=keydots.index(p2)
        if  index1>index2:
            indexNext=index1
        else:
            indexNext=index2
    elif (p1 in keydots):
        indexNext=keydots.index(p1)
    elif (p1 in keydots):
        indexNext=keydots.index(p2)
    else:
        return keydots
	
    x1=p1[0]
    x2=p2[0]
    y1=p1[1]
    y2=p2[1]
    dx=(x1-x2)
    dy=(y1-y2)
    dl=sqrt((dx*dx)+(dy*dy))
    L=get_length(keydots)
    W=width_ellipse(L,area)

    if dl>(alpha*W) and len(keydots)<100:
        coefficients=define_perpendicular_bisector(p1,p2)
        p3x,p3y,p3Flag=choose_near_point(coefficients[0],coefficients[1],image_pixels,image,method,near_distance)
        p3=[p3x,p3y]
        if (p3Flag==0):		
            if (not(p3 in keydots)):
                keydots.insert(indexNext,p3)
                keydots=find_keydots(p1,p3,image_pixels,image,keydots,area, method,alpha,near_distance)
                keydots=find_keydots(p3,p2,image_pixels,image,keydots,area, method,alpha,near_distance)
        else:
            xmed=float(x1+x2)/2.
            ymed=float(y1+y2)/2.
            pmed=[xmed,ymed]
            if p1 in keydots: 
                keydots=find_keydots(p1,pmed,image_pixels,image,keydots,area, method,alpha,near_distance)
            if p2 in keydots:
                keydots=find_keydots(pmed,p2,image_pixels,image,keydots,area, method,alpha,near_distance)

    return keydots

def mediatrix_decomposition(image_name,image_dir='', method="medium",alpha=1,near_distance=(sqrt(2)/2)):
    """
    Function to perform the mediatrix decomposition method on a given object. 

    Input:
     - image_name   <str> : the image file name.
     - image_dir   <str> : the image directory. If it is on the same directory, directory=''.
     - method   <string> : possible values are 'medium'  or 'brightest'.
     - alpha      <float> : the factor alpha=l_i/w to stop the bisection.
     - near_distance      <float> : the distance to consider a point near to the perpendicular bisector.
     
    Output:
     - <dic> :  Dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has two extra keys 'id' wich contains the image_name and 'center' that keeps the objet center defined by the first mediatrix point in the first mediatrix level.
         
    """
    image,hdr = getdata(image_dir+image_name, header = True )
    pixels=where(image>0)
    E1,E2=get_extrema_2loops( pixels[0], pixels[1], 0 )
    Area=len(pixels[1])
    p1=[pixels[0][E1],pixels[1][E1]] # the extreme points p_1 and p_2
    p2=[pixels[0][E2],pixels[1][E2]]
    keydots=[p1,p2]
    keydots=find_keydots(p1,p2,pixels,image,keydots,Area, method=method,alpha=alpha,near_distance=near_distance)
    mediatrix_vectors=find_mediatrix_vectors(keydots)
    mediatrix_vectors['id']=image_name
    medium=int(float(len(keydots))/2)
    mediatrix_vectors['center']=keydots[medium]
    L=get_length(keydots)
    W=len(pixels[0])/(atan(1)*L)
    mediatrix_vectors['L/W']=L/W
    mediatrix_vectors['L']=L
    x=[pixels[0][E1],mediatrix_vectors['center'][0],pixels[0][E2]]
    y=[pixels[1][E1],mediatrix_vectors['center'][1],pixels[1][E2]]
    x_c,y_c,r=three_points_to_circle(x,y)
    circle_center=[x_c,y_c]
    mediatrix_vectors['circle_params']=[circle_center,p1,p2]

    return mediatrix_vectors


def get_length_from_mediatrix(image_name,image_dir='', method="medium",alpha=1,near_distance=(sqrt(2)/2)):
    """
    Function to calculate length of an object using the Mediatrix Decomposition.

    Input:
     - image_name   <array> : the image file name.
     - method   <string> : possible values are 'medium'  or 'brightest'.
     - alpha      <float> : the factor alpha=l_i/w.
     - near_distance      <float> : the distance to consider a point near to the perpendicular bisector.
     
    Output:
     - <Float> : the object length  
     
    """
    image,hdr = getdata(image_dir+image_name, header = True )
    pixels=where(image>0)
    E1,E2=get_extrema_2loops( pixels[0], pixels[1], 0 )
    Area=len(pixels[1])
    p1=[pixels[0][E1],pixels[1][E1]]
    p2=[pixels[0][E2],pixels[1][E2]]
    keypoints=[p1,p2]
    keypoints=find_keydots(p1,p2,pixels,image,keypoints,Area, method="brightest",alpha=alpha,near_distance=near_distance)
    L_f=get_length(keypoints)
    W=width_ellipse(L_f,Area)
    
    return L_f, W


def find_mediatrix_vectors(points): 
    """
    From a given set of points, this function returns the mediatrix decomposition vector between those points.
  
    Input:
     - points      <list> : list  of p_i points where p_i=(x_i,y_i).
     
    Output:
     - <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus.  
     
    """
    mediatrix_vectors= {'origin': [] , 'end': [], }
    vectors=[]
    t=0
    theta_ext,c_ext=two_points_to_line(points[t],points[len(points)-1])
    for t in range(0,len(points)-1):
        coefficients=define_perpendicular_bisector(points[t],points[t+1])
        origin_x=float(points[t][0]+points[t+1][0])/2
        origin_y=float(points[t][1]+points[t+1][1])/2
        origin=[origin_x,origin_y]
        modulus=(points[t][0] - points[t+1][0])**2 + (points[t][1] - points[t+1][1])**2
        modulus=sqrt(modulus)
	if(coefficients[0]!=2*atan(1)):
	    a=tan(coefficients[0])
            sq_a=a*a
            x_end1=origin_x+(modulus/sqrt(1.+sq_a))
            x_end2=origin_x-(modulus/sqrt(1.+sq_a))
            y_end1=a*x_end1+coefficients[1]
            y_end2=a*x_end2+coefficients[1]
            
        else:
            y_end1=origin_y+modulus
            y_end2=origin_y-modulus
            x_end1=origin_x
            x_end2=origin_x
        end1=[x_end1,y_end1]
        end2=[x_end2,y_end2]
        Dend1=get_distance_from_line_to_point(end1,theta_ext,c_ext)
        Dend2=get_distance_from_line_to_point(end2,theta_ext,c_ext)
        if Dend1<Dend2:
	    end=end1
        else:
            end=end2

        #mediatrix_vector = {'theta': coefficients[0], 'linear_coefficient': coefficients[1], 'origin': origin, 'end': end, 'modulus': modulus }
        mediatrix_vectors['origin'].append(origin) 
        mediatrix_vectors['end'].append(end)
        #vectors.append(mediatrix_vector)
                
        
    return mediatrix_vectors

def print_mediatrix_Object_graph(mediatrix_data,image_dir='', keydots=False, colors= {'object': "g", 'vector': "b", 'keydots': "k"}, save=True, save_dir=''):
    """
    Make a plot presenting the object, keydots and mediatrix vectors. 

    Input:
    - mediatrix_data <list> : the output from mediatrix_decomposition.
    - image_dir   <str> : the image directory. If it is on the same directory, directory=''.
    - keydots   <bool> : 'True' if you want to display the keydots and 'False' if you do not. 
    - colors   <dic> : set the plot colors. The possible keys are 'object', 'vector' and 'keydots'.       
    Output:
     <bool>
         
    """
    
    image_name=mediatrix_data['id']
    image,hdr = getdata(image_dir+image_name, header = True )
    pixels=where(image>0)    
    A = subplot(111)
    for i in range (0,len(pixels[0])):
        xy=[pixels[1][i]-0.5,pixels[0][i]-0.5]
        rec=Rectangle(xy, 1, 1, ec=colors['object'], fc=colors['object'], zorder=100)
        A.add_patch(rec)
    #A.scatter(pixels[1], pixels[0], s=200, c='b', marker='s', edgecolors='none')
    #A.plot(pixels[1],pixels[0],colors['object'])
    Length=0
    
      
    if keydots==True:
        E1,E2=get_extrema_2loops(pixels[0], pixels[1], 0 )
        Area=len(pixels[1])
        p1=[pixels[0][E1],pixels[1][E1]] # the extreme points p_1 and p_2
        p2=[pixels[0][E2],pixels[1][E2]]
        keydots=[p1,p2]
        keydots=find_keydots(p1,p2,pixels,image,keydots,Area, method="brightest",alpha=1)
        keyX=[]
        keyY=[]
        for j in range(0,len(keydots)):
            keyX.append(keydots[j][0])
            keyY.append(keydots[j][1])
        
        A.plot(keyY,keyX,colors['keydots']+'.',zorder=500)
        #A.scatter(keyY, keyX, s=20, c='b', marker='s')

    
    for i in range(0,len(mediatrix_data['origin'])):
        origin_x=mediatrix_data['origin'][i][0]
        origin_y=mediatrix_data['origin'][i][1]
        end_x=mediatrix_data['end'][i][0]
        end_y=mediatrix_data['end'][i][1]
        Length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
        Length=Length+ sqrt(Length_aux)
        d_x= end_x - origin_x
        d_y= mediatrix_data['end'][i][1] - mediatrix_data['origin'][i][1]
        arr = Arrow(origin_y, origin_x, d_y, d_x, width=0.05*Length, fc=colors['vector'], ec='none',zorder=1000)
        A.add_patch(arr)
   
    xmin, xmax = xlim()
    ymin, ymax = ylim()
    min_inc_axis=40
    #x_axis_length=(xmax+1*Length)-(xmin-1*Length)
    #y_axis_length=(ymax+1*Length)-(ymin-1*Length)
    #if  x_axis_length<min_inc_axis
    A.axis("equal")
    A.set_xlim(xmin-1*Length,xmax+1*Length)
    A.set_ylim(ymin-1*Length,ymax+1*Length)    
    ylabel("Y")
    xlabel("X")
    #A.axis("equal")
    title("Mediatrix Decomposition applied") 
    
    if save==True:
        savefig(save_dir+image_name+"_mediatrixGraph.png")
        A.clear()
        return True
    else:
        return A



def print_mediatrix_Object_circle_graph(mediatrix_data,image_dir='', keydots=False, colors= {'object': "g", 'vector': "b", 'keydots': "k"}, mediatrix_vectors=False, save=True, save_dir=''):
    """
    Make a plot presenting the object, keydots and mediatrix vectors. 

    Input:
    - mediatrix_data <list> : the output from mediatrix_decomposition.
    - image_dir   <str> : the image directory. If it is on the same directory directory=''.
    - keydots   <bool> : 'True' if you want to display the keydots and 'False' if you do not. 
    - colors   <dic> : set the plot colors. The possible keys are 'object', 'vector' and 'keydots'.       
    Output:
     <bool>
         
    """
    
    image_name=mediatrix_data['id']
    image,hdr = getdata(image_dir+image_name, header = True )
    pixels=where(image>0)    
    A = subplot(111)
    for i in range (0,len(pixels[0])):
        xy=[pixels[1][i]-0.5,pixels[0][i]-0.5]
        rec=Rectangle(xy, 1, 1, ec=colors['object'], fc=colors['object'], zorder=100)
        A.add_patch(rec)
    #A.scatter(pixels[1], pixels[0], s=200, c='b', marker='s', edgecolors='none')
    #A.plot(pixels[1],pixels[0],colors['object'])
    Length=0
    for i in range(0,len(mediatrix_data['origin'])):
        origin_x=mediatrix_data['origin'][i][0]
        origin_y=mediatrix_data['origin'][i][1]
        end_x=mediatrix_data['end'][i][0]
        end_y=mediatrix_data['end'][i][1]
        Length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
        Length=Length+ sqrt(Length_aux)
        if mediatrix_vectors==True:
            d_x= end_x - origin_x
            d_y= mediatrix_data['end'][i][1] - mediatrix_data['origin'][i][1]
            arr = Arrow(origin_y, origin_x, d_y, d_x, width=0.05*Length, fc=colors['vector'], ec='none',zorder=1000)
            A.add_patch(arr)
      
    if keydots==True:
        E1,E2=get_extrema_2loops(pixels[0], pixels[1], 0 )
        Area=len(pixels[1])
        p1=[pixels[0][E1],pixels[1][E1]] # the extreme points p_1 and p_2
        p2=[pixels[0][E2],pixels[1][E2]]
        keydots=[p1,p2]
        keydots=find_keydots(p1,p2,pixels,image,keydots,Area, method="brightest",alpha=1)
        keyX=[]
        keyY=[]
        for j in range(0,len(keydots)):
            keyX.append(keydots[j][0])
            keyY.append(keydots[j][1])
        
        A.plot(keyY,keyX,colors['keydots']+'.',zorder=500)
        #A.scatter(keyY, keyX, s=20, c='b', marker='s')

    
    last=len(mediatrix_data['origin'])-1
    x=[pixels[0][E1],mediatrix_data['center'][0],pixels[0][E2]]
    y=[pixels[1][E1],mediatrix_data['center'][1],pixels[1][E2]]
    x_c,y_c,r=three_points_to_circle(x,y)
    if r>0:
        xy=[y_c,x_c]
        cir=Circle(xy,r,fc='none',ec='m', zorder=501)
        A.add_patch(cir)
    else:
        print "impossible to define a circle in" +str(mediatrix_data['id'])
      

   
    xmin, xmax = xlim()
    ymin, ymax = ylim()
    min_inc_axis=40
    #x_axis_length=(xmax+1*Length)-(xmin-1*Length)
    #y_axis_length=(ymax+1*Length)-(ymin-1*Length)
    #if  x_axis_length<min_inc_axis
    A.axis("equal")
    A.set_xlim(xmin-1*Length,xmax+1*Length)
    A.set_ylim(ymin-1*Length,ymax+1*Length)    
    ylabel("Y")
    xlabel("X")
    #A.axis("equal")
    title("Mediatrix Decomposition applied") 
    
    if save==True and r>0:
        savefig(save_dir+image_name+"_mediatrixGraph_Circle.png")
        A.clear()
        return True
    else:
        return A





def choose_near_point(theta,c,object_pixels,image,method,near_distance): 
    """
    This function choose the Mediatrix points from a perpendicular bisector parameters in an object.

    Input:
     - theta          <float> : straight line angle with x axis.
     - c              <float> : linear coeficient the line.
     - object_pixels   <list> : list of object's points coordidates.
     - image_name   <array> : the image matrix.
     - method   <string> : possible values are 'medium'  or 'brightest'.
     - near_distance      <float> : the distance to consider a point near to the perpendicular bisector.
     
    Output:
     - <float> : the chosen point x coordinate.  
     - <float> : the chosen point y coordinate.
     - <int>   : a Flag error. if Flag=1, it was not possible to choose a point.
    """

    nearXs= []
    nearYs= []

    for i in range(0,len(object_pixels[1])):
        pixel=[object_pixels[0][i],object_pixels[1][i]]
	D=get_distance_from_line_to_point(pixel,theta,c)
        if(near_distance>=D):
            nearXs.append(object_pixels[0][i])
            nearYs.append(object_pixels[1][i])
           
    if (len(nearXs)<1):
        FlagErr=1
        chosenX=0
        chosenY=0
        return chosenX, chosenY, FlagErr
    if method=='brightest':
        chosenX=nearXs.pop()
        chosenY=nearYs.pop()
        FlagErr=0
        while len(nearXs)>=1:
            chosenAuxX=nearXs.pop()
            chosenAuxY=nearYs.pop()
            if (image[chosenX][chosenY])<(image[chosenAuxX][chosenAuxY]):
                chosenX=chosenAuxX
                chosenY=chosenAuxY
            
    elif method=='medium':
        i,j=get_extrema_2loops( nearX, nearY, 0 )
        chosenX=float(nearX[i]+nearX[j])/2.
        chosenY=float(nearY[i]+nearY[j])/2.
        FlagErr=0
    else:
        FlagErr=1
        chosenX=0
        chosenY=0
            
            
    return chosenX,chosenY,FlagErr




def get_length(points): 
    """
    Function to calculate the length in a path determined by several points.

    Input:
     - points      <list> : list  of p_i points where p_i=(x_i,y_i).
     
     
    Output: 
     - <float> : the length.  
     
    """
    X=[]
    Y=[]
    for i in range(0,len(points)):
        X.append(points[i][0])
        Y.append(points[i][1])
    L=0
    err=0
    for j in range(0,len(X)-1):
        dx=X[j]-X[j+1]
        dy=Y[j]-Y[j+1]
        dx=dx*dx
        dy=dy*dy
        dl=sqrt(dx+dy)
        L=L+dl

    return float(L)	

def width_ellipse(L,Area):
	W=Area/(atan(1)*L)
	return W

