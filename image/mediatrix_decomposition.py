from geometry import *
from get_extrema_2loops import *
import os
import numpy
from math import sin, cos ,sqrt, fabs, atan 
from pylab import *
import pyfits


def FindKeydots (p1,p2,image_pixels,image,keydots,Area, method="medium",alpha=1,nearDistance=(sqrt(2)/2)):
    """
    Function to calculate the keydots Points in Mediatrix Decomposition.
    
    Input:
     - p1      <list> : coordinates (x,y) of the first extreme point
     - p2      <list> : coordinates (x,y) of the second extreme point
     - image_pixels   <list> : list of points coordidates fitting the object
     - image   <array> : the image with a single object
     - keydots  <list> : list with the two extreme points e.g. [p_1,p_2]. And p_i=[p_i_x,p_i_y]
     - Area  <list> : number of points belonging to the object
     - method   <string> : 'medium' if you want to pick as Mediatrix Point the point at same distance from the two most distant points near the bisector. Or 'brightest' if you want to pick as Mediatrix Point the brightest point near the bisector. DEFAULT 'medium'
     - alpha      <float> : the factor alpha=l_i/w to stop the bisection, w is the estimated width and l_i the distance from two neighbours keydots. if l_i<alpha*w the function will not bisect to next Mediatrix level. DEFAULT 1
     - nearDistance      <list> : the hightest distance to consider a point near to the perpendicular bisector. DEFAULT sqrt(2)/2
     
    Output:
     - <list> : list  with the keydots points  
     
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
    L=lenght(keydots)
    W=width_ellipse(L,Area)

    if dl>(alpha*W):
        coefficients=Perpendicular_Bisector(p1,p2)
        p3x,p3y,p3Flag=ChooseNearPoint(coefficients[0],coefficients[1],image_pixels,image,method,nearDistance)
        p3=[p3x,p3y]
        if (p3Flag==0):		
            if (not(p3 in keydots)):
                keydots.insert(indexNext,p3)
                keydots=FindKeydots(p1,p3,image_pixels,image,keydots,Area, method,alpha,nearDistance)
                keydots=FindKeydots(p3,p2,image_pixels,image,keydots,Area, method,alpha,nearDistance)
        else:
            xmed=float(x1+x2)/2.
            ymed=float(y1+y2)/2.
            pmed=[xmed,ymed]
            if p1 in keydots: 
                keydots=FindKeydots(p1,pmedimage_pixels,image,keydots,Area, method,alpha,nearDistance)
            if p2 in keydots:
                keydots=FindKeydots(pmed,p2,image_pixels,image,keydots,Area, method,alpha,nearDistance)

    return keydots

def Mediatrix_Decomposition(image_name,PATH='', method="medium",alpha=1,nearDistance=(sqrt(2)/2)):
    """
    Function to perform the mediatrix decomposition method on a given object. 

    Input:
     - image_name   <str> : the image name. This image is a POST STAMP with a single object
     - PATH   <str> : the image Path. If it is on the same directory PATH=''
     - method   <string> : 'medium' if you want to pick as Mediatrix Point the point at same distance from the two most distant points near the bisector. Or 'brightest' if you want to pick as Mediatrix Point the brightest point near the bisector. DEFAULT 'medium'
     - alpha      <float> : the factor alpha=l_i/w to stop the bisection, w is the estimated width and l_i the distance from two neighbours keydots. if l_i<alpha*w the function will not bisect to next Mediatrix level. DEFAULT 1
     - nearDistance      <list> : the hightest distance to consider a point near to the perpendicular bisector. DEFAULT sqrt(2)/2
     
    Output:
     - <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image_name.
         
    """
    image,hdr = pyfits.getdata(PATH+image_name, header = True )
    pixels=numpy.where(image>0)
    E1,E2=get_extrema_2loops( pixels[0], pixels[1], 0 )
    Area=len(pixels[1])
    p1=[pixels[0][E1],pixels[1][E1]] # the extreme points p_1 and p_2
    p2=[pixels[0][E2],pixels[1][E2]]
    keydots=[p1,p2]
    keydots=FindKeydots(p1,p2,pixels,image,keydots,Area, method="brightest",alpha=1,nearDistance=nearDistance)
    mediatrix_vectors=Find_Mediatrix_vectors(keydots)
    mediatrix_vectors[0]['id']=image_name
    return mediatrix_vectors


def LenghtMediatrix(image, method="medium",alpha=1,nearDistance=(sqrt(2)/2)):
    """
    Function to calculate lenght of an object using the Mediatrix Decomposition.

    Input:
     - image   <array> : the image with a single object
     - method   <string> : 'medium' if you want to pick as Mediatrix Point the point at same distance from the two most distant points near the bisector. Or 'brightest' if you want to pick as Mediatrix Point the brightest point near the bisector. DEFAULT 'medium'
     - alpha      <float> : the factor alpha=l_i/w to stop the bisection, w is the estimated width and l_i the distance from two neighbours keydots. if l_i<alpha*w the function will not bisect to next Mediatrix level. DEFAULT 1
     - nearDistance      <list> : the hightest distance to consider a point near to the perpendicular bisector. DEFAULT sqrt(2)/2
     
    Output:
     - <Float> : the object lenght  
     
    """
    pixels=numpy.where(image>0)
    E1,E2=get_extrema_2loops( pixels[0], pixels[1], 0 )
    Area=len(pixels[1])
    p1=[pixels[0][E1],pixels[1][E1]]
    p2=[pixels[0][E2],pixels[1][E2]]
    keypoints=[p1,p2]
    keypoints=FindKeydots(p1,p2,pixels,image,keypoints,Area, method="brightest",alpha=1,nearDistance=nearDistance)
    L_f=lenght(keypoints)

    return L_f


def Find_Mediatrix_vectors(p): 
    """
    From a given set of points, this function returns the mediatrix decomposition vector between those points.
  
    Input:
     - p      <list> : list  with the keydots points
     
    Output:
     - <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus.  
     
    """
    vectors=[]
    t=0
    theta_ext,c_ext=Twopoints_line(p[t],p[len(p)-1])
    for t in range(0,len(p)-1):
        coefficients=Perpendicular_Bisector(p[t],p[t+1])
        origin_x=float(p[t][0]+p[t+1][0])/2
        origin_y=float(p[t][1]+p[t+1][1])/2
        origin=[origin_x,origin_y]
        modulus=(p[t][0] - p[t+1][0])**2 + (p[t][1] - p[t+1][1])**2
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
        Dend1=Distance_from_line_to_point(end1,theta_ext,c_ext)
        Dend2=Distance_from_line_to_point(end2,theta_ext,c_ext)
        if Dend1<Dend2:
	    end=end1
        else:
            end=end2

        mediatrix_vector = {'theta': coefficients[0], 'linear_coefficient': coefficients[1], 'origin': origin, 'end': end, 'modulus': modulus }
        vectors.append(mediatrix_vector)
                
        
    return vectors

def print_mediatrix_Object_graph(mediatrix_data,image_Path, keydots=False, colors= {'object': "m.", 'vector': "g", 'keydots': "kD"}):
    """
    Make a plot presenting the object, keydots and mediatrix vectors. 

    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.
    - image_Path   <str> : the image Path. If it is on the same directory PATH=''
    - keydots   <bool> : 'True' if you want to display the keydots and 'False' if you do not. 
    - colors   <dic> : set the plot colors and dot style. The possible keys are 'object', 'vector' and 'keydots'       
    Output:
     none
         
    """
    
    image_name=mediatrix_data[0]['id']
    image,hdr = pyfits.getdata(image_Path+image_name, header = True )
    pixels=numpy.where(image>0)    
    A = subplot(111)
    A.plot(pixels[1],pixels[0],colors['object'])
    Lenght=0
    
      
    if keydots==True:
        E1,E2=get_extrema_2loops(pixels[0], pixels[1], 0 )
        Area=len(pixels[1])
        p1=[pixels[0][E1],pixels[1][E1]] # the extreme points p_1 and p_2
        p2=[pixels[0][E2],pixels[1][E2]]
        keydots=[p1,p2]
        keydots=FindKeydots(p1,p2,pixels,image,keydots,Area, method="brightest",alpha=1)
        keyX=[]
        keyY=[]
        for j in range(0,len(keydots)):
            keyX.append(keydots[j][0])
            keyY.append(keydots[j][1])
        
        A.plot(keyY,keyX,colors['keydots'])
    
    for i in range(0,len(mediatrix_data)):
        origin1=mediatrix_data[i]['origin'][0]
        origin2=mediatrix_data[i]['origin'][1]
        Lenght=Lenght+mediatrix_data[i]['modulus']
        d1= mediatrix_data[i]['end'][0] - mediatrix_data[i]['origin'][0]
        d2= mediatrix_data[i]['end'][1] - mediatrix_data[i]['origin'][1]
        arr = Arrow(origin2, origin1, d2, d1, width=0.05*Lenght, fc=colors['vector'], ec='none',zorder=1000)
        A.add_patch(arr)
    print (1.3/Lenght)
    xmin, xmax = xlim()
    ymin, ymax = ylim()
    A.set_xlim(xmin-1*Lenght,xmax+1*Lenght)
    A.set_ylim(ymin-1*Lenght,ymax+1*Lenght)    
    ylabel("Y")
    xlabel("X")
    #A.axis("equal")
    title("Mediatrix Decomposition applied") 
    #A.axis("equal")
    savefig(image_Path+image_name+"_mediatrixGraph.png")
    A.clear()			
    return True


def ChooseNearPoint(theta,c,object_pixels,image,method,maxDistance): 
    """
    This function choose the Mediatrix points from a perpendicular bisector parameters in an object.

    Input:
     - theta          <float> : angle defined by the line and the x axis
     - c              <float> : linear coeficient from the line in the line equation y=ax+c
     - object_pixels   <list> : list of points coordidates fitting the object
     - image          <array> : the image with a single object
     - method         <string> : 'medium' if you want to pick as Mediatrix Point the point at same distance from the two most distant points near the bisector. Or 'brightest' if you want to pick as Mediatrix Point the brightest point near the bisector. DEFAULT 'medium'
     - maxDistance    <list> : the hightest distance to consider a point near to the perpendicular bisector. DEFAULT sqrt(2)/2
     
    Output:
     - <float> : the chosen point x coordinate  
     - <float> : the chosen point y coordinate
     - <int>   : a Flag error. if Flag=1, it was not possible to choose a point.
    """

    nearX= []
    nearY= []

    for i in range(0,len(object_pixels[1])):
        pixel=[object_pixels[0][i],object_pixels[1][i]]
	D=Distance_from_line_to_point(pixel,theta,c)
        if(maxDistance>=D):
            nearX.append(object_pixels[0][i])
            nearY.append(object_pixels[1][i])
           
    if (len(nearX)<1):
        FlagErr=1
        chosenX=0
        chosenY=0
        return chosenX, chosenY, FlagErr
    if method=='brightest':
        chosenX=nearX.pop()
        chosenY=nearY.pop()
        FlagErr=0
        while len(nearX)>=1:
            chosenAuxX=nearX.pop()
            chosenAuxY=nearY.pop()
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




def lenght(points): 
    """
    Function to calculate the lenght in a path determined by several points.

    Input:
     - points      <list> : coordinates p_i(x_i,y_i) of the points. e.g. points=[p_1,p_2,p_3...,p_n]
     
     
    Output: 
     - <float> : the lenght  
     
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

