from geometry import *
from get_extrema_2loops import *
import os
import numpy
from math import sin, cos ,sqrt, fabs, atan 
from pylab import *


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

def Mediatrix_Decomposition(image, method="medium",alpha=1,nearDistance=(sqrt(2)/2)):
   """
    Function to perform the mediatrix decomposition method on a given object. 
    
    Input:
     - image   <array> : the image with a single object
     - method   <string> : 'medium' if you want to pick as Mediatrix Point the point at same distance from the two most distant points near the bisector. Or 'brightest' if you want to pick as Mediatrix Point the brightest point near the bisector. DEFAULT 'medium'
     - alpha      <float> : the factor alpha=l_i/w to stop the bisection, w is the estimated width and l_i the distance from two neighbours keydots. if l_i<alpha*w the function will not bisect to next Mediatrix level. DEFAULT 1
     - nearDistance      <list> : the hightest distance to consider a point near to the perpendicular bisector. DEFAULT sqrt(2)/2
     
    Output:
     - <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'modulus' for vector modulus
         
    """
	pixels=numpy.where(image>0)
	E1,E2=get_extrema_2loops( pixels[0], pixels[1], 0 )
	Area=len(pixels[1])
	p1=[pixels[0][E1],pixels[1][E1]] # the extreme points p_1 and p_2
	p2=[pixels[0][E2],pixels[1][E2]]
	keydots=[p1,p2]
	keydots=FindKeydots(p1,p2,pixels,image,keydots,Area, method="brightest",alpha=1,nearDistance)
	mediatrix_vectors=Find_Mediatrix_vectors(keydots)
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
	keypoints=FindKeydots(p1,p2,pixels,image,keypoints,Area, method="brightest",alpha=1,nearDistance)
	L_f=lenght(keypoints)
	
	return L_f

def Find_Mediatrix_vectors(keydots): 
"""
    From a given set of points, this function returns the mediatrix decomposition vector between those points.
    
    Input:
     - p1      <list> : list  with the keydots points
     - p2      <list> : coordinates (x,y) of the second extreme point
     - image_pixels   <list> : list of points coordidates fitting the object
     - image   <array> : the image with a single object
     - keydots  <list> : list with the two extreme points e.g. [p_1,p_2]. And p_i=[p_i_x,p_i_y]
     - Area  <list> : number of points belonging to the object
     - method   <string> : 'medium' if you want to pick as Mediatrix Point the point at same distance from the two most distant points near the bisector. Or 'brightest' if you want to pick as Mediatrix Point the brightest point near the bisector. DEFAULT 'medium'
     - alpha      <float> : the factor alpha=l_i/w to stop the bisection, w is the estimated width and l_i the distance from two neighbours keydots. if l_i<alpha*w the function will not bisect to next Mediatrix level. DEFAULT 1
     - nearDistance      <list> : the hightest distance to consider a point near to the perpendicular bisector. DEFAULT sqrt(2)/2
     
    Output:
     - <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'modulus' for vector modulus  
     
    """
	vectors=[]
	t=0
	for t in range(0,len(p)-1):
		coefficients=Perpendicular_Bisector(p[t],p[t+1])
		origin_x=float(p[t][0]+p[t+1][0])/2
		origin_y=float(p[t][1]+p[t+1][1])/2
		origin=[origin_x,origin_y]
		modulus=(p[t][0] - p[t+1][0])**2 + (p[t][1] - p[t+1][1])**2
		modulus=sqrt(modulus)
		mediatrix_vector = {'theta': coefficients[0], 'linear_coefficient': coefficients[1], 'origin': origin, 'modulus': modulus }
		vectors.append(mediatrix_vector)
				
		
	return vectors



def ChooseNearPoint(theta,c,pixels,image,method,maxDistance): 
	"""
    This function choose the Mediatrix points from a perpendicular bisector parameters in an object.

    Input:
     - theta          <float> : angle defined by the line and the x axis
     - c              <float> : linear coeficient from the line in the line equation y=ax+c
     - image_pixels   <list> : list of points coordidates fitting the object
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
	if(theta!=2*atan(1)) and (theta!=0.):
		a=tan(theta)
		m=-1./a
		for i in range(0,len(pixels[1])):
			b=pixels[1][i]-m*pixels[0][i]
			x=(c-b)/(m-a)
			y=a*x+c
			dx=pixels[0][i]-x
			dy=pixels[1][i]-y
			D=sqrt(dx*dx + dy*dy)
			if(maxDistance>D):
				nearX.append(pixels[0][i])
				nearY.append(pixels[1][i])
				
	elif (theta==2*atan(1)): # case of the perpendicular bisector is a vertical line a-->inf.
		for i in range(0,len(pixels[1])):
			dx=fabs(pixels[0][i]-c)
			D=dx
			if(maxDistance>D):
				nearX.append(pixels[0][i])
				nearY.append(pixels[1][i])
	elif theta==0.: # case of the perpendicular bisector is a horizontal line a-->0
		for i in range(0,len(pixels[1])):
			dy=fabs(pixels[1][i]-c)
			D=dy
			if(maxDistance>D):
				nearX.append(pixels[0][i])
				nearY.append(pixels[1][i])
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
	#X=points[0]
	#Y=points[1]	
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

