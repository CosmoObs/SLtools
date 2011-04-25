from geometry import *
from get_extrema_2loops import *
import os
import numpy
from math import sin, cos ,sqrt, fabs, atan 
from pylab import *
import pyfits
from scipy.optimize import fmin
from mediatrix_decomposition import lenght

def Sum_inner_product(mediatrix_data):

    """
    Function to calculate the sum of scalar products from mediatrix vectors to its 'mirrored' mediatrix vector.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.

     
    Output:
     - <float> : sum of inner product between a mediatrix vector and its neighbour over the number of vectors 
     
    """
    total_inner_prod=0
    for i in range(0,len(mediatrix_data)-1):
        
        component_1=mediatrix_data[i]['end'][0]-mediatrix_data[i]['origin'][0]
        component_2=mediatrix_data[i]['end'][1]-mediatrix_data[i]['origin'][1]
        vector1=[component_1,component_2]

        component_1=mediatrix_data[i+1]['end'][0]-mediatrix_data[i+1]['origin'][0]
        component_2=mediatrix_data[i+1]['end'][1]-mediatrix_data[i+1]['origin'][1]

        vector2=[component_1,component_2]
        
        inner_prod=numpy.inner(vector1,vector2)
        total_inner_prod=total_inner_prod+inner_prod
    total_inner_prod=total_inner_prod/len(mediatrix_data)
    return total_inner_prod

def Sum_oposite_inner_product(mediatrix_data):

    """
    Function to calculate the sum of scalar products from mediatrix vectors to its 'mirrored' mediatrix vector.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.

     
    Output:
     - <list> : list  with the keydots points  
     
    """
    return 0


def Evaluate_S_Statistic(mediatrix_data, sigma_out=True,sigma=1,sigma_pre=0.5,Area_out=True,image_Path=''):

    """
    Function to calculate the S estatistic measurements.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.

     
    Output:
     - <dic> : dictionary with the S statistic measurements. The keys are  
     
    """
    guess=mediatrix_data[0]['center']
    L=0
    theta=[]
    linear=[]
    for i in range(0,len(mediatrix_data)):
        theta.append(mediatrix_data[i]['theta'])
        linear.append(mediatrix_data[i]['linear_coefficient'])
        L=L+mediatrix_data[i]['modulus']
        
    MinM=fmin(M_function,guess,args=(theta,linear),maxiter=1000, disp=0)
    MinM_val=M_function(MinM,theta,linear)
    R=(guess[0]-MinM[0])**2 + (guess[0]-MinM[0])**2
    R=sqrt(R)
    S_output={'id': mediatrix_data[0]['id']}
    try:
        S_output['MinM_norm']=MinM_val/(L*L)
    except:
        print " Impossible to define a lenght for "+str(mediatrix_data[0]['id'])
        S_output['MinM_norm']=-1
    try:
        S_output['L/R']=L/R
    except:
        print "Impossible to define a radius and curvature for "+str(mediatrix_data[0]['id'])
        S_output['L/R']=-1
    if sigma_out==True:
        lim=round((L+2.)/2.,2)
	sigma_X=[]
        sigma_Y=[]
        X_search=arange(round(MinM[0],2)-lim, round(MinM[0],2)+lim,sigma_pre*1)
	Y_search=arange(round(MinM[1],2)-lim, round(MinM[1],2)+lim,sigma_pre*1)
	for i in range(0,len(X_search)):
		for j in range(0,len(Y_search)):
				p=[X_search[i],Y_search[j]]
				M=  M_function(p,theta,linear) - MinM_val
				if M<=sigma:
					X_aux=round(X_search[i],0)
					Y_aux=round(Y_search[j],0)
					p_aux=[X_aux,Y_aux]
					rep=0
					for k in range(0,len(sigma_X)):
						if sigma_X[k]==X_aux and sigma_Y[k]==Y_aux:
							rep+=1
							
					if rep==0:
						sigma_X.append(X_aux)
						sigma_Y.append(Y_aux)
	E1,E2=get_extrema_2loops(sigma_X, sigma_Y, 0 )
        Extreme_sigma1=[sigma_X[E1],sigma_Y[E1]]
        Extreme_sigma2=[sigma_X[E2],sigma_Y[E2]]
        Extreme=[Extreme_sigma1,Extreme_sigma2]
        sigma_lenght=lenght(Extreme)
        try:
            S_output['sigma_lenght']=(sigma_lenght)/L
        except:
            print "Impossible to define sigma_lenght for "+str(mediatrix_data[0]['id'])
        sigma_minor_axis=len(sigma_X)/(4*atan(1)*sigma_lenght)
        try:
            S_output['sigma_exc']=(sigma_lenght)/sigma_minor_axis
        except:
            print "Impossible to define sigma_exc for "+str(mediatrix_data[0]['id'])
        if Area_out==True:
            image_name=mediatrix_data[0]['id']
            image,hdr = pyfits.getdata(image_Path+image_name, header = True )
            pixels=numpy.where(image>0)
            sigma_points=[sigma_X,sigma_Y]
            intersection=Area_intersection(pixels,sigma_points)
            S_output['intersection']=intersection
                
    return S_output
    


def M_function(p,theta,b):
    M=0
    for i in range(0,len(theta)):
	aux=Distance_from_line_to_point(p,theta[i],b[i])
        M=M+(aux*aux)
    M=M/float(len(theta))
    return M

def Area_intersection(a_coords,b_coords):
	common=0
	for i in range(0,len(a_coords[0])):
		for k in range(0,len(b_coords[0])):
			if (a_coords[0][i]==b_coords[0][k]) and (a_coords[1][i]==b_coords[1][k]):
				 common+=1
	
	return common



