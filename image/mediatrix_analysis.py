from sltools.geometry.elementary_geometry import get_distance_from_line_to_point
from sltools.geometry.get_extrema_pts import get_extrema_2loops
from time import time
import os
import numpy as np
from math import sin, cos ,sqrt, fabs, atan 
from pylab import *
import pyfits
from scipy.optimize import fmin
from mediatrix_decomposition import get_length, get_length_c,plot_mediatrix, plot_mediatrixapl
import imcp

def Sum_inner_product(mediatrix_data):

    """
    Function to calculate the sum of scalar products from mediatrix vectors to its 'mirrored' mediatrix vector.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.

     
    Output:
     - <float> : sum of inner product between a mediatrix vector and its neighbour over the number of vectors 
     
    """
    total_inner_prod=0
    for i in range(0,len(mediatrix_data['end'])-1):
        
        component_1=mediatrix_data['end'][i][0]-mediatrix_data['origin'][i][0]
        component_2=mediatrix_data['end'][i][1]-mediatrix_data['origin'][i][1]
        vector1=[component_1,component_2]
        #print vector1
        component_1=mediatrix_data['end'][i+1][0]-mediatrix_data['origin'][i+1][0]
        component_2=mediatrix_data['end'][i+1][1]-mediatrix_data['origin'][i+1][1]

        vector2=[component_1,component_2]
        
        inner_prod=np.inner(vector1,vector2)
        total_inner_prod=total_inner_prod+inner_prod
    total_inner_prod=total_inner_prod/len(mediatrix_data['end'])
    return total_inner_prod


def Sum_inner_product_c(mediatrix_data):

    """
    Function to calculate the sum of scalar products from mediatrix vectors to its 'mirrored' mediatrix vector.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.   
    
        
    Output:
     - <float> : sum of inner product between a mediatrix vector and its neighbour over the number of vectors 
     
    """
    total_inner_prod=0
    for i in range(0,len(mediatrix_data['end'])-1):
        
        
        component=mediatrix_data['end'][i]-mediatrix_data['origin'][i]
        vector1=[(component.real/abs(component)),(component.imag/abs(component))]
        #print "complex"
        #print vector1
        
        component=mediatrix_data['end'][i+1]-mediatrix_data['origin'][i+1]
        vector2=[(component.real/abs(component)),(component.imag/abs(component))]
        
        
        inner_prod=np.inner(vector1,vector2)
        total_inner_prod=total_inner_prod+inner_prod
    total_inner_prod=total_inner_prod/len(mediatrix_data['end'])
    return total_inner_prod


def Sum_inner_product_cNNormalized(mediatrix_data):

    """
    Function to calculate the sum of scalar products from mediatrix vectors to its 'mirrored' mediatrix vector.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information of corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.   
    
        
    Output:
     - <float> : sum of inner product between a mediatrix vector and its neighbour over the number of vectors 
     
    """
    total_inner_prod=0
    for i in range(0,len(mediatrix_data['end'])-1):
        
        
        component=mediatrix_data['end'][i]-mediatrix_data['origin'][i]
        vector1=[component.real,component.imag]
        #print "complex"
        #print vector1
        
        component=mediatrix_data['end'][i+1]-mediatrix_data['origin'][i+1]
        vector2=[component.real,component.imag]
        
        
        inner_prod=np.inner(vector1,vector2)
        total_inner_prod=total_inner_prod+inner_prod
    total_inner_prod=total_inner_prod/len(mediatrix_data['end'])
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


def Evaluate_S_Statistic(mediatrix_data, sigma_out=True,sigma=1,sigma_pre=0.5,Area_out=True,image_path='', PS=False):

    """
    Function to calculate the S estatistic measurements.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.

     
    Output:
     - <dic> : dictionary with the S statistic measurements. The keys are  
     
    """
    guess=mediatrix_data['center']
    L=0
    theta=[]
    linear=[]
    for i in range(0,len(mediatrix_data['origin'])):
       origin_x=mediatrix_data['origin'][i][0]
       origin_y=mediatrix_data['origin'][i][1]
       end_x=mediatrix_data['end'][i][0]
       end_y=mediatrix_data['end'][i][1]
       Length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
       L=L+ sqrt(Length_aux)
       delta_x=float((end_x-origin_x ))
       if delta_x!=0:
           a=float((end_y-origin_y ))/delta_x
           theta.append(atan(a))
           b=end_y-a*(end_x)
           linear.append(b)
       else:
           theta.append(2*atan(1))
           linear.append(end_y-origin_y)
    MinM=fmin(M_function,guess,args=(theta,linear),maxiter=1000, disp=0)
    MinM_val=M_function(MinM,theta,linear)
    R=(guess[0]-MinM[0])**2 + (guess[1]-MinM[1])**2
    R=sqrt(R)
    center_comparison=(mediatrix_data['circle_params'][0][0]-MinM[0])**2 + (mediatrix_data['circle_params'][0][1]-MinM[1])**2
    center_comparison=sqrt(center_comparison)
    alpha=Find_angle_from_circle_section(mediatrix_data['circle_params'][0],mediatrix_data['circle_params'][1],mediatrix_data['circle_params'][2])
    S_output={'id': mediatrix_data['id']}
    try:
        S_output['MinM_norm']=MinM_val/(L*L)
    except:
        print " Impossible to define a lenght for "+str(mediatrix_data['id'])
        S_output['MinM_norm']=-1
    try:
        S_output['L/R']=L/R
    except:
        print "Impossible to define a radius and curvature for "+str(mediatrix_data['id'])
        S_output['L/R']=-1
    try:
        S_output['curvature_comparison']=(L/R)/alpha
    except:
        print "Impossible to compare curvature from circle section and mediatrix S for "+str(mediatrix_data['id'])
        S_output['curvature_comparison']=-1
    try:
        S_output['center_comparison']=center_comparison/L
    except:
        print "Impossible to compare center from circle method and mediatrix S for "+str(mediatrix_data['id'])
        S_output['curvature_comparison']=-1
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
        sigma_lenght=get_length(Extreme)
        try:
            S_output['sigma_lenght']=(sigma_lenght)/L
        except:
            print "Impossible to define sigma_lenght for "+str(mediatrix_data['id'])
        sigma_minor_axis=len(sigma_X)/(4*atan(1)*sigma_lenght)
        try:
            S_output['sigma_exc']=(sigma_lenght)/sigma_minor_axis
        except:
            print "Impossible to define sigma_exc for "+str(mediatrix_data['id'])
        if Area_out==True:
            
            if PS==True:
                #this is gambiarra
                image_name_list= mediatrix_data['file_seg_name'].split("/")
                image_name_aux=image_name_list[len(image_name_list)-1]
                image_name_aux=image_name_aux.replace("_seg.fits", "")
                #end of gambiarra the next line is correct
                #image_name_aux=mediatrix_data['file_seg_name'].replace("_seg.fits", "")
                
                image_name=image_name_aux+"_stamp_"+str(mediatrix_data['id'])+".fits"
                image,hdr = pyfits.getdata(image_path+image_name, header = True )
                pixels=np.where(image>0)
            else:
                image_name=mediatrix_data['file_seg_name']
                image,hdr = pyfits.getdata(image_path+image_name, header = True )
                pixels=np.where(image==int(mediatrix_data['id']))
            
            sigma_points=[sigma_X,sigma_Y]
            intersection=Area_intersection(pixels,sigma_points)
            S_output['intersection']=intersection
                
    return S_output



def Evaluate_S_Statistic_on_matrix(mediatrix_data,obj_stamp, sigma_out=True,sigma=1,sigma_pre=0.5,Area_out=True):

    """
    Function to calculate the S estatistic measurements.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.

     
    Output:
     - <dic> : dictionary with the S statistic measurements. The keys are  
     
    """
    guess=mediatrix_data['center']
    L=0
    theta=[]
    linear=[]
    for i in range(0,len(mediatrix_data['origin'])):
       origin_x=mediatrix_data['origin'][i][0]
       origin_y=mediatrix_data['origin'][i][1]
       end_x=mediatrix_data['end'][i][0]
       end_y=mediatrix_data['end'][i][1]
       Length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
       L=L+ sqrt(Length_aux)
       delta_x=float((end_x-origin_x ))
       if delta_x!=0:
           a=float((end_y-origin_y ))/delta_x
           theta.append(atan(a))
           b=end_y-a*(end_x)
           linear.append(b)
       else:
           theta.append(2*atan(1))
           linear.append(end_y-origin_y)
    MinM=fmin(M_function,guess,args=(theta,linear),maxiter=1000, disp=0)
    MinM_val=M_function(MinM,theta,linear)
    R=(guess[0]-MinM[0])**2 + (guess[1]-MinM[1])**2
    R=sqrt(R)
    #print "R real"
    #print R
    center_comparison=(mediatrix_data['circle_params'][0][0]-MinM[0])**2 + (mediatrix_data['circle_params'][0][1]-MinM[1])**2
    center_comparison=sqrt(center_comparison)
    alpha=Find_angle_from_circle_section(mediatrix_data['circle_params'][0],mediatrix_data['circle_params'][1],mediatrix_data['circle_params'][2])
    S_output={'id': mediatrix_data['id']}
    try:
        S_output['MinM_norm']=MinM_val/(L*L)
    except:
        print " Impossible to define a lenght for "+str(mediatrix_data['id'])
        S_output['MinM_norm']=-1
    try:
        S_output['L/R']=L/R
    except:
        print "Impossible to define a radius and curvature for "+str(mediatrix_data['id'])
        S_output['L/R']=-1
    try:
        S_output['curvature_comparison']=(L/R)/alpha
    except:
        print "Impossible to compare curvature from circle section and mediatrix S for "+str(mediatrix_data['id'])
        S_output['curvature_comparison']=-1
    try:
        S_output['center_comparison']=center_comparison/L
    except:
        print "Impossible to compare center from circle method and mediatrix S for "+str(mediatrix_data['id'])
        S_output['curvature_comparison']=-1
    if sigma_out==True:
        lim=round((L+2.)/2.,2)
	sigma_X=[]
        sigma_Y=[]
        X_search=arange(round(MinM[0],2)-lim, round(MinM[0],2)+lim,sigma_pre*1)
	Y_search=arange(round(MinM[1],2)-lim, round(MinM[1],2)+lim,sigma_pre*1)
        #print "Real lenght for sigma search x and y"
        #print len(X_search)
        #print len(Y_search)
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
        sigma_lenght=get_length(Extreme)
        try:
            S_output['sigma_lenght']=(sigma_lenght)/L
        except:
            print "Impossible to define sigma_lenght for "+str(mediatrix_data['id'])
        sigma_minor_axis=len(sigma_X)/(4*atan(1)*sigma_lenght)
        #print "Area de sigma, real"
        #for i in range(0,len(sigma_X)):
        #    for j in range(0,len(sigma_Y)):
        #        print sigma_X[i],sigma_Y[j]
        try:
            S_output['sigma_exc']=(sigma_lenght)/sigma_minor_axis
        except:
            print "Impossible to define sigma_exc for "+str(mediatrix_data['id'])
        if Area_out==True:
            
            pixels=np.where(obj_stamp>0)
            sigma_points=[sigma_X,sigma_Y]
            intersection=Area_intersection(pixels,sigma_points)
            S_output['intersection']=intersection
                
    return S_output




def Evaluate_S_Statistic_on_matrix_c(mediatrix_data,obj_stamp, sigma_out=True,sigma=1,sigma_pre=0.5,Area_out=True):

    """
    Function to calculate the S estatistic measurements.
    
    Input:
     - mediatrix_data <list> : a list of dictionary structure. Each list item is a dictionary with information corresponding to a mediatrix vector. The keys are 'theta' for the angle with x axis, 'linear_coefficient' for the linear coefficient from the line in the vector direction, 'origin' the point (x,y) of the vector origin, 'end' the point (x,y) of the vector, 'modulus' for vector modulus. The first item from the list has an extra key 'id' wich contains the image file name. It is the output from Mediatrix_Decomposition.

     
    Output:
     - <dic> : dictionary with the S statistic measurements. The keys are  
     
    """
    guess=[mediatrix_data['center'].real,mediatrix_data['center'].imag]
    L=0
    theta=[]
    linear=[]
    L=mediatrix_data['L']
    #print " E o L eh"
    #print L
    for i in range(0,len(mediatrix_data['origin'])):
       origin=mediatrix_data['origin'][i]
       end=mediatrix_data['end'][i]
       delta=end-origin
       if delta.real!=0:
           a=float((delta.imag ))/delta.real
           theta.append(atan(a))
           b=end.imag-a*(end.real)
           linear.append(b)
       else:
           theta.append(2*atan(1))
           linear.append(end.imag-origin.imag)
    minM_r=fmin(M_function,guess,args=(theta,linear),maxiter=1000, disp=0)
    minM_val=M_function(minM_r,theta,linear)
    minM=minM_r[0]+minM_r[1]*1j
    R=abs(mediatrix_data['center']-minM)
    #print "esse e o ponto em que m eh minima"
    #print minM
    circle_center=mediatrix_data['circle_params'][0]
    #print "Este e o centro do circulo"
    #print circle_center
    center_comparison=abs(circle_center-minM)
    alpha=Find_angle_from_circle_section_c(mediatrix_data['circle_params'][0],mediatrix_data['circle_params'][1],mediatrix_data['circle_params'][2])
    S_output={'init': 0}
    try: 
        S_output['MinM_norm']=minM_val/(L*L)
    except:
        print " Impossible to define a lenght for "+str(mediatrix_data['id'])
        S_output['MinM_norm']=-1
    if R!=0:
        S_output['L/R']=L/R
    else:
        print "Impossible to define a radius and curvature for "+str(mediatrix_data['id'])
        S_output['L/R']=-1
    if R!=0 and alpha!=0:
        S_output['curvature_comparison']=(L/R)/alpha
    else:
        print "Impossible to compare curvature from circle section and mediatrix S for "+str(mediatrix_data['id'])
        S_output['curvature_comparison']=-1
    try:
        S_output['center_comparison']=center_comparison/L
    except:
        print "Impossible to compare center from circle method and mediatrix S for "+str(mediatrix_data['id'])
        S_output['curvature_comparison']=-1
    if sigma_out==True:
        
	search=np.array([])
        
        rc_center=minM
 
        
        sigma_points=[]
        step=sigma_pre
        
        lim=step
        sigma_points=find_rc_region_c(sigma_points,rc_center,step,lim=lim,theta=theta,linear=linear,minM=minM_val,sigma=sigma) 
        
        
        sigma_points=np.array(sigma_points)
        sigma_points=np.unique(np.round(sigma_points.real,0)+np.round(sigma_points.imag,0)*1j)


        E1,E2=get_extrema_2loops(sigma_points.real, sigma_points.imag, 0 )
        Extreme=[sigma_points[E1],sigma_points[E2]]
        sigma_lenght=get_length_c(Extreme)
        try:
            S_output['sigma_lenght']=(sigma_lenght)/L
        except:
            print "Impossible to define sigma_lenght for "+str(mediatrix_data['id'])
        sigma_minor_axis=len(sigma_points)/(4*atan(1)*sigma_lenght)
        #verificar essa eq.
        #print "Area de sigma, complex"
        #print sigma_points
        try:
            S_output['sigma_exc']=(sigma_lenght)/sigma_minor_axis
        except:
            print "Impossible to define sigma_exc for "+str(mediatrix_data['id'])
        if Area_out==True:
            
            pixels=np.where(obj_stamp>0)
            pixels_c=np.array(pixels[0]+pixels[1]*1j)
            intersection_points=np.intersect1d(pixels_c,sigma_points)
            S_output['intersection']=len(intersection_points)
                
    return S_output


def find_rc_region_c(region,center,step,lim,theta,linear,minM,sigma):
    
    nex_level=0
    ini_region=len(region)
    X_search=np.arange(round(center.real,2)-lim, round(center.real,2)+lim,step)
    Y_search=np.arange(round(center.imag,2)-lim, round(center.imag,2)+lim,step)
    
    search=np.array([])
    #search_len=len(X_search)*len(Y_search)
   
    #for i in range(0,len(X_search)):
    #    piece=np.repeat(X_search[i],len(Y_search))+Y_search*1j
    #    search=np.concatenate((piece,search),axis=1)
    len_search=0
    piece=np.repeat(X_search[0],len(Y_search))+Y_search*1j
    search=np.concatenate((piece,search),axis=1)
    len_search+=len(piece)

    piece=np.repeat(X_search[-1],len(Y_search))+Y_search*1j
    search=np.concatenate((piece,search),axis=1)
    len_search+=len(piece)
    
    piece=X_search+np.repeat(Y_search[0],len(X_search))*1j
    search=np.concatenate((piece,search),axis=1)
    len_search+=len(piece)    

    piece=X_search+np.repeat(Y_search[-1],len(X_search))*1j
    search=np.concatenate((piece,search),axis=1)
    len_search+=len(piece)

    #search_len=len(X_search)*4
    search=np.reshape(search,len_search)
    #search=np.unique(search)     
    m_values_c=np.apply_along_axis(M_function_c,0,search,theta,linear)
    m_values_c=m_values_c- np.repeat(minM,len(search))
    for i in range(0,len(search)):
            if m_values_c[i]<=sigma:
               region.append(search[i])
    
    if len(region)>ini_region and len(region) < 2000:
        region=find_rc_region_c(region,center,step,lim+step,theta,linear,minM,sigma)
    return region
    



    
def Find_angle_from_circle_section(circle_center,circle_point1,circle_point2):

    x_aux=float(circle_point1[0]+circle_point2[0])/2
    y_aux=float(circle_point1[1]+circle_point2[1])/2
    m_point=[x_aux,y_aux]
    dx=m_point[0]-circle_center[0]
    dy=m_point[1]-circle_center[1]
    adjacent=sqrt((dx*dx)+(dy*dy))
    dx=circle_point1[0]-circle_center[0]
    dy=circle_point1[1]-circle_center[1]
    opposite=sqrt((dx*dx)+(dy*dy))
    alpha=atan(opposite/adjacent)
    alpha=2*alpha 
    return alpha


def Find_angle_from_circle_section_c(circle_center,circle_point1,circle_point2):

    m_point=(circle_point1+circle_point2)/2
    delta_adjacent=abs(m_point-circle_center)
    delta_opposite=abs(circle_point1-circle_center)
    alpha=atan(delta_opposite/delta_adjacent)
    alpha=2*alpha 
    return alpha




def plot_S_Statistic(mediatrix_data, ps_name, sigma_out=True,sigma=1,sigma_pre=0.5,Area_out=True, save=True,out_image='',title=''):

    """
    Function to make plot of S estatistic interesting measurements.
    
    Input:
     - mediatrix_data <dic> : the output from mediatrix decomposition.

     
    Output:
     - <bool> :   
     
    """
 
    if out_image=='':
        out_image=ps_name.replace(".fits","")+"_mediatrixS_plot.png"
    if title=='':
        title="Mediatrix S"
    mediatrix_plot=plot_mediatrix(mediatrix_data,ps_name, keydots=True, colors= {'object': "g", 'vector': "b", 'keydots': "k"}, save=False)
    guess=[mediatrix_data['center'].real,mediatrix_data['center'].imag]
    Length=0
    theta=[]
    linear=[]
    for i in range(0,len(mediatrix_data['origin'])):
       origin_x=mediatrix_data['origin'][i].real
       origin_y=mediatrix_data['origin'][i].imag
       end_x=mediatrix_data['end'][i].real
       end_y=mediatrix_data['end'][i].imag
       Length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
       Length=Length+ sqrt(Length_aux)
       delta_x=float((end_x-origin_x ))
       if delta_x!=0:
           a=float((end_y-origin_y ))/delta_x
           theta.append(atan(a))
           b=end_y-a*(end_x)
           linear.append(b)
       else:
           theta.append(2*atan(1))
           linear.append(end_y-origin_y)
        
    MinM=fmin(M_function,guess,args=(theta,linear),maxiter=1000, disp=0)
    MinM_val=M_function(MinM,theta,linear)
    MinM_plot_X=[MinM[0],MinM[0]]
    MinM_plot_Y=[MinM[1],MinM[1]]
    mediatrix_plot.plot(MinM_plot_Y,MinM_plot_X,"b.", zorder=2000)
    lim=round((Length/8.)+5,2)
    sigma_X=[]
    sigma_Y=[]
    X_search=arange(round(MinM[0],2)-lim, round(MinM[0],2)+lim,sigma_pre*1)
    Y_search=arange(round(MinM[1],2)-lim, round(MinM[1],2)+lim,sigma_pre*1)
    for j in range(0,len(Y_search)):
        for i in range(0,len(X_search)):
            p=[X_search[i],Y_search[j]]
            M=  M_function(p,theta,linear) - MinM_val
            
            if M<=sigma:
                X_aux=round(X_search[i],2)
                Y_aux=round(Y_search[j],2)
                p_aux=[X_aux,Y_aux]
                rep=0
                for k in range(0,len(sigma_X)):
                    if sigma_X[k]==X_aux and sigma_Y[k]==Y_aux:
                        rep+=1
                if rep==0:
                    sigma_X.append(X_aux)
                    sigma_Y.append(Y_aux)
    duplicates=[]
    x_up=[]
    x_down=[]
    for i in range(0,len(sigma_Y)):
        if len(duplicates)!=sigma_Y.count(sigma_Y[i]):
            duplicates.append(sigma_X[i])
        if len(duplicates)==sigma_Y.count(sigma_Y[i]):
            duplicates=array(duplicates)
            x_up.append(duplicates.max())
            x_down.append(duplicates.min())
            duplicates=[]
    sigma_Y.sort()
    last = sigma_Y[-1]
    for i in range(len(sigma_Y)-2, -1, -1):
        if last == sigma_Y[i]:
            del sigma_Y[i]
        else:
            last = sigma_Y[i]


   # mediatrix_plot.plot(sigma_Y,x_up,"m.")
  #  mediatrix_plot.plot(sigma_Y,x_down,"m.")
    mediatrix_plot.fill_between(sigma_Y, x_up, x_down, where=None, color="m", zorder=999)
    
           
            #mediatrix_plot.plot(sigma_Y,sigma_X,"m.")

    if save==True:
        xmin, xmax = xlim()
        ymin, ymax = ylim()
        mediatrix_plot.set_xlim(xmin-1*Length,xmax+1*Length)
        mediatrix_plot.set_ylim(ymin-1*Length,ymax+1*Length)
        savefig(out_image)
        mediatrix_plot.clear()
        return True
    else:
        return mediatrix_plot




def plot_S_Statistic_apl(image_name,_id='',keydots=False,circle=False,rc=True, save=True,out_image='', args={}):

    """
    Function to make plot of S estatistic interesting measurements.
    
    Input:
     - mediatrix_data <dic> : the output from mediatrix decomposition.

     
    Output:
     - <bool> :   
     
    """
    opt={'increase': 2, 'relative_increase': True,'connected': False,'object_centered':True, 'type':'cutout', 'pmin':0.25 , 'pmax':99.75 , 'invert':True ,'out_title': 'Mediatrix Method', 'keys_color': "r" ,'alpha': 1 ,'max_level': 1000, 'near_distance': sqrt(2)/2, 'max_level': 1000, 'method':"brightest",'sigma':1,'sigma_pre':0.5, 'rc_color': 'm'}
    opt.update(args)
     
    if out_image=='':
        out_image=image_name.replace(".fits","")+"_mediatrixS_plot.png"
 

    

    mediatrix_plot,mediatrix_data,image_ps=plot_mediatrixapl(image_name,_id=_id, keydots=keydots,circle=circle, save=False, args=opt)  
    
    

    guess=[mediatrix_data['center'].real,mediatrix_data['center'].imag]
    Length=0
    theta=[]
    linear=[]
    for i in range(0,len(mediatrix_data['origin'])):
       origin_x=mediatrix_data['origin'][i].real
       origin_y=mediatrix_data['origin'][i].imag
       end_x=mediatrix_data['end'][i].real
       end_y=mediatrix_data['end'][i].imag
       Length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
       Length=Length+ sqrt(Length_aux)
       delta_x=float((end_x-origin_x ))
       if delta_x!=0:
           a=float((end_y-origin_y ))/delta_x
           theta.append(atan(a))
           b=end_y-a*(end_x)
           linear.append(b)
       else:
           theta.append(2*atan(1))
           linear.append(end_y-origin_y)
        
    MinM=fmin(M_function,guess,args=(theta,linear),maxiter=1000, disp=0)
    MinM_val=M_function(MinM,theta,linear)
    MinM_plot_X=[MinM[0],MinM[0]]
    MinM_plot_Y=[MinM[1],MinM[1]]
    
    mediatrix_plot.show_markers(MinM_plot_Y,MinM_plot_X,c='red',marker='.',zorder=2000)
    print 'min m'
    print MinM_plot_X
    print MinM_plot_Y
    lim=round((Length/8.)+5,2)
    sigma_X=[]
    sigma_Y=[]
    ps=[]
    Ms=[]
    X_search=arange(round(MinM[0],2)-lim, round(MinM[0],2)+lim,opt['sigma_pre']*1)
    Y_search=arange(round(MinM[1],2)-lim, round(MinM[1],2)+lim,opt['sigma_pre']*1)
    for j in range(0,len(Y_search)):
        for i in range(0,len(X_search)):
            p=[X_search[i],Y_search[j]]
            M=  M_function(p,theta,linear) - MinM_val
            ps.append(p)
            Ms.append(M)
            if M<=opt['sigma']:
                X_aux=round(X_search[i],2)
                Y_aux=round(Y_search[j],2)
                p_aux=[X_aux,Y_aux]
                rep=0
                for k in range(0,len(sigma_X)):
                    if sigma_X[k]==X_aux and sigma_Y[k]==Y_aux:
                        rep+=1
                if rep==0:
                    sigma_X.append(X_aux)
                    sigma_Y.append(Y_aux)
    duplicates=[]
    x_up=[]
    x_down=[]
    for i in range(0,len(sigma_Y)):
        if len(duplicates)!=sigma_Y.count(sigma_Y[i]):
            duplicates.append(sigma_X[i])
        if len(duplicates)==sigma_Y.count(sigma_Y[i]):
            duplicates=array(duplicates)
            x_up.append(duplicates.max())
            x_down.append(duplicates.min())
            duplicates=[]
    sigma_Y.sort()
    last = sigma_Y[-1]
    for i in range(len(sigma_Y)-2, -1, -1):
        if last == sigma_Y[i]:
            del sigma_Y[i]
        else:
            last = sigma_Y[i]
    
    
    if rc==True:
        #image_segname=image_name.replace(".fits","")+"_seg.fits"
        #image_objname=image_name.replace(".fits","")+"_obj.fits"

        #image_seg,hdr = pyfits.getdata(image_segname, header = True )
        #image_obj,hdr = pyfits.getdata(image_objname, header = True )
        #image_ps,hdr=imcp.segstamp(segimg=image_seg, objID=_id, objimg=image_obj, hdr=hdr, increase=increase, relative_increase=relative_increase, connected=False)


        pixels=where(image_ps>0)
        image_ps[pixels]=0
        pixels_up=[x_up,sigma_Y]
        pixels_down=[x_down,sigma_Y]
        lines=[]
        for j in range(0,len(sigma_Y)):
            lines.append([sigma_Y[j],x_up[j]]) 
        #for j in range(0,len(ps)):
            #image_ps[ps[j]]=Ms[j]   
        #image_ps[pixels_up]=1
        #image_ps[pixels_down]=1
        #mediatrix_plot.show_contour(image_ps,colors='red',levels=[1.0,1.5], smooth=1)#,filled=False, smooth=1)
    #A = subplot(111)
        mediatrix_plot.show_markers(sigma_Y,x_up,c=opt['rc_color'],marker=".", linewidths='0',zorder=999)
        mediatrix_plot.show_markers(sigma_Y,x_down,c=opt['rc_color'],marker=".",linewidths='0', zorder=999)
        for j in range(0,len(sigma_Y)):
           fill_lim = round(abs(x_up[j]-x_down[j]),1)
           for k in range(0,int(fill_lim)):
               mediatrix_plot.show_markers(sigma_Y[j],x_down[j]+k,c=opt['rc_color'],marker=".",linewidths='0', zorder=999)
               mediatrix_plot.show_markers(sigma_Y[j],x_down[j]+k-0.25,c=opt['rc_color'],marker=".",linewidths='0', zorder=999)
               mediatrix_plot.show_markers(sigma_Y[j],x_down[j]+k+0.25,c=opt['rc_color'],marker=".",linewidths='0', zorder=999) 
    
    #A.fill_between(sigma_Y, x_up, x_down, where=None, color="m", zorder=999)
    
           
    #mediatrix_plot.plot(sigma_Y,sigma_X,"m.")

    if save==True:
        #xmin, xmax = xlim()
        #ymin, ymax = ylim()
        #mediatrix_plot.set_xlim(xmin-1*Length,xmax+1*Length)
        #mediatrix_plot.set_ylim(ymin-1*Length,ymax+1*Length)
        savefig(out_image)
        #mediatrix_plot.clear()
        return True
    else:
        return mediatrix_plot










def M_function(p,theta,b):
    M=0
    for i in range(0,len(theta)):
        #print "este e o theta real"+str(theta[i])
        aux=get_distance_from_line_to_point(p,theta[i],b[i])
        M=M+(aux*aux)
        #print aux
    M=M/float(len(theta))
    return M


def M_function_c(p,theta,b):
    M=0
       
    for i in range(0,len(theta)):
        D=0
        #print "este e o theta imag"+str(theta[i])
        if(theta[i]!=2*atan(1)) and (theta[i]!=0.):
            
            a=tan(theta[i])
            m=-1./a
            c=p.imag-m*p.real
            x=(b[i]-c)/(m-a)
            y=a*x+b[i]
            p_line=x+y*1j
            D=abs(p-p_line)
        elif (theta[i]==2*atan(1)): # case of the perpendicular bisector is a vertical line a-->inf.
            D=abs(p.real-b[i])
        elif theta[i]==0.: # case of the perpendicular bisector is a horizontal line a-->0
           D=abs(p.imag-b[i])
        #print D
        M+=(D*D)
    M=M/float(len(theta))
    return M



 

def Area_intersection(a_coords,b_coords):
	common=0
	for i in range(0,len(a_coords[0])):
		for k in range(0,len(b_coords[0])):
			if (a_coords[0][i]==b_coords[0][k]) and (a_coords[1][i]==b_coords[1][k]):
				 common+=1
	
	return common



