from geometry import *
from sltools.geometry.get_extrema_pts import *
import os
import numpy
from math import sin, cos ,sqrt, fabs, atan 
from pylab import *
import pyfits
from scipy.optimize import fmin
from mediatrix_decomposition import get_length, print_mediatrix_Object_graph

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

        component_1=mediatrix_data['end'][i+1][0]-mediatrix_data['origin'][i+1][0]
        component_2=mediatrix_data['end'][i+1][1]-mediatrix_data['origin'][i+1][1]

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


def Evaluate_S_Statistic(mediatrix_data, sigma_out=True,sigma=1,sigma_pre=0.5,Area_out=True,image_path=''):

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
    center_comparisson=(mediatrix_data['circle_params'][0][0]-MinM[0])**2 + (mediatrix_data['circle_params'][0][1]-MinM[0])**2
    center_comparisson=sqrt(center_comparisson)
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
        S_output['curvature_comparisson']=(L/R)/alpha
    except:
        print "Impossible to compare curvature from circle section and mediatrix S for "+str(mediatrix_data['id'])
        S_output['curvature_comparisson']=-1
    try:
        S_output['center_comparisson']=center_comparisson/L
    except:
        print "Impossible to compare center from circle method and mediatrix S for "+str(mediatrix_data['id'])
        S_output['curvature_comparisson']=-1
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
            image_name=mediatrix_data['id']
            image,hdr = pyfits.getdata(image_path+image_name, header = True )
            pixels=numpy.where(image>0)
            sigma_points=[sigma_X,sigma_Y]
            intersection=Area_intersection(pixels,sigma_points)
            S_output['intersection']=intersection
                
    return S_output
    
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
def Plot_S_Statistic(mediatrix_data, sigma_out=True,sigma=1,sigma_pre=0.5,Area_out=True,image_Path='', save=True,save_dir=''):

    """
    Function to make plot of S estatistic interesting measurements.
    
    Input:
     - mediatrix_data <dic> : the output from mediatrix decomposition.

     
    Output:
     - <bool> :   
     
    """


    mediatrix_plot=print_mediatrix_Object_graph(mediatrix_data,image_dir=image_Path, keydots=True, colors= {'object': "g", 'vector': "b", 'keydots': "k"}, save=False)
    guess=mediatrix_data['center']
    Length=0
    theta=[]
    linear=[]
    for i in range(0,len(mediatrix_data['origin'])):
       origin_x=mediatrix_data['origin'][i][0]
       origin_y=mediatrix_data['origin'][i][1]
       end_x=mediatrix_data['end'][i][0]
       end_y=mediatrix_data['end'][i][1]
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
        savefig(save_dir+mediatrix_data['id']+"_mediatrix_S_Estatistic_Graph.png")
        mediatrix_plot.clear()
        return True
    else:
        return mediatrix_plot




def M_function(p,theta,b):
    M=0
    for i in range(0,len(theta)):
	aux=get_distance_from_line_to_point(p,theta[i],b[i])
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



