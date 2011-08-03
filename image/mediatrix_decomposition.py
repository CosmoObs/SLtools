"""
Set of functions used in Mediatrix Decomposition.

"""



##@package mediatrix_decomposition 
#
# This module includes all functions developed to perform mediatrix decomposition method.
#
#



from sltools.geometry.get_extrema_pts import get_extrema_2loops
from numpy import where
from math import sin, cos ,sqrt, fabs, atan, tan 
from pylab import subplot, Rectangle, Arrow, xlim, ylim, ylabel, xlabel, title, savefig, Circle
from pyfits import getdata
from sltools.geometry.elementary_geometry import define_perpendicular_bisector, get_distance_from_line_to_point, length_from_connected_dots, two_points_to_line, three_points_to_circle, width_ellipse

def find_keydots (p1,p2,image_pixels,image,keydots,area, method="medium",alpha=1,n=-1,near_distance=(sqrt(2)/2)):
    """
    Function to calculate the keydot points in Mediatrix Decomposition.

    This function performs the full Mediatrix Decomposition method and returns 
    a numpy array with the [x,y] coordinates of all the keydot points. The 
    input 'method' allows the user to choose how to define the mediatrix 
    point (medium between width extrema or the brightest point along the 
    perpendicular bisector). 

    To find all the keydots, a recursion is used. To choose between the two 
    possible criteria to stop the recursion, the user must either set the 
    input 'alpha' to zero and choose a value for the input 'n' or leave 'n' 
    at its default value -1 and choose alpha. The default of the function 
    is to work with 'alpha' = 1.
    
    Input:
     - p1            <ndarray> : coordinates (x,y) of the first extreme point.
     - p2            <ndarray> : coordinates (x,y) of the second extreme point.
     - image_pixels     <list> : list of points' coordinates fitting the object.
     - image         <ndarray> : the image matrix.
     - keydots       <ndarray> : array of coordinates [p_i_x,p_i_y] of the two 
                                 extreme points.
     - area          <ndarray> : the object area.
     - method         <string> : possible values are 'medium' or 'brightest'.
     - alpha           <float> : the factor alpha=l_i/w to stop the bisection.
     - n                 <int> : the number of level iterations.
     - near_distance   <float> : the distance (in pixels) to consider a point close 
                                 to the perpendicular bisector.
     
    Output:
     - keydots       <ndarray> : array with the coordinate pairs of all the keydots.  
     
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
    L=length_from_connected_dots(keydots)
    W=width_ellipse(L,area)
    coefficients=define_perpendicular_bisector(p1,p2)
    p3x,p3y,p3Flag=choose_near_point(coefficients[0],coefficients[1],image_pixels,image,method,near_distance)
    p3=[p3x,p3y]
    if n==-1:
        if dl>(alpha*W) and len(keydots)<100:
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
    else:
        if n>0:
            n=n-1
            if (p3Flag==0):		
                if (not(p3 in keydots)):
                    keydots.insert(indexNext,p3)
                    keydots=find_keydots(p1,p3,image_pixels,image,keydots,area, method,alpha=0,n,near_distance)
                    keydots=find_keydots(p3,p2,image_pixels,image,keydots,area, method,alpha=0,n,near_distance)
            else:
                xmed=float(x1+x2)/2.
                ymed=float(y1+y2)/2.
                pmed=[xmed,ymed]
                if p1 in keydots: 
                    keydots=find_keydots(p1,pmed,image_pixels,image,keydots,area, method,alpha=0,n,near_distance)
                if p2 in keydots:
                    keydots=find_keydots(pmed,p2,image_pixels,image,keydots,area, method,alpha=0,n,near_distance)



    return keydots

def run_mediatrix_decomposition(image_name,image_dir='', method="medium",alpha=1,n=-1,near_distance=(sqrt(2)/2)):
    """
    Function to perform the mediatrix decomposition method on a given object. 

    Input:
     - image_name        <str> : the image file name.
     - image_dir         <str> : the image directory. If it is on the same 
                                 directory, directory=''.
     - method         <string> : possible values are 'medium' or 'brightest'.
     - alpha           <float> : the factor alpha=l_i/w to stop the bisection.
     - n                 <int> : the number of level iterations.
     - near_distance   <float> : the distance (in pixels) to consider a point 
                                 close to the perpendicular bisector.
     
    Output:
     - <dic> : Dictionary structure. The two keys 'origin' and 'end' represents each mediatrix vector. The 'origin' key is a list of points (x,y) of the vector origin,  and 'end'  key is a list of points (x,y) of the vector end. The 'id' key contains the image_name and 'center' key the object center (x,y) defined by the first mediatrix point in the first mediatrix level. The 'circle_params' key is a list where '0' index is the center (x,y) of the circle defined by two extrema and object center, the '1' and '2' index are the two extrema (x,y)
         
    """
    image,hdr = getdata(image_dir+image_name, header = True )
    pixels=where(image>0)
    e_1,e_2=get_extrema_2loops(pixels[0], pixels[1], 0 )
    Area=len(pixels[1])
    p1=[pixels[0][e_1],pixels[1][e_1]] # the extreme points p_1 and p_3
    p3=[pixels[0][e_2],pixels[1][e_2]]
    keydots=[p1,p2]
    keydots=find_keydots(p1,p3,pixels,image,keydots,Area, method=method,alpha=alpha,n,near_distance=near_distance)
    mediatrix_vectors=find_mediatrix_vectors(keydots)
    mediatrix_vectors['id']=image_name
    medium=int(float(len(keydots))/2)
    mediatrix_vectors['center']=keydots[medium]
    L=length_from_connected_dots(keydots)
    W=len(pixels[0])/(atan(1)*L)
    mediatrix_vectors['L/W']=L/W
    mediatrix_vectors['L']=L
    p2=[mediatrix_vectors['center'][0],mediatrix_vectors['center'][1]]
    x_c,y_c,r=three_points_to_circle(p1,p2,p3)
    circle_center=[x_c,y_c]
    mediatrix_vectors['circle_params']=[circle_center,p1,p2]

    return mediatrix_vectors


def get_length_from_mediatrix(image_name,image_dir='', method="medium",alpha=1,near_distance=(sqrt(2)/2)):
    """
    Function to calculate the length of an object using the Mediatrix Decomposition.

    The length is defined as the sum of the distances between the keydots.

    Input:
     - image_name      <array> : the image file name.
     - method         <string> : possible values are 'medium' or 'brightest'.
     - alpha           <float> : the factor alpha=l_i/w.
     - near_distance   <float> : the distance to consider a point 
                                 close to the perpendicular bisector.
     
    Output:
     - L <float> : the object length. 
     
    """
    image,hdr = getdata(image_dir+image_name, header = True )
    pixels=where(image>0)
    E1,E2=get_extrema_2loops( pixels[0], pixels[1], 0 )
    Area=len(pixels[1])
    p1=[pixels[0][E1],pixels[1][E1]]
    p2=[pixels[0][E2],pixels[1][E2]]
    keypoints=[p1,p2]
    keypoints=find_keydots(p1,p2,pixels,image,keypoints,Area, method="brightest",alpha=alpha,near_distance=near_distance)
    L_f=length_from_connected_dots(keypoints)
    W=width_ellipse(L_f,Area)
    
    return L_f, W


def find_mediatrix_vectors(points): 
    """
    From a given set of points, this function returns the mediatrix decomposition 
    vector between those points.

    The output is a dictionary with keys 'origin' and 'end'. Each key corresponds
    to a list of tuples, each tuple containing the x and y coordinate-values of a 
    given point. Therefore, the 'origin' key gives all vector origin points and 
    respectively for 'end'.
  
    Input:
     - points <list> : list of tuples of float - p_i mediatrix points where 
                       p_i=(x_i,y_i).
     
    Output:
     - mediatrix_vectors <dic> : Dictionary with two <string> keys for lists
                                 of tuples of floats.
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

        mediatrix_vectors['origin'].append(origin) 
        mediatrix_vectors['end'].append(end)
        
                
        
    return mediatrix_vectors

def print_mediatrix_Object_graph(mediatrix_data,image_dir='', keydots=False, colors= {'object': "g", 'vector': "b", 'keydots': "k"},alpha=1, plot_title="Mediatrix Decomposition graph", save=True, save_dir=''):

    """
    Make a plot presenting the object, keydots and mediatrix vectors. 

    Input:
    - mediatrix_data     <list> : the output from mediatrix_decomposition.
    - image_dir           <str> : the image directory. If it is on the same directory, 
                                  directory=''.
    - keydots            <bool> : 'True' -> display keydots, 'False' -> do not display. 
    - colors              <dic> : set the plot colors. The possible keys are 'object', 
                                  'vector' and 'keydots'.  
     
    Output:
     - <bool> or <subplot>: If save=='True' it returns a <bool> that indicates if the graph was su  if the graph was saved or 'false' if it was not. If sava=='False' retunrs the subplot instance. 
         
    """
    
    image_name=mediatrix_data['id']
    image,hdr = getdata(image_dir+image_name, header = True )
    pixels=where(image>0)    
    A = subplot(111)
    for i in range (0,len(pixels[0])):
        xy=[pixels[1][i]-0.5,pixels[0][i]-0.5]
        rec=Rectangle(xy, 1, 1, ec=colors['object'], fc=colors['object'], zorder=100)
        A.add_patch(rec)
   
    length=0
    
      
    if keydots==True:
        e_1,e_2=get_extrema_2loops(pixels[0], pixels[1], 0 )
        area=len(pixels[1])
        p1=[pixels[0][e_1],pixels[1][e_1]] # the extreme points p_1 and p_2
        p2=[pixels[0][e_2],pixels[1][e_2]]
        keydots=[p1,p2]
        keydots=find_keydots(p1,p2,pixels,image,keydots,area, method="brightest",alpha)
        key_x=[]
        key_y=[]
        for j in range(0,len(keydots)):
            key_x.append(keydots[j][0])
            key_y.append(keydots[j][1])
        
        A.plot(key_y,key_x,colors['keydots']+'.',zorder=500)
        

    
    for i in range(0,len(mediatrix_data['origin'])):
        origin_x=mediatrix_data['origin'][i][0]
        origin_y=mediatrix_data['origin'][i][1]
        end_x=mediatrix_data['end'][i][0]
        end_y=mediatrix_data['end'][i][1]
        length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
        length=length+sqrt(length_aux)
        d_x= end_x - origin_x
        d_y= mediatrix_data['end'][i][1] - mediatrix_data['origin'][i][1]
        arr = Arrow(origin_y, origin_x, d_y, d_x, width=0.05*length, fc=colors['vector'], ec='none',zorder=1000)
        A.add_patch(arr)
   
    xmin, xmax = xlim()
    ymin, ymax = ylim()
    min_inc_axis=40
    A.axis("equal")
    A.set_xlim(xmin-1*length,xmax+1*length)
    A.set_ylim(ymin-1*length,ymax+1*length)    
    ylabel("Y")
    xlabel("X")
    title(plot_title) 
    
    if save==True:
        savefig(save_dir+image_name+"_mediatrixGraph.png")
        A.clear()
        return True
    else:
        return A



def print_mediatrix_Object_circle_graph(mediatrix_data,image_dir='', keydots=False, colors= {'object': "g", 'vector': "b", 'keydots': "k"}, mediatrix_vectors=False, alpha=1, plot_title="Mediatrix Decomposition graph", save=True, save_dir=''):
    """
    Make a plot presenting the object, keydots and mediatrix vectors. 

    Input:
    - mediatrix_data  <list> : the output from mediatrix_decomposition.
    - image_dir        <str> : the image directory. If it is on the same directory directory=''.
    - keydots         <bool> : 'True' if you want to display the keydots and 'False' if you do not. 
    - colors           <dic> : set the plot colors. The possible keys are 'object', 'vector' and 'keydots'.       
    Output:
    - <bool> or <subplot>: If save=='True' it returns a <bool> that indicates if the graph was su  if the graph was saved or 'false' if it was not. If sava=='False' retunrs the subplot instance.   
         
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
    length=0
    for i in range(0,len(mediatrix_data['origin'])):
        origin_x=mediatrix_data['origin'][i][0]
        origin_y=mediatrix_data['origin'][i][1]
        end_x=mediatrix_data['end'][i][0]
        end_y=mediatrix_data['end'][i][1]
        length_aux=(origin_x - end_x)**2 + (origin_y - end_y)**2
        length=length+ sqrt(length_aux)
        if mediatrix_vectors==True:
            d_x= end_x - origin_x
            d_y= mediatrix_data['end'][i][1] - mediatrix_data['origin'][i][1]
            arr = Arrow(origin_y, origin_x, d_y, d_x, width=0.05*Length, fc=colors['vector'], ec='none',zorder=1000)
            A.add_patch(arr)
      
    if keydots==True:
        e_1,e_2=get_extrema_2loops(pixels[0], pixels[1], 0 )
        area=len(pixels[1])
        p1=[pixels[0][e_1],pixels[1][e_1]] # the extreme points p_1 and p_2
        p2=[pixels[0][e_2],pixels[1][e_2]]
        keydots=[p1,p2]
        keydots=find_keydots(p1,p2,pixels,image,keydots,area, method="brightest",alpha=1)
        key_x=[]
        key_y=[]
        for j in range(0,len(keydots)):
            key_x.append(keydots[j][0])
            key_y.append(keydots[j][1])
        
        A.plot(key_y,key_x,colors['keydots']+'.',zorder=500)
        

    
    #last=len(mediatrix_data['origin'])-1
    p1=[pixels[0][e_1],pixels[1][e_1]] # the extreme points p_1 and p_3
    p2=[mediatrix_data['center'][0],mediatrix_data['center'][1]]
    p3=[pixels[0][e_2],pixels[1][e_2]]
    x_c,y_c,r=three_points_to_circle(p1,p2,p3)
    if r>0:
        xy=[y_c,x_c]
        cir=Circle(xy,r,fc='none',ec='m', zorder=501)
        A.add_patch(cir)
    else:
        print "impossible to define a circle in" +str(mediatrix_data['id'])
      

   
    xmin, xmax = xlim()
    ymin, ymax = ylim()
    min_inc_axis=40
    A.axis("equal")
    A.set_xlim(xmin-1*length,xmax+1*length)
    A.set_ylim(ymin-1*length,ymax+1*length)    
    ylabel("Y")
    xlabel("X")
    #A.axis("equal")
    title("Mediatrix Decomposition graph") 
    
    if save==True and r>0:
        savefig(save_dir+image_name+"_mediatrixGraph_Circle.png")
        A.clear()
        return True
    else:
        return A


def choose_near_point(nearsX,nearsY,image,method='brightest'):
    
    """
    Function to choose the mediatrix point from a sample of selected points near a perpendicular bisector.

    Input:
     - near_xs        <list> : straight-line angle with respect to the x-axis.
     - nears_ys       <list> : linear coefficient of the line.
     - image       <ndarray> : the image matrix.
     - method       <string> : possible methods are 'medium'  or 'brightest'.
          
    Output:
     - <list> : the chosen point x coordinate.  
     - <list> : the chosen point y coordinate.
     -  <int> : a Flag error. if Flag=1, there is not any close point.
    """

    if (len(near_xs)<1):
        flag_err=1
        chosen_x=0
        chosen_y=0
        return chosen_x, chosen_y, flag_err
        if method=='brightest':
            chosen_x=near_xs.pop()
            chosen_y=near_ys.pop()
            flag_err=0
            while len(near_xs)>=1:
                chosen_aux_x=near_xs.pop()
                chosen_aux_y=near_ys.pop()
                if (image[chosen_x][chosen_y])<(image[chosen_aux_x][chosen_aux_y]):
                    chosen_x=chosen_aux_x
                    chosen_y=chosen_aux_y
            
        elif method=='medium':
            i,j=get_extrema_2loops(near_xs, near_ys, 0 )
            chosenX=float(near_xs[i]+near_xs[j])/2.
            chosenY=float(near_ys[i]+near_ys[j])/2.
            flagErr=0
    else:
        flagErr=1
        chosenX=0
        chosenY=0
    return chosen_x,chosen_y,flag_err

def pick_near_points(theta,b,object_pixels,near_distance): 
    """
    From a sample of points this function returns the ones that are close to a perpendicular bisector.

    Input:
     - theta           <float> : straight-line angle with respect to the x axis.
     - b               <float> : linear coefficient of the line.
     - object_pixels    <list> : a list of lists of floats, e.g. [x[i],y[i]].
     - near_distance   <float> : maximum distance to the perpendicular bisector.
     
    Output:
     - <list> : x-coordinates of the selected points.  
     - <list> : y-coordinates of the selected points.
     -  <int> : a Flag error. If Flag=1, there is no close point.
    """

    near_xs= []
    near_ys= []

    for i in range(0,len(object_pixels[1])):
        pixel=[object_pixels[0][i],object_pixels[1][i]]
	d=get_distance_from_line_to_point(pixel,theta,b)
        if(near_distance>=d):
            nearXs.append(object_pixels[0][i])
            nearYs.append(object_pixels[1][i])
    if len(nearXs)==0:
        flag_err=1
    else:
        flag_err=0

    return nearsX,nearsY,flag_err







