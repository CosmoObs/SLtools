#!/usr/bin/env python

#
# Example boxplot code
#

from matplotlib import pyplot as plt;
import numpy as np;
import math;

nyticks = 4;

# ---
def scattergrid(points,sizes=[],colors=[]):
    """
    Plots a scatter graphic in a X,Y axes with sized symbols.

    Input:
     - points  <[(,)]> : a list of tuples of X,Y datapoints
     - sizes   <[]>    : a list of numbers for symbols sizing
     - colors  <[]>    : a list of numbers for symbols colors
    
    Output:
     - fig  <pyplot.figure()>  : A matplotlib figure instance
        Save (fig.savefig()) or show (fig.show()) the object.

    ---
    """

    if len(points) != len(sizes)  and  len(sizes) != 0:
        return False;

    RA,DEC = zip(*points);
    RA_min = min(RA);
    RA_max = max(RA);
    DEC_min = min(DEC);
    DEC_max = max(DEC);

    ra_ = abs(RA_max-RA_min)/20
    dec_ = abs(DEC_max-DEC_min)/20

    sizes = np.asarray(sizes)
    sizes = 100 * (sizes-sizes.min()) / (sizes.max()-sizes.min())
    sizes = sizes.astype('uint32').tolist();

    x,y = zip(*points);
    
    # Sort the data just to get it(the plot) right done when 
    # 'scatter' run through datasets.
    tosort_ = zip(sizes,colors,x,y);
    tosort_.sort();
    tosort_.reverse();
    sizes,colors,x,y = zip(*tosort_);
    #
    
    plt.set_cmap('hot');
    fig = plt.figure();
    ax = fig.add_subplot(111);
    ax.patch.set_facecolor('0.9');
    ax.grid(True);
    ax.set_xlabel('RA');
    ax.set_ylabel('DEC');
#    ax.set_title('Massive objects scatter (position) plot')

    ax.scatter(x,y,s=sizes,c=colors,marker='o',edgecolor='black',alpha=0.75)


    s2 = [sizes[0],sizes[-1]]
    x2 = [1,1]
    y2 = [1,2]
#    fig.show();

    return fig;


# ---
def histbox(data,nbins=None,outliers=True,legend=None,xlabel="X",ylabel="Normed counts",title="Distribution plot"):
    """
    Plots data distribution in two different ways: a boxplot and (over) a histogram.
    
    
    Input:
     - data       < [] > : A list of data values or (for comparative plots) a list of datasets(vectors)
     - nbins     < int > : Number of bins for histogram. If None (default) nbins will be in somewhere between 10 and 100
     - outliers < bool > : If False, outliers (points) will not be shown
     - legend     < [] > : Data label for each dataset
     - xlabel    < str > : Label for X axis
     - ylabel    < str > : Label for Y axis
     - title     < str > : Figure title

    Output:
     - fig < plt.figure() > : Matplotlib figure instance

    ---
    """

#    if type(data) is not np.ndarray:
#        data = np.array(data);
    try:
        ND = data.ndim;
        if ND == 2:
            ND = data.shape[1];
    except:
        ND = len(data);
    
    colors=['b', 'r', 'y', 'k', 'g', 'm', 'c'];
    
    # Bins width:
    #
    def binswidth(dado):
        width = 3.49 * np.std(dado) * math.pow(dado.size,-1/3.);
        return width;
    ##
    if not nbins:
        _dat = np.asarray(_dat);
        _w = binswidth(_dat.ravel());
        _min = _dat.ravel().min();
        _max = _dat.ravel().max();
        nbins = (_max-_min)/_w;
        nbins = min(100,max(10,nbins));
        del _dat,_w,_min,_max;
        
    # Initialize the plot environment:
    fig,(ax1,ax2) = plt.subplots(2,1,sharex=True);

    # BOXPLOT
    bpD = ax1.boxplot(data,vert=0,sym='*',patch_artist=True,notch=1)#,positions=bp_positions[:-1],bootstrap=None);
    _box = bpD['boxes'];
    _med = bpD['medians'];
    _cap = bpD['caps'];
    _wsk = bpD['whiskers'];
    _flr = bpD['fliers'];
    
    # Fix the colors for each boxplot plotted..
    xmin = [];
    xmax = [];
    for i in range(ND):
        _box[i].set(color=colors[i],linewidth=2,alpha=0.5);
        _wsk[2*i].set(color=colors[i])#,linestyle=':');
        _wsk[2*i+1].set(color=colors[i])#,linestyle=':');
        #_med[i].set(color=colors[i]);
        #_cap[2*i].set(color=colors[i],linewidth=2);
        #_cap[2*i+1].set(color=colors[i],linewidth=2);
        xmin.extend(_wsk[i].get_xdata());
        xmax.extend(_wsk[i+1].get_xdata());
    
    # HISTOGRAM
    num,bin,p = ax2.hist(data,bins=nbins,normed=True,rwidth=0.8,label=legend,color=colors[:ND],align='mid',histtype='bar',alpha=0.5,linewidth=1)#,range=xrange);
    
    # Adjust the plots
    fig.subplots_adjust(hspace=0);
    ax1.yaxis.set_visible(False);
    #ax1.xaxis.set_visible(False);
    ax1.xaxis.set_ticks_position('top');
    ax2.xaxis.set_ticks_position('bottom');
    ax1.xaxis.grid(True)
    ax2.yaxis.grid(True)
    ax1.patch.set_facecolor('0.9');
    ax2.patch.set_facecolor('0.9');

    if not outliers:
        ax2.set_xlim((min(xmin),max(xmax)));
    
    fig.suptitle(title)
    ax2.set_ylabel(ylabel)
    ax2.set_xlabel(xlabel)
    ax2.legend(loc='upper right')
    
    return fig;
    
