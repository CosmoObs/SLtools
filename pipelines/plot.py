#!/usr/bin/env python

#
# Example boxplot code
#

from matplotlib import pyplot as plt;
import numpy as np;

nbins = 10;
nyticks = 4;

def histbox_plot(data=[],data_label=["data"],xlabel="X",ylabel="#",title="Title",xrange=None,yrange=None):
    """
    Input:
     - data : list of vectors
    """

    if len(data)==0:
        return False;
        
    ND = len(data);
    
    colors=['b', 'r', 'k', 'm', 'g', 'c', 'y'];
    
    if xrange == None:
        maxs = [ np.max(data[i]) for i in range(ND) ];
        maximo = max(maxs);
        mins = [ np.min(data[i]) for i in range(ND) ];
        minimo = min(mins);
        xrange=(int(minimo),int(maximo+1));
    if yrange == None:
        pass;
    
    f,ax = plt.subplots(1,1);

    # Hold on figure to plot Histogram and Boxplot together..
    #
    f.hold(True);

    # Histogram
    num,bin,p = ax.hist(data,bins=nbins,normed=True,rwidth=0.8,range=xrange,align='mid',color=colors[:ND],label=data_label);
    
    # Y limits definition
    _lim = int(np.max(num));
    positions = (np.linspace(0,1,ND+1)*ND/4)+(_lim+1);
    y_ticks = np.linspace(0,_lim,nyticks+1);

    # Error bars
#    y_error = np.sqrt(num[0])/num[0];
#    x_mid = (bin[1:]+bin[:-1])/2.;
#    ax.errorbar(x_mid,num[0],yerr=y_error);
    
    # Boxplot
    bp = ax.boxplot(data,vert=0,notch=1,positions=positions[:-1],sym='');

    ax.legend(loc='upper right')

    # Fix the colors for each boxplot plotted..
    xmin = 0;
    xmax = 0;
    for i in range(ND):
        plt.setp(bp['boxes'][i],color=colors[i],linewidth=2);
        plt.setp(bp['medians'][i],color=colors[i]);
        plt.setp(bp['caps'][2*i:2*(i+1)],color=colors[i],linewidth=2);
        plt.setp(bp['whiskers'][2*i:2*(i+1)],color=colors[i],linestyle=':');
        xmin = min(xmin,plt.getp(bp['whiskers'][2*i:2*(i+1)]),'xdata');
        xmax = max(xmax,plt.getp(bp['whiskers'][2*i:2*(i+1)]),'xdata');
        
    ax.set_ylim(0,(positions[-1]+ND*0.5));
    ax.set_yticks(y_ticks)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xmin,xmax)
    ax.set_xticks(bin)
    ax.set_xlabel(xlabel)
    
    f.suptitle(title)
    
    plt.savefig('Ell_ArcsAll.png')
    plt.show();
    
