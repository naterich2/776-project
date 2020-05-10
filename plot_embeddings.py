
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as Normalize
import matplotlib.gridspec as gridspec
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.lines import Line2D

def plot_embedding(data,meta,figtitle,outfile='out.png',labels='DC'):
    """Create 3D plot for embedding

    Parameters
    ----------
    data : data to plot
    meta : metadata
    figtitle : title of figure
    outfile : file to save image in
    labels : axes label prefix, ex. 'DC' leads to axes labels of 'DC1','DC2','DC3'

    Returns
    -------
    n/a
    """
    gs = gridspec.GridSpec(2,2)
    fig = plt.figure(figsize=(15,15))
    gs = fig.add_gridspec(nrows=2,ncols=2)
    ax1 = fig.add_subplot(gs[0,0],projection='3d')
    ages = meta['age']
    ax1.scatter(data[:,0],data[:,1],data[:,2],c=meta['sample'])
    ax1.set_xlabel(labels+'1')
    ax1.set_ylabel(labels+'2')
    ax1.set_zlabel(labels+'3')
    ax1.set_title('Colored by dataset')

    ax2 = fig.add_subplot(gs[0,1],projection='3d')
    ages = meta['age']
    ax2.scatter(data[:,0],data[:,1],data[:,2],c=meta['diagnosis'].str.lower())
    ax2.set_xlabel(labels+'1')
    ax2.set_ylabel(labels+'2')
    ax2.set_zlabel(labels+'3')
    ax2.set_title('Colored by diagnosis')
    control=Line2D([0],[0],color='c',marker='o',markerfacecolor='c',
            markersize=6,lw=0,label='Control')
    mdd=Line2D([0],[0],lw=0,color='m',marker='o',markerfacecolor='m',
            markersize=6,label='MDD')
    ax2.legend(handles=[control,mdd])

    gs01 = gs[1,:].subgridspec(1,5)
    ax3 = fig.add_subplot(gs01[1:4],projection='3d')
    ages = meta['age']
    age_scatter = ax3.scatter(data[:,0],data[:,1],data[:,2],
            c=ages,norm=Normalize(np.min(ages),np.max(ages)))
    ax3.set_xlabel(labels+'1')
    ax3.set_ylabel(labels+'2')
    ax3.set_zlabel(labels+'3')
    ax3.set_title('Colored by age')

    fig.colorbar(age_scatter,ax=ax3,shrink=0.7)
    fig.suptitle(figtitle,fontsize=16)

    fig.tight_layout()
    fig.savefig(outfile)
