#!/usr/bin/env python

import math
import requests

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as Normalize
from mpl_toolkits.mplot3d import Axes3D

from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram

from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

import pydiffmap as dm
from pydiffmap.visualization import embedding_plot, data_plot


### Load datasets
gse101521 = pd.read_csv('data/GSE101521_totalRNA_counts.csv',index_col=0).T
gse101521_meta = pd.read_csv('data/GSE101521_samples_metadata.csv',index_col=0)
gse101521_meta['diagnosis'] = gse101521_meta['diagnosis'].replace(r'CON','C',regex=True)
gse101521_meta['diagnosis'] = gse101521_meta['diagnosis'].replace(r'MDD-S','M',regex=True)
gse101521_meta['diagnosis'] = gse101521_meta['diagnosis'].replace(r'MDD','M',regex=True)
gse101521_meta.index = gse101521_meta['name']

# Find indices of control or MDD
gse101521_c_ids = pd.Index(gse101521_meta[gse101521_meta['diagnosis'] == 'C']['name']).intersection(gse101521.index)
gse101521_m_ids = pd.Index(gse101521_meta[gse101521_meta['diagnosis'] == 'M']['name']).intersection(gse101521.index)

#GSE80655 dataset
gse80655 = pd.read_csv('data/GSE80655_data_norm.csv',index_col=0).T
gse80655_meta = pd.read_csv('data/GSE80655_meta.csv')

gse80655 = gse80655.loc[gse80655_meta[gse80655_meta['brain_region'] == 'DLPFC'].index.intersection(gse80655.index),:]

gse80655_c_ids = gse80655_meta[gse80655_meta['diagnosis'] == 'C'].index.intersection(gse80655.index)
gse80655_m_ids = gse80655_meta[gse80655_meta['diagnosis'] == 'M'].index.intersection(gse80655.index)

### Normalizatoin
# Normalize
ss = StandardScaler()
gse80655 = pd.DataFrame(data=ss.fit_transform(gse80655),
        index=gse80655.index,
        columns=gse80655.columns)
gse101521 = pd.DataFrame(data=ss.fit_transform(gse101521),
        index=gse101521.index,
        columns=gse101521.columns)
gse80655_age_info = gse80655_meta[((gse80655_meta['diagnosis'] == 'C') | (gse80655_meta['diagnosis'] == 'M')) & (gse80655_meta['brain_region'] == 'DLPFC')].loc[:,('age','diagnosis')].reindex(index=gse80655.index)
gse101521_age_info = gse101521_meta[(gse101521_meta['diagnosis'] == 'C') | (gse101521_meta['diagnosis'] == 'M')].loc[:,('age','diagnosis')].reindex(index=gse101521.index)
age_info = gse80655_age_info.append(gse101521_age_info)
age_bins = np.array([10,20,30,40,50,60,70,80])
age_binned = np.digitize(age_info['age'],age_bins) - 1
age_info['age_binned'] = age_binned
age_info.to_csv('data/total_meta.csv',index_label=False)


#Find common genes in datasets
common_genes = gse80655.columns.intersection(gse101521.columns)
#Make one large control dataset
control = gse80655.reindex(index=gse80655_c_ids,columns=common_genes).append(gse101521.reindex(index=gse101521_c_ids,columns=common_genes))
mdd = gse80655.reindex(index=gse80655_m_ids,columns=common_genes).append(gse101521.reindex(index=gse101521_m_ids,columns=common_genes))
total = control.append(mdd)
total.to_csv('data/total_data.csv',index_label=False)

age_info['colors'] = age_info['diagnosis'].replace('M','0').replace('C','1')
age_info['diag_cat'] = 0
age_info.loc[age_info['diagnosis'] == 'M','diag_cat'] = 1


### DIFFUSIONMAP/TSNE
diffusion = pd.read_csv('data/diffusion_map.csv',index_col=0)
tsne = TSNE()
pca = PCA(n_components=10)
pca_trans = pca.fit_transform(control.append(mdd))
#colors = age_info.loc[total.index,'colors'].values
plt.style.use('ggplot')
points_tsne = tsne.fit_transform(pca_trans)
#diffmap = dm.diffusion_map.DiffusionMap.from_sklearn(n_evecs=3,k=100)
#points = diffmap.fit_transform(total)
#embedding_plot(diffmap)
#data_plot(diffmap, dim=3, scatter_kwargs = {'c': ages,'norm': Normalize(np.min(ages),np.max(ages))})
#fig,[ax1,ax2,ax3] = plt.subplots(3,1,figsize=(10,15))

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1,projection='3d')
ages = age_info.loc[diffusion.index,'age_binned']
diffusion_scatter = ax1.scatter(diffusion.iloc[:,0],diffusion.iloc[:,1],diffusion.iloc[:,2],
        c=ages,
        label=age_bins[ages],
        norm=Normalize(np.min(ages),np.max(ages)))
ax1.legend(handles=diffusion_scatter.legend_elements()[0],labels=('10-19','20-29','30-39','40-49','50-59','60-69','70-79','80-69'))
ax1.set_xlabel('DC1')
ax1.set_ylabel('DC2')
ax1.set_zlabel('DC3')
#fig2,ax2 = plt.subplots(1,1)
#tsne_scatter = ax2.scatter(points_tsne[:,0],points_tsne[:,1],c=ages)#,norm=Normalize(np.min(ages),np.max(ages)))
#pca_scatter = ax3.scatter(pca_trans[:,0],pca_trans[:,1],c=colors)#,norm=Normalize(np.min(ages),np.max(ages)))
fig2 = plt.figure()
hasarray=[]
ax2 = fig2.add_subplot(1,1,1,projection='3d')
colors = age_info.loc[diffusion.index,'colors']
diffusion_scatter2 = ax2.scatter(diffusion.iloc[:,0],diffusion.iloc[:,1],diffusion.iloc[:,2],
        c=colors.values)
ax2.legend(numpoints=2,labels=('Control','MDD'))
ax2.set_xlabel('DC1')
ax2.set_ylabel('DC2')
ax2.set_zlabel('DC3')

#fig.colorbar(tsne_scatter)
plt.show()

# Mutual Information
def calc_mi(gene,age,diagnosis,how='all'):
    """TODO: Docstring for calc_mi.

    Parameters
    ----------
    gene : TODO
    age : TODO
    diagnosis : TODO
    diagnosis : TODO, optional

    Returns
    -------
    TODO

    """

meta = age_info.loc[total.index,['age','diag_cat']]
diag_bins = np.array([0,1,2])
for gene in total.columns:
e   joint_df = meta
    joint_df[gene] = total.loc[:,gene]

    join = np.histogramdd()







