#!/usr/bin/env python

import math
import requests
import mi

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize as Normalize
from mpl_toolkits.mplot3d import Axes3D

from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram
from scipy.stats import norm
from scipy.stats import beta

from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler
from sklearn.decomposition import PCA

import pydiffmap as dm
from pydiffmap.visualization import embedding_plot, data_plot


### Load datasets
gse101521 = pd.read_csv('data/GSE101521_totalRNA_counts.csv',index_col=0)
gse101521_meta = pd.read_csv('data/GSE101521_samples_metadata.csv',index_col=0)
gse101521_meta['diagnosis'] = gse101521_meta['diagnosis'].replace(r'CON','C',regex=True)
gse101521_meta['diagnosis'] = gse101521_meta['diagnosis'].replace(r'MDD-S','M',regex=True)
gse101521_meta['diagnosis'] = gse101521_meta['diagnosis'].replace(r'MDD','M',regex=True)
gse101521_meta.index = gse101521_meta['name']

# Find indices of control or MDD
gse101521_c_ids = pd.Index(gse101521_meta[gse101521_meta['diagnosis'] == 'C']['name']).intersection(gse101521.columns)
gse101521_m_ids = pd.Index(gse101521_meta[gse101521_meta['diagnosis'] == 'M']['name']).intersection(gse101521.columns)

#GSE80655 dataset
gse80655 = pd.read_csv('data/GSE80655_data_norm.csv',index_col=0)
gse80655_meta = pd.read_csv('data/GSE80655_meta.csv')

gse80655 = gse80655.loc[:,gse80655_meta[gse80655_meta['brain_region'] == 'DLPFC'].index.intersection(gse80655.columns)]

gse80655_c_ids = gse80655_meta[gse80655_meta['diagnosis'] == 'C'].index.intersection(gse80655.columns)
gse80655_m_ids = gse80655_meta[gse80655_meta['diagnosis'] == 'M'].index.intersection(gse80655.columns)

### Normalizatoin
# Normalize
ss = StandardScaler()
rs = RobustScaler()
ms = MinMaxScaler()
gse80655_scaled = ss.fit_transform(gse80655)
gse80655 = pd.DataFrame(data=gse80655_scaled,
        index=gse80655.index,
        columns=gse80655.columns)
gse80655 = gse80655.T
gse101521_scaled = ss.fit_transform(gse101521)
gse101521 = pd.DataFrame(data=gse101521_scaled,
        index=gse101521.index,
        columns=gse101521.columns)
gse101521 = gse101521.T
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

### Mutual Information
meta = age_info.loc[total.index,['age','diag_cat']]
meta.columns = ['age','diagnosis']
# Run MI algorithm on all different types
my_mi = mi.MI(total,meta)
mis_age = my_mi.run('age')
mis_diag = my_mi.run('diagnosis')
mis_kwii = my_mi.run()
mis_tci = my_mi.run('TCI')
#relabel columns
mis_age.columns = ['by_age']
mis_diag.columns = ['by_diagnosis']
mis_kwii.columns = ['by_kwii']
mis_tci.columns = ['by_tci']
total_mis = mis_age.join(mis_diag).join(mis_kwii).join(mis_tci)
total_mis.to_csv('data/total_values.csv',index_label=False)

top_100_age = total.loc[:,total_mis.sort_values('by_age').iloc[-100:,].index]
top_100_diag = total.loc[:,total_mis.sort_values('by_diagnosis').iloc[-100:,].index]
top_100_kwii = total.loc[:,total_mis.sort_values('by_kwii').iloc[-100:,].index]
top_100_tci = total.loc[:,total_mis.sort_values('by_tci').iloc[-100:,].index]

top_100_age.to_csv('data/top_100_age.csv',index_label=False)
top_100_diag.to_csv('data/top_100_diag.csv',index_label=False)
top_100_kwii.to_csv('data/top_100_kwii.csv',index_label=False)
top_100_tci.to_csv('data/top_100_tci.csv',index_label=False)


### DIFFUSIONMAP/TSNE
diffusion_age = pd.read_csv('data/top_100_age_dm.csv',index_col=0)
diffusion_diag = pd.read_csv('data/top_100_diag_dm.csv',index_col=0)
diffusion_kwii = pd.read_csv('data/top_100_kwii_dm.csv',index_col=0)
diffusion_tci = pd.read_csv('data/top_100_tci_dm.csv',index_col=0)
diffusion_all = pd.read_csv('data/p_all_sig_dm.csv',index_col=0)
diffusions = (diffusion_age,diffusion_diag,diffusion_kwii,diffusion_tci)
names = ['By Age','By Diagnosis','By KWII','By TCI']
plt.style.use('ggplot')
fig,ax = plt.subplots(2,2,figsize=(10,10))
for i,axes in enumerate(ax.ravel()):
    ages = age_info.loc[diffusions[i].index,'age_binned']
    axes.scatter(diffusions[i].iloc[:,0],diffusions[i].iloc[:,1],
            c=ages,
            norm=Normalize(np.min(ages),np.max(ages)))
    axes.set_xlabel('DC1')
    axes.set_ylabel('DC2')
    axes.set_title(names[i])
fig,ax = plt.subplots(1,1)
total_meta = pd.read_csv('data/total_meta.csv')
ages = total_meta['age']
ax.scatter(diffusion_all.iloc[:,0],diffusion_all.iloc[:,1],
        c=total_meta['diagnosis'])
#        norm=Normalize(np.min(ages),np.max(ages)))
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

## Mutual information analysis
p_diag = pd.read_csv('data/p_vals_diag.csv',index_col=0)
p_all = pd.read_csv('data/p_vals_all.csv',index_col=0)
plt.style.use('ggplot')

x = np.arange(1,p_all.shape[0]+1)
x_cont = np.linspace(0,1)
fig,[[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2,figsize=(15,15))
ax1.scatter(x,p_diag['MI'])
ax1.set_title('MI values')
ax2.scatter(x,p_diag['means'])
ax2.set_title('MI permutation means')
ax3.scatter(x,p_diag['stds'])
ax3.set_title('MI permutation standard deviations')
ax4.plot(x_cont,beta.pdf(x_cont,alpha,beta))
ax4.plot(x_cont,norm.pdf(x_cont,p_mean,np.sqrt(p_var)))
ax4.hist(p_diag['p'],density=True)
ax4.legend(labels=['Beta distribution','Normal distribution'])
ax4.set_title('Distribution of p-values')
fig.savefig('out.png')



fig,[[ax1,ax2],[ax3,ax4]] = plt.subplots(2,2,figsize=(15,15))
p_mean = np.mean(p_all['p'])
p_var = np.var(p_all['p'])
ax1.scatter(x,p_all['MI'])
ax1.set_title('MI values')
ax2.scatter(x,p_all['means'])
ax2.set_title('MI permutation means')
ax3.scatter(x,p_all['stds'])
ax3.set_title('MI permutation standard deviations')
#ax4.plot(x_cont,beta.pdf(x_cont,alpha,beta))
ax4.plot(x_cont,norm.pdf(x_cont,p_mean,np.sqrt(p_var)))
ax4.hist(p_all['p'],density=True)
ax4.legend(labels=['Beta distribution','Normal distribution'])
ax4.set_title('Distribution of p-values')
plt.plot(x,p_all['p'])




