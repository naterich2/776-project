#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# My files
import preprocessing
import mi

## Preprocess
#total,meta = preprocessing.preprocess()

## MI
#my_mi = mi.MI(total,meta)

plt.switch_backend('Qt5Agg')
FDR = 0.5
n_iter = 100
#my_mi.run(n_iter)
#p_vals_all.to_csv('data/p_vals_all.csv',index_label=False)
#p_vals_diag = my_mi.run(n_iter,how='diag')
#p_vals_diag.to_csv('data/p_vals_diag.csv',index_label=False)
dist_diag = pd.read_csv('data/dist.csv',index_col=0)
dist_all = pd.read_csv('data/dist_all.csv',index_col=0)


diag_means = np.mean(dist_diag.iloc[:,1:],axis=1)
all_means = np.mean(dist_all.iloc[:,1:],axis=1)

diag_stds = np.std(dist_diag.iloc[:,1:],axis=1)
all_stds = np.std(dist_all.iloc[:,1:],axis=1)

diag_t = (dist_diag.iloc[:,0] - diag_means)/diag_stds
all_t = (dist_all.iloc[:,0] - all_means)/all_stds

diag_p = stats.t.sf(diag_t,df=99)
all_p = stats.t.sf(all_t,df=99)

p_diag = pd.DataFrame(data=diag_p,index=dist_diag.index,columns=['p_val'])
p_all = pd.DataFrame(data=all_p,index=dist_all.index,columns=['p_val'])

plt.style.use('ggplot')
fig,[ax1,ax2] = plt.subplots(2,1,figsize=(10,15))

#ax1.hist(p_diag.values,density=True,bins=np.linspace(0,1))
#ax1.set_xlabel('p-val')
#ax1.set_ylabel('Frequency')
#ax1.set_title('MI by Diagnosis')
#ax2.hist(p_all.values,density=True,bins=np.linspace(0,1))
#ax2.set_title('KWII')
#ax2.set_xlabel('p-val')
#ax2.set_ylabel('Frequency')
#plt.show()

plt.hist(dist_diag.loc['ENSG00000160963',:])
plt.show()
#fig.savefig('output.png')
##Diffusion map/clustering/pseudotime

##GSEA




