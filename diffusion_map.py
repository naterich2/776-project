import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import subprocess
# Import total data
data = pd.read_csv('data/total_data.csv',index_col=0)
meta = pd.read_csv('data/total_meta.csv',index_col=0)

# Import p-values from permutations
p_all = pd.read_csv('data/p_vals_all.csv',index_col=0)
p_diag = pd.read_csv('data/p_vals_diag.csv',index_col=0)

# Select significant genes from KWII and export
p_all_sig_genes = p_all[p_all['p'] < 0.01]
p_all_sig = data.loc[:,p_all_sig_genes.index]
p_all_sig.to_csv('data/p_all_sig.csv',index_label=False)

# Select significant genes from diagnosis+expression MI and export

# Run Rscript
rscript = subprocess.Popen(['Rscript','diffusion_map.R'])
rscript.wait()

# Import data
p_all_dm = pd.read_csv('data/p_all_sig_dm.csv')
p_diag_dm = pd.read_csv('data/p_diag_sig_dm.csv')

## Plot data
fig,[ax1,ax2] = plt.subplots(1,2,figsize=(16,8),sharey=True)
ages = meta['age']
diagnosis = meta['diagnosis'].str.lower()
age_scatter = ax1.scatter(p_all_dm.iloc[:,0],p_all_dm.iloc[:,1],
        c=ages,
        norm=Normalize(np.min(ages),np.max(ages)),
        cmap='spring')
ax1.set_ylabel('DC2')
ax1.set_xlabel('DC1')
ax1.set_title('Diffusion map embedding vs. age')
fig.colorbar(age_scatter,ax=ax1)
ax2.scatter(p_all_dm.iloc[:,0],p_all_dm.iloc[:,1],
    c=diagnosis)
control=Line2D([0],[0],color='c',marker='o',markerfacecolor='c',
        markersize=6,lw=0,label='Control')
mdd=Line2D([0],[0],lw=0,color='m',marker='o',markerfacecolor='m',
        markersize=6,label='MDD')
ax2.legend(handles=[control,mdd])
ax2.set_xlabel('DC1')
ax2.set_title('Diffusion map embedding vs. diagnosis')




