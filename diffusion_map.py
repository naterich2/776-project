import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import subprocess

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

from plot_embeddings import plot_embedding

plt.switch_backend('Agg')
# Import total data
data = pd.read_csv('data/total_data.csv',index_col=0)
meta = pd.read_csv('data/total_meta.csv',index_col=0)
depr = pd.read_csv('data/depr_row.csv',index_col=0)

# Import p-values from permutations
p_all = pd.read_csv('data/p_vals_all.csv',index_col=0)
p_diag = pd.read_csv('data/p_vals_diag.csv',index_col=0)

# Select significant genes from KWII and export
p_all_sig_genes = p_all[p_all['p'] < 0.01]
p_all_sig = data.loc[:,p_all_sig_genes.index]
p_all_sig.to_csv('data/p_all_sig.csv',index_label=False)

# Select significant genes from diagnosis+expression MI and export
p_diag_sig_genes = p_diag[p_diag['p'] <= p_diag['bounds']]
p_diag_sig = data.loc[:,p_diag_sig_genes.index]
p_diag_sig.to_csv('data/p_diag_sig.csv',index_label=False)

#
p_depr = data.loc[:,data.columns.intersection(depr['name'])]
p_depr.to_csv('data/p_depr.csv')
# Run Rscript
rscript = subprocess.Popen(['Rscript','diffusion_map.R'])
rscript.wait()

# Import data
p_all_dm = pd.read_csv('data/p_all_sig_dm.csv',index_col=0)
p_diag_dm = pd.read_csv('data/p_diag_sig_dm.csv',index_col=0)
p_depr_dm = pd.read_csv('data/p_depr_dm.csv',index_col=0)

## Plot data
tsne_diag = TSNE(n_components=3)
tsne_all = TSNE(n_components=3)
tsne_depr = TSNE(n_components=3)

pca_diag = PCA(n_components=3)
pca_all = PCA(n_components=3)
pca_depr = PCA(n_components=3)

tsne_diag_points = tsne_diag.fit_transform(p_diag_sig)
tsne_all_points = tsne_all.fit_transform(p_all_sig)
tsne_depr = tsne_all.fit_transform(p_depr_dm)

plot_embedding(tsne_diag_points,meta,'TSNE embeddings (expression,diagnosis)',
        outfile='figs/tsne_diag.png',labels='TSNE')
plot_embedding(tsne_all_points,meta,'TSNE embeddings (expression,age,diagnosis)',
        outfile='figs/tsne_all.png',labels='TSNE')
plot_embedding(tsne_depr,meta,'TSNE embeddings depression genes',
        outfile='figs/tsne_depr.png',labels='TSNE')

pca_diag_points = pca_diag.fit_transform(p_diag_sig)
pca_all_points = pca_all.fit_transform(p_all_sig)
pca_depr_points = pca_depr.fit_transform(p_depr)

plot_embedding(pca_diag_points,meta,'PCA scores (expression,diagnosis)',
        outfile='figs/pca_diag.png',labels='PC')
plot_embedding(pca_all_points,meta,'PCA scores (expression,age,diagnosis)',
        outfile='figs/pca_all.png',labels='PC')
plot_embedding(pca_depr_points,meta,'PCA scores depression genes',
        outfile='figs/pca_depr.png',labels='PC')

plot_embedding(p_diag_dm.values,meta,'Diffusionmap embeddings (expression,diagnosis)',
        outfile='figs/dm_diag.png',labels='DC',by_sample=False)
plot_embedding(p_all_dm.values,meta,'Diffusionmap embeddings (expression,age,diagnosis)',
        outfile='figs/dm_all.png',labels='DC',by_sample=False)
plot_embedding(p_depr.values,meta,'Diffusionmap embeddings depression genes',
        outfile='figs/dm_depr.png',labels='DC',by_sample=False)





