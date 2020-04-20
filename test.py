#!/usr/bin/env python

import math
import requests

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import dendrogram

from sklearn.manifold import TSNE


total_data = pd.read_csv('data/GSE101521_totalRNA_counts.csv',index_col=0)
total_metadata = pd.read_csv('data/GSE101521_samples_metadata.csv')

control_ind = np.where(total_metadata['diagnosis'] == 'CON')[0]
control_ind = control_ind[control_ind < len(total_data.columns)]

mdd_ind = np.where((total_metadata['diagnosis'] == 'MDD') | (total_metadata['diagnosis'] == 'MDD-S'))[0]
mdd_ind = mdd_ind[mdd_ind < len(total_data.columns)]

control = total_data.iloc[:,control_ind]
mdd = total_data.iloc[:,mdd_ind]

host = "https://rest.ensembl.org"
endpoint = "/lookup/id?expand=0&format=condensed&species=homo_sapiens"
headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
for page in math.ceil(total_data.index.shape[0] / 1000):





