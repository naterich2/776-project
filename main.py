#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


##### Preprocessing

gse101521 = pd.read_csv('data/GSE101521_totalRNA_counts.csv',index_col=0)
gse101521_meta = pd.read_csv('data/GSE101521_samples_metadata.csv')

# Find indices of control or MDD
control_ind = np.where(gse101521_meta['diagnosis'] == 'CON')[0]
mdd_ind = np.where((gse101521_meta['diagnosis'] == 'MDD') | (gse101521_meta['diagnosis'] == 'MDD-S'))[0]

# metadata file contains more samples than full data, so only take metadata thats in full data file
control_ind = control_ind[control_ind < len(gse101521.columns)]
mdd_ind = mdd_ind[mdd_ind < len(gse101521.columns)]

gse101521_c = gse101521.iloc[:,control_ind]
gse101521_m = gse101521.iloc[:,mdd_ind]


## GSE80655
gse80655 = pd.read_csv('data/GSE80655_GeneExpressionData_Updated_3-26-2018.csv',sep='\t',header=0,index_col=0)
gse80655_meta = pd.read_csv('data/GSE80655_meta.csv',index_col=0)

# This dataset also contains bipolar and schizophrenia so only select mdd
gse80655c_or_m_meta = gse80655_meta[(gse80655_meta['diagnosis'] == 'C') | (gse80655_meta['diagnosis'] == 'M')]
gse80655c_or_m = gse80655.loc[:,gse80655c_or_m_meta.index]
gse80655c_or_m.to_csv('data/GSE80655_data.csv',index_label=False)

# Data needs to be normalized, used with DeSeq package in R.:wq
## GSE42546
gse42546 = pd.read_csv('data/GSE42546_CleanedRawCounts.csv',sep='\t',header=0,index_col=0)
gse42546_meta = pd.read_csv('data/GSE42546_sample_metadata.csv',index_col=0)

# This dataset also contains bipolar and schizophrenia so only select mdd
gse42546c_or_m_meta = gse42546_meta[(gse42546_meta['diagnosis'] == 'control') | (gse42546_meta['diagnosis'] == 'depression')]
gse42546c_or_m_meta.to_csv('data/GSE42546_meta.csv',index_label=False)
gse42546c_or_m = gse42546.loc[:,gse42546c_or_m_meta.index]
gse42546c_or_m.to_csv('data/GSE42546_data.csv',index_label=False)



