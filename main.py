#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats

# My files
import preprocessing
import mi

## Preprocess
total,meta = preprocessing.preprocess()

## MI
my_mi = mi.MI(total,meta)

plt.switch_backend('Qt5Agg')
FDR = 0.5
n_iter = 100
#my_mi.run(n_iter)
p_vals_all.to_csv('data/p_vals_all.csv',index_label=False)
p_vals_diag = my_mi.run(n_iter,how='diag')
p_vals_diag.to_csv('data/p_vals_diag.csv',index_label=False)
##Diffusion map/clustering/pseudotime


##GSEA




