#!/usr/bin/env python

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# My files
import preprocessing
import mi

## Preprocess
total,meta = preprocessing.preprocess()

## MI
my_mi = mi.MI(total,meta)

FDR = 0.5
n_iter = 100
p_vals_all = mi.run(n_iter)
p_vals_all.to_csv('data/p_vals_all.csv',index_label=False)
p_vals_diag = mi.run(n_iter,how='diag')
p_vals_diag.to_csv('data/p_vals_diag.csv',index_label=False)

##Diffusion map/clustering/pseudotime

##GSEA




