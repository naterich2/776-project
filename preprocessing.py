import pandas as pd
import numpy as np

from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler
from sklearn.preprocessing import MinMaxScaler

def preprocess():
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

    # Read gene annotations and only keep coding genes
    genes = pd.read_csv('data/gene_annotation.csv',index_col=0)
    genes_intersection = genes[genes['gene_type'] == 'protein_coding'].index.intersection(total.columns)

    total_coding = total.loc[:,genes_intersection]

    age_info['colors'] = age_info['diagnosis'].replace('M','0').replace('C','1')
    age_info['diag_cat'] = 0
    age_info.loc[age_info['diagnosis'] == 'M','diag_cat'] = 1

    ### Mutual Information
    meta = age_info.loc[total.index,['age','diag_cat']]
    meta.columns = ['age','diagnosis']
    print('Finished preprocessing, {} samples of {} genes'.format(*total_coding.shape))
    return (total_coding,meta)
