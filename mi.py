#!/usr/bin/env python
import numpy as np
import pandas as pd
import math

class MI(object):

    """Calculates single and multivariate MI for multiple genes/ages/phenotype"""

    def __init__(self,data,meta,pseudocount=0.1,num_bins=10):
        """Initialize class with total data and metadata

        Parameters
        ----------
        data : pandas.DataFrame
            total data in the form of (samples,genes)
        meta : pandas.DataFrame
            total metadata in the form of (samples,meta_feature), needs to have 'age','diagnosis' as part of it
        pseudocount : float
            Default: 0.1
            Pseudocount to add to each cell in the joint histogram
        num_bins : int
            Default: 10
            Number of bins to use for the gene expression data

        """
        if data.shape[0] != meta.shape[0]:
            raise ValueError("Number of rows not equal, please make sure data and metadata have same rows")
            return
        if meta.shape[1] > 3 or meta.shape[1] < 2:
            raise ValueError("Metadata needs to have between 2 and 3 features.")
            return
        self._data = data
        self._meta = meta.loc[data.index,:]
        self._pseudo = pseudocount
        self._gene_bins = num_bins
        # Create bins for age and diagnosis
        age_min = np.floor(np.min(self._meta['age']) / 10.0) * 10 # Round down to nearest multiple of 10
        age_max = np.ceil(np.max(self._meta['age']) / 10.0) * 10 # Round up to nearest multiple of 10

        self._age_bins = np.round(np.linspace(age_min,
            age_max,
            math.ceil((age_max-age_min + 1)/10.0))) # Bins for approximately 1 bin every 10 years
        self._diag_bins = np.array([0,1,2]) # Diagnosis is categorical with 0 being control, and 1 being positive

    def _calc_mi(self,gene,how='all'):
        """Calculate MI for a specific gene with other features

        Parameters
        ----------
        how : String
            Default: 'all'
            Options: 'all', 'age', 'diagnosis'
            Name of column to compute mutual information with gene against.

        Returns
        -------
        KWII: k-wise Ineraction Information value
            OR
        MI: Mutual Information between two variables.

        """

        gene_min = np.min(self._data.loc[:,gene])
        gene_max = np.max(self._data.loc[:,gene])
        gene_bins = np.linspace(gene_min,gene_max,self._gene_bins + 1)

        if how == 'all' or how == 'TCI':
            joint_df = self._meta.loc[:,['age','diagnosis']]
            joint_df[gene] = self._data.loc[:,gene]
            joint_counts,bins = np.histogramdd(joint_df.values,
                    bins=(self._age_bins,self._diag_bins,gene_bins))

            # Find counts for C(age,diag,gene), C(diag,gene) and C(age,gene)
            joint_counts += self._pseudo
            diag_gene_counts = np.sum(joint_counts,axis=0) #Sum out age
            age_gene_counts = np.sum(joint_counts,axis=1) #Sum out diagnosis
            gene_counts = np.sum(np.sum(joint_counts,axis=0),axis=0)

            # Find dist for P(age,diag,gene), P(diag,gene), and P(age,gene)
            joint_dist = joint_counts/np.sum(joint_counts)
            diag_gene_dist = diag_gene_counts/np.sum(diag_gene_counts)
            age_gene_dist = age_gene_counts/np.sum(age_gene_counts)
            gene_dist = gene_counts/np.sum(gene_counts)

            # Find entropy: H(age,diag,gene), H(diag,gene), and H(age,gene)
            H_joint = - np.sum(joint_dist*np.log2(joint_dist))
            H_diag_gene = - np.sum(diag_gene_dist*np.log2(diag_gene_dist))
            H_age_gene = - np.sum(age_gene_dist*np.log2(age_gene_dist))
            H_gene = - np.sum(gene_dist*np.log2(gene_dist))

            # Calculate K-way Interaction Information:
            #   KWII(age,diag,gene) =
            #       - H(age) - H(diag) - H(gene) + H(diag,gene) + H(age,gene) - H(age,diag,gene)
            if how == 'all':
                KWII =  H_gene + self._H_diag + self._H_age - \
                        H_diag_gene - H_age_gene - self._H_diag_age + H_joint
                return KWII
            else:
                TCI = H_gene + self._H_diag + self._H_age - H_joint
                return TCI

        elif how == 'age':
            joint_df = self._meta.loc[:,'age'].to_frame()
            joint_df[gene] = self._data.loc[:,gene]
            joint_counts,bins = np.histogramdd(joint_df.loc[:,['age',gene]].values,
                    bins=(self._age_bins,gene_bins))

            joint_counts += self._pseudo
            gene_counts = np.sum(joint_counts,axis=0)

            joint_dist = joint_counts/np.sum(joint_counts)
            gene_dist = gene_counts/np.sum(gene_counts)

            H_joint = - np.sum(joint_dist*np.log2(joint_dist))
            H_gene = - np.sum(gene_dist*np.log2(gene_dist))

            MI =  H_gene + self._H_age - H_joint
            return MI

        elif how == 'diagnosis':
            joint_df = self._meta.loc[:,'diagnosis'].to_frame()
            joint_df[gene] = self._data.loc[:,gene]
            joint_counts,bins = np.histogramdd(joint_df.loc[:,['diagnosis',gene]].values,
                    bins=(self._diag_bins,gene_bins))

            joint_counts += self._pseudo
            gene_counts = np.sum(joint_counts,axis=0)

            joint_dist = joint_counts/np.sum(joint_counts)
            gene_dist = gene_counts/np.sum(gene_counts)

            H_joint = - np.sum(joint_dist*np.log2(joint_dist))
            H_gene = - np.sum(gene_dist*np.log2(gene_dist))

            MI =  H_gene + self._H_diag - H_joint
            return MI

        else:
            raise ValueError("How needs to be one of {'all','age','diagnosis'}")

    def run(self,how='all'):
        """Run MI for each gene

        Parameters
        ----------
        how : String
            Default: 'all'
            Options: 'all', 'age', 'diagnosis','TCI'
            Name of column to compute mutual information with gene against.

        Returns
        -------
        TODO

        """
        total_counts_size = self._gene_bins*(self._age_bins.shape[0] - 1)*(self._diag_bins.shape[0] - 1 )
        self._obs = self._data.shape[0] + self._pseudo*total_counts_size

        age_counts,bins = np.histogram(self._meta['age'],bins=self._age_bins)
        diag_counts,bins = np.histogram(self._meta['diagnosis'],bins=self._diag_bins)
        diag_age_counts,bins = np.histogramdd(self._meta.loc[:,['age','diagnosis']].values,
                bins=(self._age_bins,self._diag_bins))

        age_counts = age_counts + (total_counts_size/(self._age_bins.shape[0] - 1)) * self._pseudo
        diag_counts = diag_counts + (total_counts_size/(self._diag_bins.shape[0] - 1)) * self._pseudo
        diag_age_counts = diag_age_counts + \
                (total_counts_size/((self._diag_bins.shape[0] - 1)*(self._age_bins.shape[0] - 1))) * self._pseudo

        self._age_dist = age_counts/np.sum(age_counts)
        self._diag_dist = diag_counts/np.sum(diag_counts)
        diag_age_dist = diag_age_counts/np.sum(diag_age_counts)

        self._H_age = - np.sum(self._age_dist*np.log2(self._age_dist))
        self._H_diag = - np.sum(self._diag_dist*np.log2(self._diag_dist))
        self._H_diag_age = - np.sum(diag_age_dist*np.log2(diag_age_dist))

        MIs = {}
        for gene in self._data.columns:
            MI = self._calc_mi(gene,how)
            MIs[gene] = MI
        MIs_df = pd.DataFrame.from_dict(MIs,orient='index')
        MIs_df.columns = ['MI']
        return MIs_df

