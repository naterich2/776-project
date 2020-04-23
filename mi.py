#!/usr/bin/env python

class MI(object):

    """Calculates single and multivariate MI for multiple genes/ages/phenotype"""

    def __init__(self,data,meta,pseudocount=0.1):
        """Initialize class with total data and metadata

        Parameters
        ----------
        data : pandas.DataFrame
            total data in the form of (samples,genes)
        meta : pandas.DataFrame
            total metadata in the form of (samples,meta_feature), needs to have 'age','diagnosis' as part of it
        pseudocount: float
            Default: 0.1
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
        # Create bins for age and diagnosis
        age_min = np.floor(np.min(self._meta['age']) / 10) * 10 # Round down to nearest multiple of 10
        age_max = np.ceil(np.max(self._meta['age']) / 10) * 10 # Round up to nearest multiple of 10

        self._age_bins = np.round(np.linspace(age_min,age_max,10)) # 10 bins rounded to nearest age
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


        """

        gene_min = np.min(self._data.loc[:,gene])
        gene_max = np.max(self._data.loc[:,gene])
        gene_bins = np.linspace(gene_min,gene_max,10)
        total_counts_size = gene_bins.shape[0]*self._age_bins.shape[0]*self._diag_bins.shape[0] # TODO move to above function
        age_counts = np.histogram(self._meta['age'])

        if how == 'all':
            joint_df = self._meta.loc[:,['age','diagnosis']]
            joint_df[gene] = self._data.loc[:,gene]
            joint_counts,bins = np.histogramdd(joint_df.values,
                    bins(self._age_bins,self._diag_bins,gene_bins))

            joint_counts += self._pseudo

            age_counts = np.sum(np.sum(joint_counts,axis=1),axis=1)
            diag_counts = np.sum(np.sum(join_counts,axis=1))


        elif how == 'age':
            joint_df = self._meta.loc[:,'age']
            joint_df[gene] = self._data.loc[:,gene]
            joint_counts,bins = np.histogramdd(joint_df.values,
                    bins(self._age_bins,gene_bins))
            joint_counts += self._pseudo

        elif how == 'diagnosis':
            joint_df = self._meta.loc[:,'diagnosis']
            joint_df[gene] = self._data.loc[:,gene]
            joint_counts,bins = np.histogramdd(joint_df.values,
                    bins(self._diag_bins,gene_bins))
            joint_counts += self._pseudo

        else:
            raise ValueError("How needs to be one of {'all','age','diagnosis'}")
