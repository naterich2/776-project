import shlex, subprocess, os.path
import pandas as pd
import numpy as np
class gsea(object):

    """Docstring for gsea. """

    def __init__(self, data, metadata):
        """Instantiate gsea object"""
        self._data = pd.read_csv(data)
        self._meta = pd.read_csv(metadata)

    def write_gct(self,gctfile):
        """Write a GCT file from a dataframe

        Parameters
        ----------
        outfile : file to write to

        Returns
        -------
        n/a

        """
        dataframe = self._data.copy().T
        dataframe.index.name = 'NAME'
        dataframe.insert(0,column='Description',value='')
        with open(gctfile,'w') as out:
            print('#1.2',file=out)
            print(dataframe.shape[0],(dataframe.shape[1]-1),file=out)
        dataframe.to_csv(gctfile,sep='\t',mode='a',header=True)

    def write_cls(self,clsfile):
        """TODO: Docstring for write_cls.

        Parameters
        ----------
        meta : TODO
        outfile : TODO

        Returns
        -------
        n/a

        """
        meta = self._meta.copy()
        diagnosis = meta.loc[self._data.index,'diagnosis']
        with open(clsfile,'w') as out:
            print(diagnosis.shape[0], np.unique(diagnosis).shape[0], '1',file=out)
            print('#',*list(np.unique(diagnosis)),file=out)
        diagnosis.to_csv(clsfile,mode='a',index=False,header=False,line_terminator='\t')


    def run_gsea(self,
            gctfile,clsfile,label,
            outdir='/home/nrichman/gsea_home/output',
            gsea_path="/home/nrichman/GSEA_Linux_4.0.3/gsea-cli.sh",
            gmx='ftp.broadinstitute.org://pub/gsea/gene_sets/c2.all.v7.1.symbols.gmt',
            metric='Signal2Noise'):
        """TODO: Docstring for run_gsea.

        Parameters
        ----------
        gctfile : relative path to gct file
        clsfile : relative path to cls file
        label : name of analysis
        outdir : absolute path to out directory
        gsea_path : absolute path to gsea-cli.sh file
        gmx : gene database to use
        metric : metric to use

        Returns
        -------
        TODO

        """
        command = [
            '/bin/bash',
            gsea_path,'GSEA',
            '-res',os.path.abspath(gctfile),
            '-cls',os.path.abspath(clsfile),
            '-gmx',gmx,
            '-collapse','Collapse',
            '-mode','Max_probe',
            '-norm','meandiv',
            '-nperm','1000',
            '-permute','phenotype',
            '-rnd_type','no_balance',
            '-scoring_scheme','weighted',
            '-rpt_label',label,
            '-metric',metric,
            '-sort','real',
            '-order','descending',
            '-chip','ftp.broadinstitute.org://pub/gsea/annotations_versioned/Human_ENSEMBL_Gene_MSigDB.v7.1.chip',
            '-create_gcts','false',
            '-create_svgs','false',
            '-include_only_symbols','true',
            '-make_sets','true',
            '-median','false',
            '-num','100',
            '-plot_top_x','20',
            '-rnd_seed','timestamp',
            '-save_rnd_lists','false',
            '-set_max','500',
            '-set_min','15',
            '-zip_report','false',
            '-out',outdir
        ]
        print('Creating GCT file')
        self.write_gct(gctfile)
        print('Creating CLS file')
        self.write_cls(clsfile)
        print('Starting GSEA')

        gsea_proc = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
        for line in gsea_proc.stdout:
            print(line.decode('utf-8').replace('\n',''))


my_gse = gsea('data/total_coding.csv','data/total_meta.csv')
my_gse.run_gsea(gctfile='total_coding.gct',clsfile='total_meta.cls',label='test1',outdir='/home/nrichman/gsea_home/output/test1',gmx='ftp.broadinstitute.org://pub/gsea/gene_sets/c2.cgp.v7.1.symbols.gmt')

