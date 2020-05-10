import pandas as pd

from gsea_utils import gsea

# Import data
p_all_sig = pd.read_csv('data/p_all_sig.csv',index_col=0)
p_diag_sig = pd.read_csv('data/p_diag_sig.csv',index_col=0)

#Import metadata and gene annotations
meta = pd.read_csv('data/total_meta.csv',index_col=0)
annotations = pd.read_csv('data/gene_annotation.csv',index_col=0)

# Create GSEA objects and run
#all_gsea = gsea('data/p_all_sig.csv','data/total_meta.csv')
#all_gsea.run_gsea(gctfile='p_all_sig.gct',clsfile='total_meta.cls',label='p_all_sig',outdir='/home/nrichman/gsea_home/output/project',gmx='ftp.broadinstitute.org://pub/gsea/gene_sets/c2.cgp.v7.1.symbols.gmt')

#diag_gsea = gsea('data/p_diag_sig.csv','data/total_meta.csv')
#diag_gsea.run_gsea(gctfile='p_diag_sig.gct',clsfile='total_meta.cls',label='p_diag_sig',outdir='/home/nrichman/gsea_home/output/project',gmx='ftp.broadinstitute.org://pub/gsea/gene_sets/c2.cgp.v7.1.symbols.gmt')
#depr_genes = pd.read_csv('data/depr_row.csv',index_col=0)

#gsea.create_gmt([p_all_sig.columns,p_diag_sig.columns,depr_genes['name']],
#        annotations='data/gene_annotation.csv',
#        names=['all_sig','diag_sig','depr_genes'],
#        descriptions=['n/a','n/a','n/a'],outfile='data/gene_sets.gmt')
coding_gsea = gsea('data/total_coding.csv','data/total_meta.csv')

coding_gsea.run_gsea(gctfile='total_coding.gct',clsfile='total_meta.cls',label='total_coding',outdir='/home/nrichman/gsea_home/output/project',gmx='data/gene_sets.gmt')

