library('destiny')

setwd("~/project")

all_sig <- read.csv('data/p_all_sig.csv')
diag_sig <- read.csv('data/p_diag_sig.csv')
depr_sig <- read.csv('data/p_depr.csv')


diffmap_all <- DiffusionMap(data = all_sig)
diffmap_diag <- DiffusionMap(data = diag_sig)
diffmap_depr <- DiffusionMap(data = depr_sig)

eigs_all <- as.data.frame(eigenvectors(diffmap_all))
eigs_diag <- as.data.frame(eigenvectors(diffmap_diag))
eigs_depr <- as.data.frame(eigenvectors(diffmap_depr))


write.csv(eigs_all,file='data/p_all_sig_dm.csv')
write.csv(eigs_diag,file='data/p_diag_sig_dm.csv')
write.csv(eigs_depr,file='data/p_depr_dm.csv')


