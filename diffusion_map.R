library('destiny')

setwd("~/project")

all_sig = read.csv('data/p_all_sig.csv')


diffmap_all <- DiffusionMap(data = all_sig)

eigs_all <- as.data.frame(eigenvectors(diffmap_all)


write.csv(eigs_all,file='data/p_all_sig_dm.csv')


