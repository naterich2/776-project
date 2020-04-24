library(destiny)
library(ggplot2)
library(plotly)
setwd("~/Documents/SPRING 20/BMI 776/Project")
#total_data <- read.csv('data/total_data.csv')
#total_data <- read.csv('data/GSE80655_data_norm.csv')
#total_meta <- read.csv('data/total_meta.csv')
#total_meta <- read.csv('data/GSE80655_meta.csv')
top_100_age <- read.csv('data/top_100_age.csv')
top_100_diag <- read.csv('data/top_100_diag.csv')
top_100_kwii <- read.csv('data/top_100_kwii.csv')
top_100_tci <- read.csv('data/top_100_tci.csv')

diffmap_age <- DiffusionMap(data = top_100_age)
diffmap_diag <- DiffusionMap(data = top_100_diag)
diffmap_kwii <- DiffusionMap(data = top_100_kwii)
diffmap_tci <- DiffusionMap(data = top_100_tci)

eigs_age <- as.data.frame(eigenvectors(diffmap_age))
eigs_diag <- as.data.frame(eigenvectors(diffmap_diag))
eigs_kwii <- as.data.frame(eigenvectors(diffmap_kwii))
eigs_tci <- as.data.frame(eigenvectors(diffmap_tci))

write.csv(eigs_age,file='data/top_100_age_dm.csv')
write.csv(eigs_diag,file='data/top_100_diag_dm.csv')
write.csv(eigs_kwii,file='data/top_100_kwii_dm.csv')
write.csv(eigs_tci,file='data/top_100_tci_dm.csv')

#plot(diffmap,
#      col=factor(total_meta$age),
#      pal=hcl.colors)
#DC1 <- eigs['DC1']
#DC2 <- eigs['DC2']
#DC3 <- eigs['DC3']
#plot(eigenvalues(diffmap))
#eigs['age'] = total_meta$age
#ggplot(eigs,aes(x=DC1,y=DC2,colors=age)) +
#  geom_point()
#plot(eigenvalues(diffmap))
#plot(diffmap,col=factor(total_meta$age),draw_legend=TRUE,pal=palette(rainbow()))

