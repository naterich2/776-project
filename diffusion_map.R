library(destiny)
library(ggplot2)
library(plotly)
setwd("~/Documents/SPRING 20/BMI 776/Project")
total_data <- read.csv('data/total_data.csv')
#total_data <- read.csv('data/GSE80655_data_norm.csv')
total_meta <- read.csv('data/total_meta.csv')
#total_meta <- read.csv('data/GSE80655_meta.csv')
diffmap <- DiffusionMap(data = total_data)
plot(diffmap,
      col=factor(total_meta$age),
      pal=hcl.colors)
eigs <- as.data.frame(eigenvectors(diffmap))
write.csv(eigs,file='data/diffusion_map.csv')
DC1 <- eigs['DC1']
DC2 <- eigs['DC2']
DC3 <- eigs['DC3']
plot(eigenvalues(diffmap))
eigs['age'] = total_meta$age
ggplot(eigs,aes(x=DC1,y=DC2,colors=age)) +
  geom_point()
plot(eigenvalues(diffmap))
plot(diffmap,col=factor(total_meta$age),draw_legend=TRUE,pal=palette(rainbow()))

