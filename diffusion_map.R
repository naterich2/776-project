library(destiny)
library(ggplot2)
setwd("~/Documents/SPRING 20/BMI 776/Project")
total_data <- read.csv('data/total_data.csv')
#total_data <- read.csv('data/GSE80655_data_norm.csv')
total_meta <- read.csv('data/total_meta.csv')
#total_meta <- read.csv('data/GSE80655_meta.csv')
diffmap <- DiffusionMap(data = total_data)
plot(diffmap,c
     ol=factor(total_meta$age),
     pal=hcl.colors)
DPT(diffmap)
ggplot(diffmap)
