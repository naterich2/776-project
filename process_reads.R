library("DESeq2")
setwd("~/Documents/SPRING 20/BMI 776/Project")
GSE80655 <- read.csv('data/GSE80655_data.csv')
GSE80655_meta <- read.csv('data/GSE80655_meta.csv')

dds1 <- DESeqDataSetFromMatrix(countData=GSE80655,
                              colData=GSE80655_meta,
                              design=~diagnosis)

dseq1 <- DESeq(dds1)

counts_data1 = as.data.frame(counts(dseq1,normalized=TRUE))
write.csv(counts_data1,file='data/GSE80655_data_norm.csv')

GSE42546 <- read.csv('data/GSE42546_data.csv')
GSE42546_meta <- read.csv('data/GSE42546_meta.csv')

dds2 <- DESeqDataSetFromMatrix(countData=GSE42546,
                               colData=GSE42546_meta,
                               design=~diagnosis)

dseq2 <- DESeq(dds2)

counts_data2 = as.data.frame(counts(dseq2,normalized=TRUE))
write.csv(counts_data2,file='data/GSE42546_data_norm.csv')
