
############## Script to Run MNN (Mutual Nearest Neighbor) analysis
############## Seurat >v3.0.0 Objects required


############ Load libraries
library(Seurat)
library(SeuratWrappers)
library(dplyr)

setwd(Directory of Interest)

load('Dataset1 Seurat Object')
Dataset1$Celltypes <- Dataset1$seurat_clusters
Dataset1 <- SCTransform(Dataset1)

load('Dataset2 Seurat Object')
Dataset2$Celltypes <- Dataset2$seurat_clusters
Dataset2 <- SCTransform(Dataset2)
 
MNN <- RunFastMNN(object.list = c(Dataset1, Dataset2), features = 3000) ########### Run MNN
head(MNN@meta.data)
save(MNN, file = 'MNNCorrected.Robj')
write.table(data.frame(MNN@tools$RunFastMNN@metadata$merge.info$pairs[[1]]), "MNNPairs.txt", sep="\t", quote=F)