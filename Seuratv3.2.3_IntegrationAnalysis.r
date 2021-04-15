
#################### Seurat v3.2.3 Integration Analysis ###############
#################### Read Counts Matrix and metadata to run Canonical correlation based Integration #######

############## Load libraries ##############
library(Biobase)
library(Seurat)
library(cowplot)
library(plotly)


########################## Read Counts Data ##################
E14KI.data <- as.matrix(read.table(file = "10X_E14Six2KI_3mm_RawCounts_Final_Matrix.txt", sep = "\t", header=T, row.names=1, check.names=F))   ###### Raw counts
E14TGC.data <- as.matrix(read.table(file = "10X_E14Six2TGC_3mm_RawCounts_Final_Matrix.txt", sep = "\t", header=T, row.names=1, check.names=F)) ###### Raw counts
E14TSC.data <- as.matrix(read.table(file = "10X_E14Six2TGC_TSCF_3mm_RawCounts_Final_Matrix.txt", sep = "\t", header=T, row.names=1, check.names=F))  ###### Raw counts
P0KI.data <- as.matrix(read.table(file = "10X_P0Six2KI_3mm_RawCounts_Final_Matrix.txt", sep = "\t", header=T, row.names=1, check.names=F))   ###### Raw counts
P0TGC.data <- as.matrix(read.table(file = "10X_P0Six2TGC_3mm_RawCounts_Final_Matrix.txt", sep = "\t", header=T, row.names=1, check.names=F))  ###### Raw counts
P0TSC.data <- as.matrix(read.table(file = "10X_P0Six2TGC_TSC1_3mm_RawCounts.txt", sep = "\t", header=T, row.names=1, check.names=F))   ###### Raw counts
##############################################################


########################### Read Meta File ####################
E14.KImeta <- as.matrix(read.table(file = "MetaFile_E14_KI.txt", sep = "\t", header=T, row.names=1, check.names=F))                          #### Sample Info File
col_d <- data.frame(E14.KImeta)
metaE14KI <- new("AnnotatedDataFrame", data = col_d)
rownames(metaE14KI) <- metaE14KI$Cells

E14.TGCmeta <- as.matrix(read.table(file = "MetaFile_E14_TGC.txt", sep = "\t", header=T, row.names=1, check.names=F))                         #### Sample Info File
col_d <- data.frame(E14.TGCmeta)
metaE14TGC <- new("AnnotatedDataFrame", data = col_d)
rownames(metaE14TGC) <- metaE14TGC$Cells

E14.TSCmeta <- as.matrix(read.table(file = "MetaFile_E14_TSC.txt", sep = "\t", header=T, row.names=1, check.names=F))                         #### Sample Info File
col_d <- data.frame(E14.TSCmeta)
metaE14TSC <- new("AnnotatedDataFrame", data = col_d)
rownames(metaE14TSC) <- metaE14TSC$Cells

P0.KImeta <- as.matrix(read.table(file = "MetaFile_P0KI.txt", sep = "\t", header=T, row.names=1, check.names=F))                          #### Sample Info File
col_d <- data.frame(P0.KImeta)
metaP0KI <- new("AnnotatedDataFrame", data = col_d)
rownames(metaP0KI) <- metaP0KI$Cells

P0.TGCmeta <- as.matrix(read.table(file = "MetaFile_P0TGC.txt", sep = "\t", header=T, row.names=1, check.names=F))                          #### Sample Info File
col_d <- data.frame(P0.TGCmeta)
metaP0TGC <- new("AnnotatedDataFrame", data = col_d)
rownames(metaP0TGC) <- metaP0TGC$Cells

P0.TSCmeta <- as.matrix(read.table(file = "MetaFile_P0TSC.txt", sep = "\t", header=T, row.names=1, check.names=F))                          #### Sample Info File
col_d <- data.frame(P0.TSCmeta)
metaP0TSC <- new("AnnotatedDataFrame", data = col_d)
rownames(metaP0TSC) <- metaP0TSC$Cells
##############################################################

setwd('Directory of Interest')
######################### Set Seurat Objects ##################
E14KI <- CreateSeuratObject(counts = E14KI.data, project = "E14_KI", min.cells = 3, min.features=100)
CellsMeta = E14KI@meta.data
CellsMeta["Cluster"] <- metaE14KI$Cluster
CellsMetaTrim <- subset(CellsMeta, select = c("Cluster"))
E14KI <- AddMetaData(E14KI, CellsMetaTrim)
CellsMeta["Type"] <- metaE14KI$Type
CellsMetaTrim <- subset(CellsMeta, select = c("Type"))
E14KI <- AddMetaData(E14KI, CellsMetaTrim)

E14TGC <- CreateSeuratObject(counts = E14TGC.data, project = "E14_TGC", min.cells = 3, min.features=100)
CellsMeta = E14TGC@meta.data
CellsMeta["Cluster"] <- metaE14TGC$Cluster
CellsMetaTrim <- subset(CellsMeta, select = c("Cluster"))
E14TGC <- AddMetaData(E14TGC, CellsMetaTrim)
CellsMeta["Type"] <- metaE14TGC$Type
CellsMetaTrim <- subset(CellsMeta, select = c("Type"))
E14TGC <- AddMetaData(E14TGC, CellsMetaTrim)

E14TSC <- CreateSeuratObject(counts = E14TSC.data, project = "E14_TSC", min.cells = 3, min.features=100)
CellsMeta = E14TSC@meta.data
CellsMeta["Cluster"] <- metaE14TSC$Cluster
CellsMetaTrim <- subset(CellsMeta, select = c("Cluster"))
E14TSC <- AddMetaData(E14TSC, CellsMetaTrim)
CellsMeta["Type"] <- metaE14TSC$Type
CellsMetaTrim <- subset(CellsMeta, select = c("Type"))
E14TSC <- AddMetaData(E14TSC, CellsMetaTrim)

P0KI <- CreateSeuratObject(counts = P0KI.data, project = "P0_KI", min.cells = 3, min.features=100)
CellsMeta = P0KI@meta.data
CellsMeta["Cluster"] <- metaP0KI$Cluster
CellsMetaTrim <- subset(CellsMeta, select = c("Cluster"))
P0KI <- AddMetaData(P0KI, CellsMetaTrim)
CellsMeta["Type"] <- metaP0KI$Type
CellsMetaTrim <- subset(CellsMeta, select = c("Type"))
P0KI <- AddMetaData(P0KI, CellsMetaTrim)

P0TGC <- CreateSeuratObject(counts = P0TGC.data, project = "P0_TGC", min.cells = 3, min.features=100)
CellsMeta = P0TGC@meta.data
CellsMeta["Cluster"] <- metaP0TGC$Cluster
CellsMetaTrim <- subset(CellsMeta, select = c("Cluster"))
P0TGC <- AddMetaData(P0TGC, CellsMetaTrim)
CellsMeta["Type"] <- metaP0TGC$Type
CellsMetaTrim <- subset(CellsMeta, select = c("Type"))
P0TGC <- AddMetaData(P0TGC, CellsMetaTrim)

P0TSC <- CreateSeuratObject(counts = P0TSC.data, project = "P0_TSC", min.cells = 3, min.features=100)
CellsMeta = P0TSC@meta.data
CellsMeta["Cluster"] <- metaP0TSC$Cluster
CellsMetaTrim <- subset(CellsMeta, select = c("Cluster"))
P0TSC <- AddMetaData(P0TSC, CellsMetaTrim)
CellsMeta["Type"] <- metaP0TSC$Type
CellsMetaTrim <- subset(CellsMeta, select = c("Type"))
P0TSC <- AddMetaData(P0TSC, CellsMetaTrim)
##################################################################


################# Subset Cells, Normalize and Find Variable Genes ########
E14KI[["percent.mt"]] <- PercentageFeatureSet(E14KI, pattern = "^mt-")
VlnPlot(E14KI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('E14KI_QC_Stats_Seuratv3.png')
E14KI <- subset(E14KI, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 5.5)
E14KI <- NormalizeData(E14KI, verbose = FALSE)
E14KI <- FindVariableFeatures(E14KI, selection.method = "vst", nfeatures = 2000)

E14TGC[["percent.mt"]] <- PercentageFeatureSet(E14TGC, pattern = "^mt-")
VlnPlot(E14TGC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('E14TGC_QC_Stats_Seuratv3.png')
E14TGC <- subset(E14TGC, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 6)
E14TGC <- NormalizeData(E14TGC, verbose = FALSE)
E14TGC <- FindVariableFeatures(E14TGC, selection.method = "vst", nfeatures = 2000)

E14TSC[["percent.mt"]] <- PercentageFeatureSet(E14TSC, pattern = "^mt-")
VlnPlot(E14TSC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('E14TSC_QC_Stats_Seuratv3.png')
E14TSC <- subset(E14TSC, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 5)
E14TSC <- NormalizeData(E14TSC, verbose = FALSE)
E14TSC <- FindVariableFeatures(E14TSC, selection.method = "vst", nfeatures = 2000)

P0KI[["percent.mt"]] <- PercentageFeatureSet(P0KI, pattern = "^mt-")
VlnPlot(P0KI, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('P0KI_QC_Stats_Seuratv3.png')
P0KI <- subset(P0KI, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 4)
P0KI <- NormalizeData(P0KI, verbose = FALSE)
P0KI <- FindVariableFeatures(P0KI, selection.method = "vst", nfeatures = 2000)

P0TGC[["percent.mt"]] <- PercentageFeatureSet(P0TGC, pattern = "^mt-")
VlnPlot(P0TGC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('P0TGC_QC_Stats_Seuratv3.png')
P0TGC <- subset(P0TGC, subset = nFeature_RNA > 100 & nFeature_RNA < 5500 & percent.mt < 5.5)
P0TGC <- NormalizeData(P0TGC, verbose = FALSE)
P0TGC <- FindVariableFeatures(P0TGC, selection.method = "vst", nfeatures = 2000)

P0TSC[["percent.mt"]] <- PercentageFeatureSet(P0TSC, pattern = "^mt-")
VlnPlot(P0TSC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave('P0TSC_QC_Stats_Seuratv3.png')
P0TSC <- subset(P0TSC, subset = nFeature_RNA > 100 & nFeature_RNA < 6500 & percent.mt < 6)
P0TSC <- NormalizeData(P0TSC, verbose = FALSE)
P0TSC <- FindVariableFeatures(P0TSC, selection.method = "vst", nfeatures = 2000)
###################################################################################


##################### Perform Integration ##########################################
Integration.anchors <- FindIntegrationAnchors(object.list = list(E14KI, E14TGC, E14TSC, P0KI, P0TGC, P0TSC), dims = 1:30)
Integration.combined <- IntegrateData(anchorset = Integration.anchors, dims = 1:30)
####################################################################################


#################### Set Assay, Regress CC, Scale, UMAP and Clustering ##############
DefaultAssay(Integration.combined) <- "integrated"
#######Cell Cycle Removal #######
cc.genes <- readLines(con = "CellCycle_Genes_Mouse.txt")
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]
Integration.combined <- CellCycleScoring(object = Integration.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Integration.combined@meta.data$CC.Difference <- Integration.combined@meta.data$S.Score - Integration.combined@meta.data$G2M.Score
##################################
Integration.combined_aftercc <- ScaleData(Integration.combined, vars.to.regress=c('CC.Difference'), verbose = TRUE)
Integration.combined_aftercc <- RunPCA(Integration.combined_aftercc, npcs = 30, verbose = TRUE)
Integration.combined_aftercc <- RunUMAP(Integration.combined_aftercc, reduction = "pca", dims = 1:30)
Integration.combined_aftercc <- FindNeighbors(Integration.combined_aftercc, reduction = "pca", dims = 1:30)
Integration.combined_aftercc <- FindClusters(Integration.combined_aftercc, resolution = 0.5)
###################################################################################


##################### Visualization and Marker Prediction ##########################
p1 <- DimPlot(Integration.combined_aftercc, reduction = "umap", group.by = "Type")
p2 <- DimPlot(Integration.combined_aftercc, reduction = "umap", label = TRUE)
plot_grid(p1, p2)
ggsave('IntegrationPlot_All6Samples_afterCCremoval.png', width=15, height=15)

### visualize the two conditions side-by-side
DimPlot(Integration.combined_aftercc, reduction = "umap", split.by = "Type")
ggsave('IntegrationPlot_All6Samples_afterCCremoval_Comparison.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'IntegrationPlot_All8Samples_afterCCremoval_Comparison.html')

DimPlot(Integration.combined_aftercc, reduction = "umap", group.by="Cluster")
ggsave('IntegrationPlot_All6Samples_afterCCremoval_WithOriginalClusters.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'IntegrationPlot_All8Samples_afterCCremoval_WithOriginalClusters.html')

DimPlot(Integration.combined_aftercc, reduction = "umap", group.by="Type")
ggsave('IntegrationPlot_All6Samples_afterCCremoval_WithType.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'IntegrationPlot_All8Samples_afterCCremoval_WithType.html')

DimPlot(Integration.combined_aftercc, reduction = "umap", group.by="Phase")
ggsave('IntegrationPlot_All6Samples_afterCCremoval_CellCyclePhase.png', width=15, height=15)
htmlwidgets::saveWidget(as.widget(ggplotly()), 'IntegrationPlot_All8Samples_afterCCremoval_CellCyclePhase.html')

Idents(Integration.combined_aftercc) <- Integration.combined_aftercc@meta.data$seurat_clusters

save(Integration.combined_aftercc,file = 'IntegrationPlot_All6Samples_afterCCremoval.Robj')
#################################################