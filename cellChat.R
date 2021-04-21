############ This script performs CellChat Analysis on merged Seurat objects, implemented here to analyze the combined samples from a single timepoint (here, P0) ###############
############ Please refer to CellChat [v1.0.0] for details on parameters and functions ############


### Set working directory and load libraries ###
setwd('/path_to_directory_where_Seurat_objects_are_stored')
library(Seurat)
library(CellChat)
library(patchwork)

###### Load and Merge Seurat Objects #######

load("/path/to/Sample1_Seurat_Object.Robj")
Sample1 <- M
Sample1= UpdateSeuratObject(object = Sample1)

load("/path/to/Sample2_Seurat_Object.Robj")
Sample2 <- M
Sample2= UpdateSeuratObject(object = Sample2)

load("/path/to/Sample3_Seurat_Object.Robj")
Sample3 <- M
Sample3= UpdateSeuratObject(object = Sample3)

P0 <- merge(Sample1, y = c(Sample2, Sample3), add.cell.ids = c("Sample1", "Sample2", "Sample3"))
P0 <- NormalizeData(object = P0)
P0 <- FindVariableFeatures(object = P0)
P0 <- ScaleData(object = P0)
P0 <- RunPCA(object = P0)
P0 <- FindNeighbors(object = P0)
P0 <- FindClusters(object = P0)
P0 <- RunTSNE(object = P0)
DimPlot(object = P0, reduction = "tsne")

P0.markers <- FindAllMarkers(P0, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


##Assign cell types to each merged Seurat cluster #####
new.cluster.ids <- c("Cell Type 1", "Cell Type 2", "Cell Type 3")
names(new.cluster.ids) <- levels(P0)
P0 <- RenameIdents(P0, new.cluster.ids)
DimPlot(object = P0, label = TRUE, pt.size = 0.5) + NoLegend()

P0$CellType <- Idents(P0)

P0 <- SetIdent(P0, value = P0@meta.data$CellType)
head(P0@active.ident)

### Perform CellChat analysis [v1.0.0] with default parameters ###
cellchatP0 <- createCellChat(P0, group.by = "CellType")
cellchat <- setIdent(cellchatP0, ident.use = "CellType")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) 
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
pathways.show <- c("WNT") ### Input pathway of interest 
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.WNT <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE) ###Input pathway of interest
LR.show <- pairLR.WNT[1,]
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

pathways.show.all <- cellchat@netP$pathways


cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
