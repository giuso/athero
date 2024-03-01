###EC analysis

#load libraries

library(monocle3)
library(velocyto.R)
library(SingleCellExperiment)
library(org.Mm.eg.db)
library(scater)
library(scran)
library(dynamicTreeCut)
library(cluster)
library(Matrix)
library(pagoda2)
library(igraph)
library(umap)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(cowplot)
library(SCENIC)
library(tidyr)
library(vctrs)
library(AUCell)
library(RcisTarget)
library(ggrepel)
library(VennDiagram)
library(network)
library(sna)
library(ggplot2)
library(ggnet)
library(loomR)
library(clustree)
library(RColorBrewer)
library(gplots)



status_list <- c("#007F24", "#BA0303")
time_list <- c("#007F24", "#7EBDC2","#00ADFA","#F2811F", "#BA0303")
gender_list <- c("196343", "AF0200")
progression_list <- c("#FF0000", "#FF7A00","#F8D502","#00621B", "#1DA70B","#00349A","#2773E8","#A566C0", "#ED0FDB","#DB81AC", "#ED0FDB", "EE7E83" )
athero_list <- c("#1e7103", "#f39805", "#de0004", "#760474")
marrow_list <- c("#069247", "#ffe206", "#ff5900", "#c30101","#9502b9" )
venn3_list <- c("#F10101", "#009102",  "#0531F3" )


#load dataset
setwd("~/Documents/single_cell_analysis_documents/seurat_mice")
aorta <- readRDS("aorta_dataset.rds", refhook = NULL)


carotid <- readRDS("carotid_dataset.rds", refhook = NULL)
mice_counts <- aorta@assays$RNA@counts
human_counts <- carotid@assays$RNA@counts

##Mouse
EC_seurat <- aorta[, which(aorta@active.ident == "EC")]
#refine neighbors and clusters (parameter must be checked and readjusted)
EC_seurat <- FindNeighbors(EC_seurat, dims = 1:20)
EC_seurat <- FindClusters(EC_seurat, resolution = 0.4)
DimPlot(EC_seurat, reduction = "umap",label = TRUE, repel = TRUE)
#create EC matrix
EC_counts <- mice_counts[, which(colnames(mice_counts) %in% colnames(EC_seurat))]

#create cds object
expression_matrix <- as.matrix(EC_counts)
gene_annotation <- data.frame(rownames(EC_counts))
colnames(gene_annotation) <- c("gene_short_name")
rownames(gene_annotation) <- rownames(expression_matrix)
EC_metadata <- data.frame(EC_seurat$time, EC_seurat$gender)
colnames(EC_metadata) <- c("time", "gender")
cds <- new_cell_data_set(expression_matrix, cell_metadata = EC_metadata, gene_metadata = gene_annotation)
#Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 30)
#check if you are using enough Principal components
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
cds = cluster_cells(cds, resolution=1e-4)
#need to start from the beginning of cds
cds@int_colData@listData$reducedDims@listData$PCA <- EC_seurat@reductions$pca@cell.embeddings
cds@int_colData@listData$reducedDims@listData$UMAP <- EC_seurat@reductions$umap@cell.embeddings
cds = cluster_cells(cds, resolution=1e-4)
EC_clean <- choose_cells(cds)
EC_clean = cluster_cells(EC_clean, resolution=1e-4)
plot_cells(EC_clean, group_label_size = 4)

##bring to seurat object

EC_seurat <- EC_seurat[, which(colnames(EC_seurat) %in% colnames(EC_clean))]
#recluster
EC_seurat <- FindNeighbors(EC_seurat, dims = 1:30)
EC_seurat <- FindClusters(EC_seurat, resolution = 0.4)
DimPlot(EC_seurat, reduction = "umap",label = TRUE, repel = TRUE)

new.cluster.ids <- c("EC1", "EC3", "EC2", "EC4")
names(new.cluster.ids) <- levels(EC_seurat)
EC_seurat <- RenameIdents(EC_seurat, new.cluster.ids)

EC_seurat@active.ident <- factor( x = EC_seurat@active.ident, levels = new.cluster.ids <- c("EC1",  "EC2","EC3", "EC4"))

EC <- EC_clean
EC = cluster_cells(EC, resolution=1e-4)
EC@clusters$UMAP$clusters <- EC_seurat@active.ident
EC <- learn_graph(EC)

plot_cells(EC, group_label_size = 4, cell_size = 1)+
NoLegend() +
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold")) +
annotate( "text", x = 13, y = -7, size= 7,fontface =2,  label = "EC1") +
annotate("text", x = 12, y = -5, size= 7,fontface =2, label =  "EC2") +
annotate("text", x = 11, y =  -7,size= 7,fontface =2,  label =  "EC3") +
annotate("text", x = 6 , y = -9, size= 7,fontface =2,  label =  "EC4")





####calculate trajectory in carotid EC
carotid_EC <- carotid[, which(carotid@active.ident == "EC")]

carotid_EC_counts <- as.matrix(carotid_EC@assays$RNA@counts)

#refine neighbors and clusters (parameter must be checked and readjusted)
carotid_EC <- FindNeighbors(carotid_EC, dims = 1:20)
carotid_EC <- FindClusters(carotid_EC, resolution = 0.5)
DimPlot(carotid_EC, reduction = "umap",label = TRUE, repel = TRUE)

new.cluster.ids <- c("cEC1", "cEC2", "cEC3", "cEC4", "cEC5", "cEC6")
names(new.cluster.ids) <- levels(carotid_EC)
carotid_EC <- RenameIdents(carotid_EC, new.cluster.ids)

DimPlot(carotid_EC, reduction = "umap",label = FALSE, repel = TRUE)+
NoLegend() +
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold"))



#carotid_EC markers
carotid_EC.markers <- FindAllMarkers(carotid_EC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
carotid_EC.markers <- carotid_EC.markers[carotid_EC.markers$p_val_adj < 0.05, ]
carotid_EC.markers <- carotid_EC.markers[carotid_EC.markers$pct.1 > 0.5, ]
carotid_EC.markers <- carotid_EC.markers[carotid_EC.markers$avg_log2FC > 0.5, ]



all.markers <- FindAllMarkers(carotid_EC, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
all.markers <- all.markers[all.markers$p_val_adj < 0.001, ]

top10 <- carotid_EC.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top5 <- carotid_EC.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(carotid_EC, features = top10$gene) + NoLegend()



###percentage plaque_pct

plaque <- table(carotid_EC@active.ident, carotid_EC$symptom)
plaque <- t(plaque)
plaque <- data.frame(plaque[1,], plaque[2,])
colnames(plaque) <- c("Stable", "Unstable")
stable_pct = plaque$Stable/ rowSums(plaque)*100
unstable_pct = plaque$Unstable/ rowSums(plaque)*100
plaque$stable_pct <- stable_pct
plaque$unstable_pct <- unstable_pct
plaque_pct <- data.frame(plaque$stable_pct, plaque$unstable_pct)
plaque_pct <- t(plaque_pct)
colnames(plaque_pct) <- rownames(plaque)
rownames(plaque_pct) <- c("Stable", "Unstable")
par(mar = c(6,6,8,0) + 1)
barplot(plaque_pct, col = status_list , cex.axis = 2, lwd = 3, cex =2, beside = FALSE)
legend(12,115, legend = rownames(plaque_pct), fill=status_list, text.font = 2, cex= 2, xpd = TRUE)


