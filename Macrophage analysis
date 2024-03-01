##analysis Macrophage data

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
aorta <- readRDS("aorta_dataset.rds", refhook = NULL)

carotid <- readRDS("carotid_dataset.rds", refhook = NULL)
mice_counts <- aorta@assays$RNA@counts
human_counts <- carotid@assays$RNA@counts

MP_seurat <- aorta[, which(aorta@active.ident == "MP")]
#refine neighbors and clusters (parameter must be checked and readjusted)
MP_seurat <- FindNeighbors(MP_seurat, dims = 1:20)
MP_seurat <- FindClusters(MP_seurat, resolution = 0.6)
DimPlot(MP_seurat, reduction = "umap",label = TRUE, repel = TRUE)

MP_counts <- mice_counts[, which(colnames(mice_counts) %in% colnames(MP_seurat))]

#create cds objects
# create CDS
#transform mice counts in matrix
expression_matrix <- as.matrix(MP_counts)
#create gene annotation
gene_annotation <- data.frame(rownames(MP_counts))
colnames(gene_annotation) <- c("gene_short_name")
rownames(gene_annotation) <- rownames(expression_matrix)
MP_metadata <- data.frame(MP_seurat$time, MP_seurat$gender)
colnames(MP_metadata) <- c("time", "gender")
cds <- new_cell_data_set(expression_matrix, cell_metadata = MP_metadata, gene_metadata = gene_annotation)

#Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 30)
#check if you are using enough Principal components
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
#try monocle clustering (control resolution, from that it depends clustering numbers)
cds = cluster_cells(cds, resolution=1e-4)

##calculate trajectory using seurat pca/umap data (to verify)
#need to start from the beginning of cds
cds@int_colData@listData$reducedDims@listData$PCA <- MP_seurat@reductions$pca@cell.embeddings
cds@int_colData@listData$reducedDims@listData$UMAP <- MP_seurat@reductions$umap@cell.embeddings
cds = cluster_cells(cds, resolution=1e-4)
MP_clean <- choose_cells(cds)
MP_clean = cluster_cells(MP_clean, resolution=1e-4)
plot_cells(MP_clean, group_label_size = 4)

##bring to seurat object

MP_seurat <- MP_seurat[, which(colnames(MP_seurat) %in% colnames(MP_clean))]
#recluster
MP_seurat <- FindNeighbors(MP_seurat, dims = 1:30)
MP_seurat <- FindClusters(MP_seurat, resolution = 0.4)
DimPlot(MP_seurat, reduction = "umap",label = TRUE, repel = TRUE)
new.cluster.ids <- c("mMP2", "mMP6", "mMP1", "mMP5", "mMP3", "mMP4")
names(new.cluster.ids) <- levels(MP_seurat)
MP_seurat <- RenameIdents(MP_seurat, new.cluster.ids)

MP_seurat@active.ident <- factor( x = MP_seurat @active.ident, levels = c("mMP1","mMP2", "mMP3","mMP4", "mMP5", "mMP6") )

MP <- MP_clean
MP = cluster_cells(MP, resolution=1e-2)

MP@clusters$UMAP$clusters <- MP_seurat@active.ident
MP <- learn_graph(MP)


get_earliest_principal_node <- function(MP, time_bin="20"){
  cell_ids <- which(colData(MP)[, "time"] == time_bin)

  closest_vertex <-
 MP@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(MP), ])
  root_pr_nodes <-
  igraph::V(principal_graph(MP)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
MP <- order_cells(MP, root_pr_nodes=get_earliest_principal_node(MP))
plot_cells(MP, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5,cell_size = 0.7)+ NoLegend() +
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold"))

DimPlot(MP_seurat, label = FALSE, repel = TRUE)+ NoLegend() +
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold"))





##DE
mMP6.markers <- FindMarkers(MP_seurat, ident.1 = "MP6", min.pct = 0.5)
mMP6.markers <- mMP6.markers[mMP6.markers$p_val_adj < 0.05, ]
mMP6.markers <- mMP6.markers[order(mMP6.markers$avg_log2FC, decreasing = TRUE), ]


mMP6.markers$STARNET <- toupper(rownames(mMP6.markers))
up_MP6.markers <- mMP6.markers[mMP6.markers$avg_log2FC > 0.5, ]
down_MP6.markers <- mMP6.markers[mMP6.markers$avg_log2FC  <  -0.5, ]
write.table(mMP6.markers, file="MP6.markers.txt", row.names = T, col.names = T, sep="\t")


mMP5.markers <- FindMarkers(MP_seurat, ident.1 = "MP5", min.pct = 0.5)
mMP5.markers <- mMP5.markers[mMP5.markers$p_val_adj < 0.001, ]
mMP5.markers <- mMP5.markers[order(MP5.markers$avg_log2FC, decreasing = TRUE), ]
mMP5.markers$mean <- rowMeans(MP_seurat[rownames(mMP5.markers) , MP_seurat@active.ident == "MP5"]@assays$RNA@counts)
mMP5.markers$STARNET <- toupper(rownames(MP5.markers))
up_MP5.markers <- MP5.markers[MP5.markers$avg_log2FC > 0.5, ]
down_MP5.markers <- MP5.markers[MP5.markers$avg_log2FC  <  -0.5, ]
write.table(MP5.markers, file="MP5.markers.txt", row.names = T, col.names = T, sep="\t")

##carotid



carotid_MP <- carotid[, which(carotid@active.ident == "MP")]
MP_human_counts <- human_counts[, which(colnames(human_counts) %in% colnames(carotid_MP))]
expression_matrix <- as.matrix(MP_human_counts)
#create gene annotation
gene_annotation <- data.frame(rownames(MP_human_counts))
colnames(gene_annotation) <- c("gene_short_name")
rownames(gene_annotation) <- rownames(expression_matrix)
MP_metadata <- data.frame(carotid_MP$symptom, carotid_MP$gender)
colnames(MP_metadata) <- c("symptom", "gender")

cds <- new_cell_data_set(expression_matrix, cell_metadata = MP_metadata, gene_metadata = gene_annotation)


#Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
#check if you are using enough Principal components
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
#try monocle clustering (control resolution, from that it depends clustering numbers)
cds = cluster_cells(cds, resolution=1e-4)
##calculate trajectory using seurat pca/umap data (to verify)
#need to start from the beginning of cds
cds@int_colData@listData$reducedDims@listData$PCA <- carotid_MP@reductions$pca@cell.embeddings
cds@int_colData@listData$reducedDims@listData$UMAP <- carotid_MP@reductions$umap@cell.embeddings
cds = cluster_cells(cds, resolution=1e-4)
MP_clean <- choose_cells(cds)
MP_clean = cluster_cells(MP_clean, resolution=1e-4)
plot_cells(MP_clean, group_label_size = 4)

##bring to seurat object

carotid_MP <- carotid_MP[, which(colnames(carotid_MP) %in% colnames(MP_clean))]
carotid_MP <- FindNeighbors(carotid_MP, dims = 1:50)
carotid_MP <- FindClusters(carotid_MP, resolution = 0.5)
DimPlot(carotid_MP, reduction = "umap",label = TRUE, repel = TRUE)
new.cluster.ids <- c("hMP6", "hMP3", "hMP7", "hMP5", "hMP4", "hMP2", "hMP8", "hMP1")
names(new.cluster.ids) <- levels(carotid_MP)
carotid_MP <- RenameIdents(carotid_MP, new.cluster.ids)

carotid_MP@active.ident <- factor( x = carotid_MP@active.ident, levels = c("hMP1","hMP2", "hMP3","hMP4", "hMP5", "hMP6", "hMP7", "hMP8") )

cds@clusters$UMAP$clusters <- carotid_MP@active.ident



#trajectory
cds <- learn_graph(cds)
plot_cells(cds, group_label_size = 4, cell_size = 1, label_cell_groups = FALSE)+
NoLegend()+
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold"))



plaque <- table(carotid_MP@active.ident, carotid_MP$symptom)
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


FeaturePlot(carotid_MP, features = "MMP12",cols = c("gray","red")) +
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold"))



##DE
MP6.markers <- FindMarkers(MP_seurat, ident.1 = "MP6", min.pct = 0.5)
MP6.markers <- MP6.markers[MP6.markers$p_val_adj < 0.001, ]
MP6.markers <- MP6.markers[order(MP6.markers$avg_log2FC, decreasing = TRUE), ]
MP6.markers$mean <- rowMeans(MP_seurat[rownames(MP6.markers) , MP_seurat@active.ident == "MP6"]@assays$RNA@counts)

MP6.markers$STARNET <- toupper(rownames(MP6.markers))
up_MP6.markers <- MP6.markers[MP6.markers$avg_log2FC > 0.5, ]
down_MP6.markers <- MP6.markers[MP6.markers$avg_log2FC  <  -0.5, ]
write.table(MP6.markers, file="MP6.markers.txt", row.names = T, col.names = T, sep="\t")


MP5.markers <- FindMarkers(MP_seurat, ident.1 = "MP5", min.pct = 0.5)
MP5.markers <- MP5.markers[MP5.markers$p_val_adj < 0.001, ]
MP5.markers <- MP5.markers[order(MP5.markers$avg_log2FC, decreasing = TRUE), ]
MP5.markers$mean <- rowMeans(MP_seurat[rownames(MP5.markers) , MP_seurat@active.ident == "MP5"]@assays$RNA@counts)
MP5.markers$STARNET <- toupper(rownames(MP5.markers))
up_MP5.markers <- MP5.markers[MP5.markers$avg_log2FC > 0.5, ]
down_MP5.markers <- MP5.markers[MP5.markers$avg_log2FC  <  -0.5, ]
write.table(MP5.markers, file="MP5.markers.txt", row.names = T, col.names = T, sep="\t")

#DE carotid

carotid_MP$exp <- carotid_MP@active.ident 

new.cluster.ids <- c("hMP1","hMP1","hMP3", "hMP4", "hMP5", "hMP6", "hMP7", "hMP8")
names(new.cluster.ids) <- levels(carotid_MP)
carotid_MP <- RenameIdents(carotid_MP, new.cluster.ids)


hMP7.markers <- FindMarkers(carotid_MP, ident.1 = "hMP7", min.pct = 0.5)
hMP7.markers <- hMP7.markers[hMP7..markers$p_val_adj < 0.001, ]
hMP7.markers <- hMP7.markers[order(hMP7.markers$avg_log2FC, decreasing = TRUE), ]
hMP7.markers$mean <- rowMeans(carotid_hMP[rownames(hMP7.markers) , carotid_MP@active.ident == "hMP7"]@assays$RNA@counts)

hMP7.markers$STARNET <- toupper(rownames(hMP7.markers))
up_hMP7.markers <- hMP7.markers[hMP7.markers$avg_log2FC > 0.5, ]
down_hMP7.markers <- hMP7.markers[hMP7.markers$avg_log2FC  <  0.5, ]

write.table(carotid_MP7.markers, file="carotid_MP7.markers.txt", row.names = T, col.names = T, sep="\t")

MP_up_VennDiagram <- venn.diagram(x = list(  a = toupper(rownames(up_mMP5.markers)) , b = toupper(rownames(up_mMP6.markers)) ,c = rownames(up_hMP7s)), filename = NULL, fill = c(progression_list[1:3]), cex = 2, cat.cex = 2, lwd = 2,  lty = 2,fontfamily = "Arial")
grid.newpage()
grid.draw(MP_up_VennDiagram)


##plot_GO
plot_GO <- function(arg1){
mrk_GO<- read.delim(paste(arg1,"_mrk_GO_alt.txt", sep = ""), header = T, sep ="\t")
mrk_GO <- mrk_GO[order(mrk_GO$Adjusted.P.value), ]

df <- mrk_GO[1:5, ]
df<- df[order(-log10(df$P.value), decreasing = TRUE), ]

df$rev <- -log10(df$P.value)

gen_num <- as.data.frame(strsplit(df$Overlap, split = "/"))
gen_num <- t(gen_num)
rownames(gen_num) <- rownames(df)

colnames(gen_num) <-  c('over', "genes_number" )
gen_num <-  apply(gen_num, 2, as.numeric)
gen_num <- as.data.frame(gen_num)
gen_num$gene.ratio <- gen_num$over/gen_num$genes_number
df <- cbind(df, gen_num)
df <-df[order(df$gene.ratio, decreasing = TRUE),]
df$log10_P.value <- df$rev


ggplot(df, aes(gene.ratio, reorder(Term, gene.ratio))) +
  geom_point(aes(col = log10_P.value, size = over)) +
  theme_bw() +
  scale_color_gradient(low = "blue", high = "red")


}

for(i in levels(SMC_seurat)){
  plot_GO(i)
}

for(i in levels(carotid_SMC)){
  plot_GO(i)
}


