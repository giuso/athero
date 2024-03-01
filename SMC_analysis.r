#### script to analyze SMC set of carotid and aorta

##load libraries
library(monocle3)
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



SMC_seurat <- aorta[, which(aorta@active.ident == "SMC")]
#refine neighbors and clusters (parameter must be checked and readjusted)
SMC_seurat <- FindNeighbors(SMC_seurat, dims = 1:20)
SMC_seurat <- FindClusters(SMC_seurat, resolution = 0.8)
DimPlot(SMC_seurat, reduction = "umap",label = FALSE, repel = TRUE)
SMC_counts <- mice_counts[, which(colnames(mice_counts) %in% colnames(SMC_seurat))]


SMC_seurat<- SMC_seurat[rowSums(SMC_seurat@assays$RNA@counts)>300,]
SMC_seurat= SMC_seurat[-grep("ERCC", rownames(SMC_seurat) ), ]


#create cds object
# create CDS
#transform mice counts in matrix
expression_matrix <- as.matrix(SMC_counts)
#create gene annotation
gene_annotation <- data.frame(rownames(SMC_counts))
colnames(gene_annotation) <- c("gene_short_name")
rownames(gene_annotation) <- rownames(expression_matrix)
SMC_metadata <- data.frame(SMC_seurat$time, SMC_seurat$gender)
colnames(SMC_metadata) <- c("time", "gender")

cds <- new_cell_data_set(expression_matrix, cell_metadata = SMC_metadata, gene_metadata = gene_annotation)


#Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
#check if you are using enough Principal components
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
#try monocle clustering (control resolution, from that it depends clustering numbers)
cds = cluster_cells(cds, resolution=1e-4)
##calculate trajectory using seurat pca/umap data (to verify)
#need to start from the beginning of cds
cds@int_colData@listData$reducedDims@listData$PCA <- SMC_seurat@reductions$pca@cell.embeddings
cds@int_colData@listData$reducedDims@listData$UMAP <- SMC_seurat@reductions$umap@cell.embeddings
cds = cluster_cells(cds, resolution=1e-4)
##clean outlier cells
SMC_clean <- choose_cells(cds)
SMC_clean = cluster_cells(SMC_clean, resolution=1e-4)
plot_cells(SMC_clean, group_label_size = 4)+
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold"))


##bring to seurat object
SMC_seurat <- SMC_seurat[, which(colnames(SMC_seurat) %in% colnames(SMC_clean))]
#recluster
SMC_seurat <- FindNeighbors(SMC_seurat, dims = 1:50)
SMC_seurat <- FindClusters(SMC_seurat, resolution = 0.6)
DimPlot(SMC_seurat, reduction = "umap",label = TRUE, repel = TRUE)
new.cluster.ids <- c( "mSMC1", "mSMC2", "mSMC6", "mSMC3", "mSMC4", "mSMC5", "mSMC7")
names(new.cluster.ids) <- levels(SMC_seurat)
SMC_seurat <- RenameIdents(SMC_seurat, new.cluster.ids)
SMC_seurat@active.ident <- factor( x = SMC_seurat@active.ident, levels = c( "mSMC1","mSMC2" , "mSMC3", "mSMC4", "mSMC5","SmMC6", "mSMC7") )
#plot
DimPlot(SMC_seurat, reduction = "umap",label = FALSE, repel = TRUE) + NoLegend() +
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold"))



get_earliest_principal_node <- function(SMC_clean, time_bin="20"){
  cell_ids <- which(colData(SMC_clean)[, "time"] == time_bin)

  closest_vertex <-
 SMC_clean@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(MP), ])
  root_pr_nodes <-
  igraph::V(principal_graph(MP)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
SMC_clean <- order_cells(SMC_clean, root_pr_nodes=get_earliest_principal_node(SMC_clean))
plot_cells(SMC_clean, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5,cell_size = 0.7)+ NoLegend() +
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold"))




##carotid



carotid_SMC <- carotid[, carotid@active.ident == "SMC"]
SMC_human_counts <- human_counts[, which(colnames(human_counts) %in% colnames(carotid_SMC))]
expression_matrix <- as.matrix(SMC_human_counts)
#create gene annotation
gene_annotation <- data.frame(rownames(SMC_human_counts))
colnames(gene_annotation) <- c("gene_short_name")
rownames(gene_annotation) <- rownames(expression_matrix)
SMC_metadata <- data.frame(carotid_SMC$symptom, carotid_SMC$gender)
colnames(SMC_metadata) <- c("symptom", "gender")

cds <- new_cell_data_set(expression_matrix, cell_metadata = SMC_metadata, gene_metadata = gene_annotation)


#Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50)
#check if you are using enough Principal components
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
#try monocle clustering (control resolution, from that it depends clustering numbers)
cds = cluster_cells(cds, resolution=1e-4)
##calculate trajectory using seurat pca/umap data (to verify)
#need to start from the beginning of cds
cds@int_colData@listData$reducedDims@listData$PCA <- carotid_SMC@reductions$pca@cell.embeddings
cds@int_colData@listData$reducedDims@listData$UMAP <- carotid_SMC@reductions$umap@cell.embeddings
cds = cluster_cells(cds, resolution=1e-4)
##clean outlier cells

carotid_SMC_clean <- choose_cells(cds)
carotid_SMC_clean = cluster_cells(SMC_clean, resolution=1e-4)
plot_cells(SMC_clean, group_label_size = 4)



get_earliest_principal_node <- function(carotid_SMC_clean, time_bin="20"){
  cell_ids <- which(colData(carotid_SMC_clean)[, "time"] == time_bin)

  closest_vertex <-
 carotid_SMC_clean@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(MP), ])
  root_pr_nodes <-
  igraph::V(principal_graph(MP)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}
carotid_SMC_clean <- order_cells(SMC_clean, root_pr_nodes=get_earliest_principal_node(SMC_clean))
plot_cells(SMC_clean, color_cells_by = "pseudotime", label_cell_groups=FALSE, label_leaves=FALSE, label_branch_points=FALSE, graph_label_size=1.5,cell_size = 0.7)+ NoLegend() +
theme(axis.text=element_text(size=20,face="bold"), axis.title=element_text(size=20,face="bold")) +
theme(legend.text=element_text(size=20,face="bold"))


##bring to seurat object
carotid_SMC <- carotid_SMC[, which(colnames(carotid_SMC) %in% colnames(SMC_clean))]
#recluster
carotid_SMC <- FindNeighbors(carotid_SMC, dims = 1:50)
carotid_SMC <- FindClusters(carotid_SMC, resolution = 0.5)
DimPlot(carotid_SMC, reduction = "umap",label = TRUE, repel = TRUE)

new.cluster.ids <- c( "hSMC7","hSMC6","hSMC4", "hSMC5","hSMC2","hSMC1","hSMC8","hSMC3")
names(new.cluster.ids) <- levels(carotid_SMC)
carotid_SMC <- RenameIdents(carotid_SMC, new.cluster.ids)
carotid_SMC@active.ident <- factor( x = carotid_SMC@active.ident, levels = c( "hSMC1","hSMC2","hSMC3", "hSMC4", "hSMC5","hSMC6","hSMC7", "hSMC8"))

plaque <- table(carotid_SMC@active.ident, carotid_SMC$symptom)
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






###markers
##mouse

markers_filt <- function(arg1,arg2){
  sub.markers <- FindMarkers(arg1, ident.1 = arg2, min.pct = 0.5)
sub.markers <- SMC7.markers[sub.markers.markers$p_val_adj < 0.05, ]
sub.markers <- SMC7.markers[order(sub.markers.markers$avg_log2FC, decreasing = TRUE), ]
sub.markers$STARNET <- toupper(rownames(sub.markers.markers))
#up_sub.markers <- sub.markers[sub.markers$avg_log2FC > 0.3, ]
#down_sub.markers.markers <- sub.markers[sub.markers$avg_log2FC  <  -0.3, ]
write.table(sub.markers, file=paste(arg2, "_markers.txt", sep=""), row.names = T, col.names = T, sep="\t")
return(sub.markers)
}

##find and write markers for datasets

for (i in levels(SMC_seurat)){
  markers_filt(SMC_seurat, i)
}

for (i in levels(carotid_SMC)){
  markers_filt(SMC_seurat, i)
}





mSMC6.markers <- markers_filt(SMC_seurat, "mSMC6")
mSMC6.markers <- mSMC6.markers[mSMC6.markers$avg_log2FC > 0.5,]
mSMC7.markers <- markers_filt(SMC_seurat, "mSMC7")
mSMC7.markers <- mSMC7.markers[mSMC7.markers$avg_log2FC > 0.5,]

hSMC7.markers <- markers_filt(SMC_seurat, "hSMC7")
hSMC7.markers <- hSMC7.markers[hSMC7.markers$avg_log2FC > 0.5,]


SMC_up_VennDiagram <- venn.diagram(x = list(  a = rownames(mSMC6.markers),b = rownames(mSMC7.markers),c = rownames(hSMC7.markers) ), filename = NULL, fill = c(progression_list[1:3]), cex = 2, cat.cex = 2, lwd = 2,  lty = 2,fontfamily = "Arial")
grid.newpage()
grid.draw(SMC_up_VennDiagram)



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

for(i in levels(Mc_seurat)){
  plot_GO(i)
}
common_genes <- read.delim("geneTrans.txt", header = T, sep=",")

up_mSMC6.markers_com <- up_mSMC6.markers[rownames(up_mSMC6.markers) %in% common_genes$Gene.name, ]
up_mSMC7.markers_com <- up_mSMC7.markers[rownames(up_mSMC7.markers) %in% common_genes$Gene.name, ]
up_hSMC7.markers_com <- up_hSMC7.markers[rownames(up_hSMC7.markers) %in% common_genes$Gene.name.1, ]

SMC_up_VennDiagram <- venn.diagram(x = list(  a = toupper(rownames(up_mSMC6.markers)) ,b = toupper(rownames(up_mSMC7.markers)),c = rownames(up_hSMC7.markers) ), filename = NULL, fill = c(progression_list[1:3]), cex = 2, cat.cex = 2, lwd = 2,  lty = 2,fontfamily = "Arial")
grid.newpage()
grid.draw(SMC_up_VennDiagram)







