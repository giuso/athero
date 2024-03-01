#code to infer single cell data with STARNET GRNs


##module load
#load libraries


library(scater)
library(scran)
library(dynamicTreeCut)
library(cluster)
library(Matrix)
library(igraph)
library(umap)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(Seurat)
library(cowplot)
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
library(RColorBrewer)
library(gplots)
library(fmsb)



###load carotid data

carotid <- readRDS('carotid_dataset.rds')

#case for GRN39
cell_color_cluster_list <- c( "#B3B3B3", progression_list[1], "#0023F6","#1A8340","#B3B3B3","#B3B3B3")
names(cell_color_cluster_list ) <- c("AOR", "SMC", "MP","EC","LIV", "MAM")
  ##load edge
edge39 <- read.csv("edge39.csv", header = T, sep=",")
node39 <- read.csv("node39.csv", header = T, sep=",")
gene39 <- read.csv("gene39.csv", header = T, sep=",")

GRN_names<- as.character(edge39$target)
GRN_names = t(data.frame(strsplit(GRN_names, split='_', fixed=TRUE)))

edge39$tissue <- GRN_names[,1]
edge39$annotation <- GRN_names[,2]
edge39$ENSEMBL <- GRN_names[,3]


GRN_source <- as.character(edge39$source)
GRN_source = t(data.frame(strsplit(GRN_source, split='_', fixed=TRUE)))
edge39$start <- GRN_source[,2]

#edge39 <- edge39[edge39$start %in%edge39$annotation, ]
 

carotid_39 <- carotid[edge39$annotation, carotid@active.ident == "SMC"|carotid@active.ident == "EC"|carotid@active.ident == "Mc"]
carotid_39  <- AverageExpression(carotid_39)

carotid_39 <- as.data.frame(carotid_39$RNA)
colnames(carotid_39) <- c("SMC", "EC", "MP")

carotid_39$Max_score = rowMaxs(as.matrix((carotid_39) ))
carotid_39$Max_score[carotid_39$Max_score == carotid_39$SMC] <- "SMC"
carotid_39$Max_score[carotid_39$Max_score == carotid_39$MP] <- "MP"
carotid_39$Max_score[carotid_39$Max_score == carotid_39$EC] <- "EC"



  ###color different SMC clusters
edge39$tissue[edge39$annotation %in% rownames(carotid_39)[carotid_39$Max_score == "SMC"]] <- "SMC"
edge39$tissue[edge39$annotation %in% rownames(carotid_39)[carotid_39$Max_score == "MP"]] <- "MP"
edge39$tissue[edge39$annotation %in% rownames(carotid_39)[carotid_39$Max_score == "EC"]] <- "EC"
tissue <- unique(edge39$tissue)


key_39 <-read.csv("key_drivers_39.csv", header = T, sep=",")

key_39_filt <-key_39

names_key <- data.frame(annotation_frame[,2],annotation_frame[,2] )
rownames(names_key) <- names_key$annotation_frame...2.
colnames(names_key) <- c("genes", "key")
names_key$key[names_key$key %in% key_39_filt$gene] <- "key"
not_key <- setdiff(names_key$genes, key_39_filt$gene)
names_key$key[names_key$key %in% not_key] <- "not_key"





GRNnode <- as.character(node39$id)
GRNnode = t(data.frame(strsplit(GRNnode, split='_', fixed=TRUE)))
node39$tissue <- GRNnode[,1]
node39$annotation <- GRNnode[,2]
node39$ENSEMBL <- GRNnode[,3]

node39$tissue[node39$annotation %in% rownames(carotid_39)[carotid_39$Max_score == "SMC"]] <- "SMC"
node39$tissue[node39$annotation %in% rownames(carotid_39)[carotid_39$Max_score == "MP"]] <- "MP"
node39$tissue[node39$annotation %in% rownames(carotid_39)[carotid_39$Max_score == "EC"]] <- "EC"

node39$key <- node39$annotation

node39$key <- 5
node39$key[node39$annotation %in% key_39_filt$gene] <- 10



#create GRN network
net <- graph.data.frame(edge39, node39, directed=T)
net <- simplify(net, remove.multiple = F, remove.loops = T)
V(net)$color <- cell_color_cluster_list[V(net)$tissue]

area=vcount(net.bg)^2.4
net.bg <- net
 
l <- layout_with_fr(net.bg)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
par(mfrow=c(1,1), mar=c(0,0,0,0)) # plot two figures - 1 row, 2 columns
plot(net, edge.arrow.size=.4, vertex.label= NA, edge.width = 2, vertex.frame.color = "black",vertex.size = as.numeric(V(net)$key) ,rescale=F, edge.color = "black", layout=l*1)


pdf("grn39_celltype_label.pdf", onefile = TRUE, width=20, height=20)

l <- layout_with_fr(net.bg)
l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
par(mfrow=c(1,1), mar=c(0,0,0,0)) # plot two figures - 1 row, 2 columns
plot(net, edge.arrow.size=.4, vertex.label= node39$annotation , edge.width = 2, vertex.frame.color = "black",vertex.size = as.numeric(V(net)$key), edge.color = "black", layout=l*1)
  
dev.off()

cell_color_cluster_list <- c( "#B39B39B39", progression_list[1], "#00239F6","#1A83940","#B39B39B39","#B39B39B39")
names(cell_color_cluster_list ) <- c("AOR", "SMC", "MP","EC","LIV", "MAM")


##calculate percentage GRN genes expressed in carotid data
data <- data.frame(node39$annotation, node39$tissue)
rownames(data) <- make.unique(node39$annotation)
colnames(data) <- c("gene", "tissue")


df <- data %>% group_by(tissue) %>% # Variable to be transformed
count() %>%ungroup() %>%
mutate(perc = `n` / sum(`n`)) %>%
arrange(perc) %>%
mutate(labels = scales::percent(perc))


ggplot(df, aes(x = "", y = perc, fill = tissue)) +
geom_col(color = "black") +
  coord_polar(theta = "y") + geom_text(aes(label = labels),
        position = position_stack(vjust = 0.5)) + scale_fill_manual(values=cell_color_cluster_list)+
   theme_void()

ggplot(df, aes(x = "", y = perc, fill = tissue)) +
geom_col(color = "black") +
coord_polar(theta = "y")  + scale_fill_manual(values=cell_color_cluster_list)+
theme_void()




##calculate percentage KDR genes expressed in carotid data

key_39$tissue[key_39$gene %in% rownames(carotid_39)[carotid_39$Max_score == "SMC"]] <- "SMC"
key_39$tissue[key_39$gene %in% rownames(carotid_39)[carotid_39$Max_score == "MP"]] <- "MP" 
key_39$tissue[key_39$gene %in% rownames(carotid_39)[carotid_39$Max_score == "EC"]] <- "EC"
   
   
data <- data.frame(key_39$gene, key_39$tissue)
colnames(data) <- c("gene", "tissue")


df <- data %>% group_by(tissue) %>% # Variable to be transformed
    count() %>%
    ungroup() %>%
    mutate(perc = `n` / sum(`n`)) %>%
    arrange(perc) %>%
    mutate(labels = scales::percent(perc))


ggplot(df, aes(x = "", y = perc, fill = tissue)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") + geom_text(aes(label = labels),
  position = position_stack(vjust = 0.5)) + scale_fill_manual(values=cell_color_cluster_list)+
  theme_void()


ggplot(df, aes(x = "", y = perc, fill = tissue)) +
  geom_col(color = "black") +
  coord_polar(theta = "y") + scale_fill_manual(values=cell_color_cluster_list)+
  theme_void()

##radar chart of GRN phenotypic data
data <- read.delim("SYNTAX.txt", header = T, sep="\t", row.names = 1)

grn_39 <- data$GRN.39
names(grn_39) <- rownames(data)
grn_39_data <- data.frame(rbind(rep(100,13), rep(0,13), grn_39))

 radarchart(grn_39_data, axistype=1 ,

    #custom polygon
    pcol=rgb(0.2,0.5,0.5,0.9) , pfcol=rgb(0.2,0.5,0.5,0.5) , plwd=4 ,

    #custom the grid
    cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,100,5), cglwd=0.8,

    #custom labels
    vlcex=0.8
    )






