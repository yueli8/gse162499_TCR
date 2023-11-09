setwd("~/gse162498/data_finally")
library(Seurat)
library(monocle3)
library(tidyverse)
library(patchwork)
library(Cairo)
#https://cole-trapnell-lab.github.io/monocle3/docs/installation/
#library(monocle)#can not use this package.error will come out.

#cd4 and cd8 cells Only_T_cluster_id_test.rds
seurat <- readRDS(file="Only_T_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)#pca analysis
#umap,tSNE降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('cds.umap')
p1
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p2
cds <- reduce_dimension(cds, reduction_method = "tSNE")
p3 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="tech")
p3
p4 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype") 
p4
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="Naive_CD4_T", "1"="nkt","2"="nkt","3"="T_Helper","4"= "Effector_Memory_CD4_T",
                                                 "5"="CD8_resident_memory",  "6"="CD4_effector_memory",
                                                 "7"="cd8_terminally_exhausted", "8"="treg","9"="CD4_effector_memory",
                                                 "10"="treg",
                                                 "11"= "cd8_terminally_exhausted", "14"="treg","15"="cd8_cytotoxic",
                                                 "16"="CD8_effector_memory",
                                                 "19"="cd8_terminally_exhausted","20"="cd8_pre_exhausted"  )

colnames(colData(cds))
cds <- learn_graph(cds)
p1= plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p1
p2 = plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster = FALSE,label_cell_groups = FALSE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p2
p3 = plot_cells(cds,color_cells_by = "tech",label_groups_by_cluster = FALSE,label_cell_groups = FALSE, 
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p3
p5 = plot_cells(cds,color_cells_by = "assigned_cell_type",label_cell_groups = FALSE, label_groups_by_cluster = TRUE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 5)
p5
#细胞按拟时排序
cds <- order_cells(cds) #dragging a rectangle,then choose,done
plot_cells(cds, color_cells_by = "pseudotime",  label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE,graph_label_size = 4)

Track<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")

plot_genes_in_pseudotime(cds[Track,], color_cells_by="tech",  min_expr=0.5, ncol = 3)
#FeaturePlot图
plot_cells(cds, genes=Track, show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)

#cd4 cells
seurat <- readRDS(file="cd4_cluster_id.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)#pca analysis
#umap,tSNE降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('cds.umap')
p1
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p2
cds <- reduce_dimension(cds, reduction_method = "tSNE")
p3 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="tech")
p3
p4 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype") 
p4
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "0"="Naive_CD4_T", "3"="T_Helper",
                                              "4"= "CD4_effector_memory","6"="CD4_effector_memory","8"="treg",
                                                 "9"="CD4_effector_memory","10"="treg","14"="treg" )

colnames(colData(cds))
cds <- learn_graph(cds)
p1= plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p1
p2 = plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster = FALSE,label_cell_groups = FALSE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p2
p3 = plot_cells(cds,color_cells_by = "tech",label_groups_by_cluster = FALSE,label_cell_groups = FALSE, 
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p3
p5 = plot_cells(cds,color_cells_by = "assigned_cell_type",label_cell_groups = FALSE, label_groups_by_cluster = TRUE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 5)
p5
#细胞按拟时排序
cds <- order_cells(cds) #dragging a rectangle,then choose,done
plot_cells(cds, color_cells_by = "pseudotime",  label_cell_groups = FALSE,
     label_leaves = FALSE, label_branch_points = FALSE,graph_label_size = 4)

Track<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")

plot_genes_in_pseudotime(cds[Track,], color_cells_by="assigned_cell_type",  min_expr=0.5, ncol = 3)
#FeaturePlot图
plot_cells(cds, genes=Track, show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)


#cd8 cells
seurat <- readRDS(file="cd8_cluster_id.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = T) 
table(Idents(seurat))
##创建CDS对象并预处理数据
data <- GetAssayData(seurat, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
#preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)
plot_pc_variance_explained(cds)#pca analysis
#umap,tSNE降维
cds <- reduce_dimension(cds, preprocess_method = "PCA")
plot_cells(cds)
colnames(colData(cds))
p1 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('cds.umap')
p1
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('cds.umap')
p2
cds <- reduce_dimension(cds, reduction_method = "tSNE")
p3 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="tech")
p3
p4 <- plot_cells(cds, reduction_method="tSNE", color_cells_by="celltype") 
p4
##从seurat导入整合过的umap坐标
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(seurat, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
p2 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="celltype") + ggtitle('int.umap')
p1|p2
p3 <- plot_cells(cds, reduction_method="UMAP", color_cells_by="tech") + ggtitle('int.umap')
p3
## Monocle3聚类分区
cds <- cluster_cells(cds,cluster_method='louvain')#bug only work with cluster_method='louvain'
p1 <- plot_cells(cds, show_trajectory_graph = FALSE) + ggtitle("label by clusterID")
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE) + 
  ggtitle("label by partitionID")
p = wrap_plots(p1, p2)
p
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$seurat_clusters,
                                                 "5"="CD8_resident_memory", "7"="cd8_terminally_exhausted",
                                                 "11"= "cd8_terminally_exhausted","15"="cd8_cytotoxic",
                                                 "16"="CD8_effector_memory",
                                                 "19"="cd8_terminally_exhausted","20"="cd8_pre_exhausted")

colnames(colData(cds))
cds <- learn_graph(cds)
p1= plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, 
               label_branch_points = FALSE)
p1
p2 = plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster = FALSE,label_cell_groups = FALSE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p2
p3 = plot_cells(cds,color_cells_by = "tech",label_groups_by_cluster = FALSE,label_cell_groups = FALSE, 
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 8)
p3
p5 = plot_cells(cds,color_cells_by = "assigned_cell_type",label_cell_groups = FALSE, label_groups_by_cluster = TRUE,
                label_leaves = TRUE, label_branch_points = TRUE,graph_label_size = 5)
p5
#细胞按拟时排序
cds <- order_cells(cds) #dragging a rectangle,then choose,done
plot_cells(cds, color_cells_by = "pseudotime",  label_cell_groups = FALSE,
           label_leaves = FALSE, label_branch_points = FALSE,graph_label_size = 4)

Track<-c("CD8A","CD8B","CD38","CD69","ENTPD1","GZMA","GZMH","MYO1F","SYNE1","TSC22D3","XCL2")

plot_genes_in_pseudotime(cds[Track,], color_cells_by="assigned_cell_type",  min_expr=0.5, ncol = 3)
#FeaturePlot图
plot_cells(cds, genes=Track, show_trajectory_graph=FALSE,label_cell_groups=FALSE,  label_leaves=FALSE)


