library(Seurat)
library(scRepertoire)
library(devtools)
library(circlize)
library(scales)
#devtools::install_github("ncborcherding/scRepertoire@dev")
#manual: https://ncborcherding.github.io/vignettes/vignette.html
#not install from bioconductor,it is from github. 
setwd("~/gse162499/scRepertoire_finalize")

P57_B <- read.csv("GSM4952973_P57_Blood_filtered_contig_annotations.csv")
P57_T <- read.csv("GSM4952972_P57_Tumor_filtered_contig_annotations.csv")
P58_B <- read.csv("GSM4952975_P58_Blood_filtered_contig_annotations.csv")
P58_T <- read.csv("GSM4952974_P58_Tumor_filtered_contig_annotations.csv")
P60_B <- read.csv("GSM4952978_P60_Blood_filtered_contig_annotations.csv")
P60_T <- read.csv("GSM4952976_P60_Tumor_filtered_contig_annotations.csv")
P61_B <- read.csv("GSM4952981_P61_Blood_filtered_contig_annotations.csv")
P61_T <- read.csv("GSM4952979_P61_Tumor_filtered_contig_annotations.csv")
contig_list <- list(P57_B,P57_T,P58_B,P58_T,P60_B,P60_T,P61_B,P61_T)
combined <- combineTCR(contig_list, 
                       samples = c("P57", "P57", "P58", "P58", "P60","P60","P61","P61"), 
                       ID = c( "B", "T", "B", "T","B","T","B","T"))
quantContig(combined, cloneCall="gene+nt", scale = T)

vizGenes(combined, gene = "V", 
         chain = "TRA", 
         plot = "bar", 
         order = "variance", 
         scale = TRUE)
hms_cluster_id <- readRDS(file="Only_T_cluster_id_test.rds")
seurat <- readRDS(file="Only_T_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
table(Idents(seurat))
#add seurat
seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", group.by = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
                                            "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
                                            "#ED7953FF","#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

names(seurat@meta.data)
head(seurat@meta.data)

DimPlot(seurat, group.by = "tech") +
  scale_color_manual(values=colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

slot(seurat, "meta.data")$cloneType <- factor(slot(seurat, "meta.data")$cloneType, 
                                              levels = c("Hyperexpanded (100 < X <= 500)", 
                                                         "Large (20 < X <= 100)", 
                                                         "Medium (5 < X <= 20)", 
                                                         "Small (1 < X <= 5)", 
                                                         "Single (0 < X <= 1)", NA))
DimPlot(seurat, group.by = "cloneType") +
  scale_color_manual(values = colorblind_vector(5), na.value="grey") + 
  theme(plot.title = element_blank())

#只能做出blood和腫瘤旁邊相近的clonotype,不能做出全部
compareClonotypes(combined, numbers = 5, samples = c("P57_B","P57_T","P58_B","P58_T",
                                                     "P60_B","P60_T","P61_B","P61_T"),
                  cloneCall="aa", graph = "alluvial")

library(ggraph)
#這裏的aa:CALSGVTSYDKVIF_CASSLLGGGNNEQFF是以前的結果
#saveRDS(seurat, file = "seurat.rds")
seurat <- highlightClonotypes(seurat,  cloneCall= "aa", sequence = c("CALSGVTSYDKVIF_CASSLLGGGNNEQFF","CAASRNAGNMLTF_CASSISGTGEIGEAFF",
                                                                    "CAEGALGSYIPTF;CAFYSSASKIIF_CSGGTSGYEQYF","CALLGGINTGNQFYF_CASISWDRETGHPLHF","CALSEAGAAGNKLTF_CASSPPLGDTEAFF",
                                                                   "CALSEVNQAGTALIF_CASSDAGVSTNEKLFF","CAMREGNTGGFKTIF_CASSQEDRVEETQYF"))

DimPlot(seurat, group.by = "highlight", pt.size = 0.5)
theme(plot.title = element_blank())
occupiedscRepertoire(seurat, x.axis = "ident")

alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("celltype", "ident", "tech"), 
                   color = "TRAV12-2.TRAJ42.TRAC_TRBV20-1.TRBD2.TRBJ2-3.TRBC2") + 
  scale_fill_manual(values = c("grey", colorblind_vector(2)[2]))

alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("celltype", "ident", "tech"), 
                   color = "ident") 
combined2 <- expression2List(seurat, split.by = "ident")
length(combined2)
clonalDiversity(combined2, cloneCall = "nt")
clonalHomeostasis(combined2, cloneCall = "nt")
clonalProportion(combined2, cloneCall = "nt")
clonalOverlap(combined2, cloneCall="aa", method="overlap")
occupiedscRepertoire(seurat, x.axis = "ident")

library(circlize)
library(scales)
circles <- getCirclize(seurat, 
                       group.by = "ident")
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)
circlize::chordDiagram(circles,
                       self.link = 1, 
                       grid.col = grid.cols)

StartracDiversity(seurat, 
                  type = "tech", 
                  sample = "celltype", 
                  by = "overall")

library(ggraph)
clonalNetwork(seurat,       reduction = "umap", 
              identity = "ident", filter.clones = NULL,
              filter.identity = NULL, cloneCall = "aa")

#only cd8+ T cells
Only_cd8<-subset(seurat, idents=c("Resident_Memory_CD8_T",
                                        "Terminally_Exhausted_CD8_T",
                                        "Cytotoxic_CD8_T","Effector_Memory_CD8_T",
                                        "Pre_Exhausted_CD8_T"))
DimPlot(Only_cd8, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Only_cd8, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(Only_cd8, file = "Only_cd8_cluster_id_test.rds")
seurat <- readRDS(file="Only_cd8_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = F)
table(Idents(seurat))
seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", group.by = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
circles <- getCirclize(seurat,  group.by = "ident")
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)
circlize::chordDiagram(circles, self.link = 1, 
                       grid.col = grid.cols)
clonalNetwork(seurat, reduction = "umap", identity = "ident",
              filter.clones = NULL,filter.identity = NULL,
              cloneCall = "aa")
StartracDiversity(seurat, type = "tech", 
                  sample = "celltype", by = "overall")
alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("celltype", "ident", "tech"),   color = "ident") 
combined2 <- expression2List(seurat, split.by = "ident")
length(combined2)
clonalDiversity(combined2, cloneCall = "nt")
clonalHomeostasis(combined2, cloneCall = "nt")
clonalProportion(combined2, cloneCall = "nt")
clonalOverlap(combined2, cloneCall="aa", method="overlap")
occupiedscRepertoire(seurat, x.axis = "ident")
occupiedscRepertoire(seurat, label= FALSE,x.axis = "ident")

#only cd4+ T cells
hms_cluster_id <- readRDS(file="Only_T_cluster_id_test.rds")
Only_cd4<-subset(hms_cluster_id, idents=c("Naive_CD4_T", "T_Helper",
                                        "Effector_Memory_CD4_T",
                                      "Regulatory_T"))
DimPlot(Only_cd4, reduction = "umap", label = TRUE, pt.size = 0.5) 
DimPlot(Only_cd4, reduction = "umap", label = FALSE, pt.size = 0.5) 
saveRDS(Only_cd4, file = "Only_cd4_cluster_id_test.rds")
seurat <- readRDS(file="Only_cd4_cluster_id_test.rds")
DimPlot(seurat, label = T) + NoLegend()
DimPlot(seurat, label = F) 
table(Idents(seurat))
seurat <- combineExpression(combined, seurat, 
                            cloneCall="gene", group.by = "sample", proportion = FALSE, 
                            cloneTypes=c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))
circles <- getCirclize(seurat,  group.by = "ident")
grid.cols <- scales::hue_pal()(length(unique(seurat@active.ident)))
names(grid.cols) <- levels(seurat@active.ident)
circlize::chordDiagram(circles, self.link = 1, 
                       grid.col = grid.cols)
clonalNetwork(seurat, reduction = "umap", identity = "ident",
              filter.clones = NULL,filter.identity = NULL,
              cloneCall = "aa")
StartracDiversity(seurat, type = "tech", 
                  sample = "celltype", by = "overall")
alluvialClonotypes(seurat, cloneCall = "gene", 
                   y.axes = c("celltype", "ident", "tech"),   color = "ident") 
combined2 <- expression2List(seurat, split.by = "ident")
length(combined2)
clonalDiversity(combined2, cloneCall = "nt")
clonalHomeostasis(combined2, cloneCall = "nt")
clonalProportion(combined2, cloneCall = "nt")
clonalOverlap(combined2, cloneCall="aa", method="overlap")
occupiedscRepertoire(seurat, x.axis = "ident")
occupiedscRepertoire(seurat, label= FALSE,x.axis = "ident")


