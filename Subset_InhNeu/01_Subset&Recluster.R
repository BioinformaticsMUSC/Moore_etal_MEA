suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  #library(speckle)
  #library(magrittr)
  #library(broom)
  #library(muscat)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  #library(scds)
})

setwd("/scratch/subramanian/Haley_MEA")

load("Haley_MEA_Signac_P6/RNA_ATAC_08_Final_Labelled.RData")

table(Haley_MEA_Int$Cell_Class)
#Astro Exc_Neurons Inh_Neurons       Micro       Oligo         OPC 
#7776       10653        4252        4949       17296        3884 
#Vascular 
#1257

#source("/Users/SuganyaSubramanian/Utils.R")
source("/home/subramanian/R/x86_64-pc-linux-gnu-library/4.2/Utils.R")

head(Haley_MEA_Int)

MEA_InhNeu <- subset(Haley_MEA_Int, 
                     subset = Cell_Class == "Inh_Neurons")
table(MEA_InhNeu$Cell_Class)

MEA_InhNeu <- processing_seurat_sctransform(MEA_InhNeu, 
                                            vars_to_regress = c("nUMI","pMito","pRibo",
                                                                "Age","Sex","Hemis","EpDur","SeqBatch"),
                                            npcs = 50, 
                                            res = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))

DefaultAssay(MEA_InhNeu) <- "RNA"
MEA_InhNeu <- NormalizeData(object = MEA_InhNeu, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)

dir.create("output_InhNeu")

pdf("output_InhNeu/01_UMAP_Labelled.pdf", width = 10, height = 6)
p1 <- DimPlot(object = MEA_InhNeu, reduction = "umap",group.by = "Cell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

library(harmony)
MEA_InhNeu_withHar <- MEA_InhNeu

MEA_InhNeu_withHar <- RunHarmony(MEA_InhNeu_withHar, assay.use="RNA", group.by.vars = "Genotype")
Reductions(MEA_InhNeu_withHar)

MEA_InhNeu_withHar <- RunUMAP(MEA_InhNeu_withHar, reduction = "harmony", dims = 1:30)
MEA_InhNeu_withHar <- FindNeighbors(MEA_InhNeu_withHar, reduction = "harmony", dims = 1:30) 
MEA_InhNeu_withHar <- FindClusters(MEA_InhNeu_withHar, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6), 
                                   algorithm = 1, n.iter = 1000,
                                   graph.name = "SCT_snn")

pdf("output_InhNeu/01_UMAP_Labelled_afterHarmony_Res02.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.2" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_InhNeu/01_UMAP_Labelled_afterHarmony_Res05.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.5" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_InhNeu/01_UMAP_Labelled_afterHarmony_Res04.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.4" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_InhNeu/01_UMAP_Labelled_afterHarmony_Res03.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.3" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_InhNeu/01_UMAP_Labelled_afterHarmony_Res06.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.6" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_InhNeu/01_UMAP_Labelled_afterHarmony_Res07.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.7" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_InhNeu/01_UMAP_Labelled_afterHarmony_Res08.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.8" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()


Idents(MEA_InhNeu_withHar) <- "SCT_snn_res.0.5"
MEA_InhNeu_withHar@meta.data$seurat_clusters <- MEA_InhNeu_withHar@meta.data$SCT_snn_res.0.5

save(MEA_InhNeu_withHar, file = "output_InhNeu/01_SeuratObj_InhNeruon_Reclust_withHarmony_Res05.RData")

all_markers_clustID <- presto::wilcoxauc(MEA_InhNeu_withHar, 'SCT_snn_res.0.5', assay = 'data')

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_InhNeu/Res05_Presto_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")


































##### Relabelling again
setwd("/Users/suganyasubramanian/Haley_MEA_B2/")
load("output_InhNeu/01_SeuratObj_InhNeruon_Reclust_withHarmony_Res04.RData")
#Idents(MEA_InhNeu_withHar) <- "seurat_clusters"

Labels <- read.table("output_InhNeu/Inh_Labels.txt",header=T,sep="\t")
#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$Sub_InhCell)

current.cluster.ids
new.cluster.ids

MEA_InhNeu_withHar@active.ident <- plyr::mapvalues(x = MEA_InhNeu_withHar@active.ident, 
                                                   from = current.cluster.ids, 
                                                   to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
MEA_InhNeu_withHar@meta.data$Sub_InhCell <- MEA_InhNeu_withHar@active.ident

table(MEA_InhNeu_withHar@meta.data$Sub_InhCell)
L3/5_SST L1/3_PVALB   L1/3_VIP   L5/6_SST L5/6_LAMP5   L2/4_VIP   L1_LAMP5 
606        635        987        514        402        238        209 
L1_ADARB2    Glia_NA     L1_VIP    L1_RELN  L1/2_PAX6 L3/6_PVALB   L4/6_SST 
195        183        194        142        105         81         35 

pdf("output_InhNeu/08_VlnPlot_nUMI.pdf", width = 10, height = 8)
VlnPlot(MEA_InhNeu_withHar, features = c("nUMI","pMito"),pt.size = 0, ncol = 1) + theme(legend.position="none")
dev.off()

save(MEA_InhNeu_withHar, file = "output_InhNeu/02_Relabelled_InhNeurons_Res04.RData")

MEA_Inh <- subset(MEA_InhNeu_withHar, subset = Sub_InhCell == "Glia_NA", invert = TRUE)
table(MEA_Inh@meta.data$Sub_InhCell)

pdf("output_InhNeu/09_FinalLabelled.pdf", width = 8, height = 6)
DimPlot(object = MEA_Inh, reduction = "umap",group.by = "Sub_InhCell" ,
        label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-23,18)
dev.off()

pdf("output_InhNeu/09b_UMAP_Labelled_Splitby_Condition.pdf", width = 16, height = 6)
DimPlot(object = MEA_Inh, reduction = "umap", label = TRUE,
        pt.size = 0.5,split.by = "Condition") + theme(legend.position="none") + xlim(-23,18)
dev.off()

order <- c('L1_VIP','L1/3_VIP','L2/4_VIP',
           'L1_RELN','L1/2_PAX6',
           'L1_LAMP5','L5/6_LAMP5',
           'L1/3_PVALB','L3/6_PVALB',
           'L3/5_SST','L4/6_SST','L5/6_SST',
           'L1_ADARB2')


Idents(MEA_Inh) <- factor(Idents(MEA_Inh), levels = order) # reorder the factor based on cell type
MEA_Inh$Sub_InhCell <- factor(MEA_Inh$Sub_InhCell, levels = order)

markers <- c('VIP','RELN','PAX6','LAMP5','PVALB','SST','ADARB2','ELAVL2', 'GAD1', 'GAD2') #"CALM1",

pdf("output_InhNeu/10_DotPlot_Final.pdf", width = 10, height = 6)
DotPlot(MEA_Inh, features = markers) + rotate_x_text(angle = 45)
dev.off()

save(MEA_Inh, file = "output_InhNeu/03_Relabelled_InhNeurons_Res04_NoGlia.RData")

pdf("output_InhNeu/11_VlnPlot_nUMI_NoGlia.pdf", width = 10, height = 4)
VlnPlot(MEA_Inh, features = c("nUMI","pMito"),pt.size = 0) + theme(legend.position="none")
dev.off()

#color <- c("#17154f","#b0799a","#e69b00","#355828","#4f2400","#b62684","#5eb0e5")
#"#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")
color <- c("#92a8d1","#0e3386",
           "#c071d0","#5c2893",
           "#db821d","#890c0a",
           "#a6ea7c","#28935c",
           "#b0b0b2","#4c3c31",
           "#e2fd07","#7f6000",
           "#fddad2","#f6471c",
           "#cdf5ff","#05cefe")


prop_Cell <-  MEA_Inh@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(Sub_InhCell) %>%
  dplyr::group_by(Sub_InhCell,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(Sub_InhCell) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(Sub_InhCell) %>%
  ggbarplot("Sub_InhCell", "percent",
            fill = "Genotype", color = "Genotype", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_InhNeu/12_Sub_InhCellProportion_By_Genotype.pdf", plot = prop_Cell, width = 18, height = 15, units = "in", dpi = 150)

all_markers_clustID <- presto::wilcoxauc(MEA_Inh, 'Sub_InhCell', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_InhNeu/Sub_InhCell-Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_InhNeu/Sub_InhCell-PrestoTop10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")




