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

load("Haley_MEA_Signac_P6/RNA_ATAC_08_Final_Labelled.RData")

source("Utils.R")

MEA_ExcNeu <- subset(Haley_MEA_Int, 
                     subset = Cell_Class == "Exc_Neurons")

MEA_ExcNeu <- processing_seurat_sctransform(MEA_ExcNeu, 
                                            vars_to_regress = c("nUMI","pMito","pRibo",
                                                                "Age","Sex","Hemis","EpDur","SeqBatch"),
                                            npcs = 50, 
                                            res = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))

DefaultAssay(MEA_ExcNeu) <- "RNA"
MEA_ExcNeu <- NormalizeData(object = MEA_ExcNeu, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)

dir.create("output_ExcNeu")

pdf("output_ExcNeu/01_UMAP_Labelled.pdf", width = 10, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu, reduction = "umap",group.by = "Cell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_ExcNeu, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

library(harmony)
MEA_ExcNeu_withHar <- MEA_ExcNeu

MEA_ExcNeu_withHar <- RunHarmony(MEA_ExcNeu_withHar, assay.use="RNA", group.by.vars = "Genotype")
Reductions(MEA_ExcNeu_withHar)

MEA_ExcNeu_withHar <- RunUMAP(MEA_ExcNeu_withHar, reduction = "harmony", dims = 1:30)
MEA_ExcNeu_withHar <- FindNeighbors(MEA_ExcNeu_withHar, reduction = "harmony", dims = 1:30) 
MEA_ExcNeu_withHar <- FindClusters(MEA_ExcNeu_withHar, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6), 
                                   algorithm = 1, n.iter = 1000,
                                   graph.name = "SCT_snn")

pdf("output_ExcNeu/01_UMAP_Labelled_afterHarmony_Res02.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.2" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_ExcNeu/01_UMAP_Labelled_afterHarmony_Res05.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.5" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_ExcNeu/01_UMAP_Labelled_afterHarmony_Res04.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.4" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_ExcNeu/01_UMAP_Labelled_afterHarmony_Res03.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.3" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_ExcNeu/01_UMAP_Labelled_afterHarmony_Res06.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.6" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_ExcNeu/01_UMAP_Labelled_afterHarmony_Res07.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.7" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_ExcNeu/01_UMAP_Labelled_afterHarmony_Res08.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_ExcNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.8" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()


Idents(MEA_ExcNeu_withHar) <- "SCT_snn_res.0.4"
MEA_ExcNeu_withHar@meta.data$seurat_clusters <- MEA_ExcNeu_withHar@meta.data$SCT_snn_res.0.4

save(MEA_ExcNeu_withHar, file = "output_ExcNeu/01_SeuratObj_ExcNeruon_Reclust_withHarmony_Res04.RData")

all_markers_clustID <- presto::wilcoxauc(MEA_ExcNeu_withHar, 'SCT_snn_res.0.4', assay = 'data')

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_ExcNeu/Res04_Presto_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")