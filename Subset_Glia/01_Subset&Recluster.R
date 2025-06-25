suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
})

load("Haley_MEA_Signac_P6/RNA_ATAC_08_Final_Labelled.RData")

source("Utils.R")

MEA_Glia <- subset(Haley_MEA_Int, 
                     subset = Cell_Class == "Glia")

MEA_Glia <- processing_seurat_sctransform(MEA_Glia, 
                                            vars_to_regress = c("nUMI","pMito","pRibo",
                                                                "Age","Sex","Hemis","EpDur","SeqBatch"),
                                            npcs = 50, 
                                            res = c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8))

DefaultAssay(MEA_Glia) <- "RNA"

MEA_Glia <- NormalizeData(object = MEA_Glia, 
                            normalization.method = "LogNormalize", 
                            scale.factor = 10000)

dir.create("output_Glia")

pdf("output_Glia/01_UMAP_Labelled.pdf", width = 10, height = 6)
p1 <- DimPlot(object = MEA_Glia, reduction = "umap",group.by = "Cell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_Glia, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

library(harmony)
MEA_Glia_withHar <- MEA_Glia

MEA_Glia_withHar <- RunHarmony(MEA_Glia_withHar, assay.use="RNA", group.by.vars = "Genotype")
Reductions(MEA_Glia_withHar)

MEA_Glia_withHar <- RunUMAP(MEA_Glia_withHar, reduction = "harmony", dims = 1:30)
MEA_Glia_withHar <- FindNeighbors(MEA_Glia_withHar, reduction = "harmony", dims = 1:30) 
MEA_Glia_withHar <- FindClusters(MEA_Glia_withHar, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6), 
                                   algorithm = 1, n.iter = 1000,
                                   graph.name = "SCT_snn")

pdf("output_Glia/01_UMAP_Labelled_afterHarmony_Res02.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "SCT_snn_res.0.2" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_Glia/01_UMAP_Labelled_afterHarmony_Res05.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "SCT_snn_res.0.5" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_Glia/01_UMAP_Labelled_afterHarmony_Res04.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "SCT_snn_res.0.4" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_Glia/01_UMAP_Labelled_afterHarmony_Res03.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "SCT_snn_res.0.3" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_Glia/01_UMAP_Labelled_afterHarmony_Res06.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "SCT_snn_res.0.6" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_Glia/01_UMAP_Labelled_afterHarmony_Res07.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "SCT_snn_res.0.7" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

pdf("output_Glia/01_UMAP_Labelled_afterHarmony_Res08.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "Genotype" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_Glia_withHar, reduction = "umap",group.by = "SCT_snn_res.0.8" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()


Idents(MEA_Glia_withHar) <- "SCT_snn_res.0.5"
MEA_Glia_withHar@meta.data$seurat_clusters <- MEA_Glia_withHar@meta.data$SCT_snn_res.0.5

save(MEA_Glia_withHar, file = "output_Glia/01_SeuratObj_Glia_Reclust_withHarmony_Res05.RData")

all_markers_clustID <- presto::wilcoxauc(MEA_Glia_withHar, 'SCT_snn_res.0.5', assay = 'data')

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_Glia/Res05_Presto_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")