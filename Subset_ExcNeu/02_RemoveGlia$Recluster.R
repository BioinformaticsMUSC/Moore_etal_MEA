suppressPackageStartupMessages({
  library(WGCNA)
  library(cluster)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(ggpubr)
  library(Seurat)
})
setwd("/Users/suganyasubramanian/Haley_IVS_RNA/")
load("output_ExcNeu/01_SeuratObj_ExcNeruon_Reclust_withHarmony_Res04.RData")

### Removing Glia cell clusters
Idents(MEA_ExcNeu_withHar) <- 'seurat_clusters'

MEA_ExcNeu_noglia <- subset(MEA_ExcNeu_withHar, subset = seurat_clusters %in% c('4','8','13','16'), invert = TRUE)

MEA_ExcNeu_noglia <- processing_seurat_sctransform(MEA_Exc, 
                                                   vars_to_regress = c("nUMI","pMito","pRibo",
                                                                       "Age","Sex","Hemis","EpDur","SeqBatch"),
                                                   npcs = 50, 
                                                   res = c(0.3,0.4,0.5,0.6))

DefaultAssay(MEA_ExcNeu_noglia) <- "RNA"
MEA_ExcNeu_noglia <- NormalizeData(object = MEA_ExcNeu_noglia, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)

library(harmony)
MEA_ExcNeu_withHar_noglia <- MEA_ExcNeu_noglia

MEA_ExcNeu_withHar_noglia <- RunHarmony(MEA_ExcNeu_withHar_noglia, assay.use="RNA", group.by.vars = "Genotype")
Reductions(MEA_ExcNeu_withHar_noglia)

MEA_ExcNeu_withHar_noglia <- RunUMAP(MEA_ExcNeu_withHar_noglia, reduction = "harmony", dims = 1:30)
MEA_ExcNeu_withHar_noglia <- FindNeighbors(MEA_ExcNeu_withHar_noglia, reduction = "harmony", dims = 1:30) 
MEA_ExcNeu_withHar_noglia <- FindClusters(MEA_ExcNeu_withHar_noglia, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6), 
                                          algorithm = 1, n.iter = 1000,
                                          graph.name = "SCT_snn")

save(MEA_ExcNeu_withHar_noglia, file = 'output_ExcNeu/02_SeuObj_ExcNeurons_NoGlia_ReclusteredwithHarmony.RData')

pdf("output_ExcNeu/03_FinalLabelled_check_NoGlia.pdf", width = 12, height = 6)
p1 <- DimPlot(object = MEA_ExcNeu_withHar_noglia, reduction = "umap",group.by = "Sub_ExcCell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p2 <- DimPlot(object = MEA_ExcNeu_withHar_noglia, reduction = "umap",group.by = "SCT_snn_res.0.5" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
plot_grid(p1,p2, ncol = 2)
dev.off()