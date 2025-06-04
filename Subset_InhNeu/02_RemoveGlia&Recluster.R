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
load("output_InhNeu/01_SeuratObj_InhNeruon_Reclust_withHarmony_Res05.RData")

### Removing Glia cell clusters
Idents(MEA_InhNeu_withHar) <- 'seurat_clusters'

MEA_InhNeu_noglia <- subset(MEA_InhNeu_withHar, subset = seurat_clusters %in% c('6','8','10','13'), invert = TRUE)

MEA_InhNeu_noglia <- processing_seurat_sctransform(MEA_Inh, 
                                                   vars_to_regress = c("nUMI","pMito","pRibo",
                                                                       "Age","Sex","Hemis","EpDur","SeqBatch"),
                                                   npcs = 50, 
                                                   res = c(0.3,0.4,0.5,0.6))

DefaultAssay(MEA_InhNeu_noglia) <- "RNA"
MEA_InhNeu_noglia <- NormalizeData(object = MEA_InhNeu_noglia, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)

library(harmony)
MEA_InhNeu_withHar_noglia <- MEA_InhNeu_noglia

MEA_InhNeu_withHar_noglia <- RunHarmony(MEA_InhNeu_withHar_noglia, assay.use="RNA", group.by.vars = "Genotype")
Reductions(MEA_InhNeu_withHar_noglia)

MEA_InhNeu_withHar_noglia <- RunUMAP(MEA_InhNeu_withHar_noglia, reduction = "harmony", dims = 1:30)
MEA_InhNeu_withHar_noglia <- FindNeighbors(MEA_InhNeu_withHar_noglia, reduction = "harmony", dims = 1:30) 
MEA_InhNeu_withHar_noglia <- FindClusters(MEA_InhNeu_withHar_noglia, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6), 
                                          algorithm = 1, n.iter = 1000,
                                          graph.name = "SCT_snn")

save(MEA_InhNeu_withHar_noglia, file = 'output_InhNeu/02_SeuObj_InhNeurons_NoGlia_ReclusteredwithHarmony.RData')

pdf("output_InhNeu/03_FinalLabelled_check_NoGlia.pdf", width = 12, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar_noglia, reduction = "umap",group.by = "Sub_InhCell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p2 <- DimPlot(object = MEA_InhNeu_withHar_noglia, reduction = "umap",group.by = "SCT_snn_res.0.5" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
plot_grid(p1,p2, ncol = 2)
dev.off()