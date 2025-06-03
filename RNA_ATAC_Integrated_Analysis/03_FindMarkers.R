suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  library(speckle)
  library(magrittr)
  library(broom)
  library(muscat)
  library(Seurat)
  library(clustree)
  #library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(garnett)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(monocle)
})

load("Integrate_RNA_ATAC/RNA_ATAC_05_Final_Integrated.RData")


############ FIND MARKERS #############
Haley_MEA_Int$seurat_clusters <- Haley_MEA_Int$SCT_snn_res.0.5
Idents(Haley_MEA_Int) <- "seurat_clusters"

all_markers_clustID <- presto::wilcoxauc(Haley_MEA_Int, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "Integrate_RNA_ATAC/Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

## top 20 markers for each group
## filter for nominally significant (p<0.05) and over-expressed (auc>0.5)
top20 <- presto::top_markers(all_markers.Sign,n = 20,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top20, 
                     file = "Integrate_RNA_ATAC/Presto_Top20.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")