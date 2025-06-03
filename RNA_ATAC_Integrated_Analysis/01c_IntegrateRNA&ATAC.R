suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(Matrix)
  library(tidyverse)
  library(ggplot2)
  library(sctransform)
  library(hdf5r) #install.packages("hdf5r")
  library(ggrastr)
  library(clustree)
  library(cowplot)
  library(ggpubr)
  library(future)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
})

setwd("/zfs/musc3/SS/Haley_MEA_B2/")

load("IntegrateATAC/ATAC_01_Integrated.RData") #integrated
load("Integrate_RNA/03_seuObject_Integrated_withNF.RData") #seuObject_integrated

#Find matching column names
rna_names<-colnames(seuObject_integrated) #merged seurat object from raw.h5
length(rna_names) #95879

atac_names<-colnames(integrated)
length(atac_names) #77759

intersect<-intersect(atac_names, rna_names)
length(intersect) #65032

subrna <- seuObject_integrated[,intersect]

subatac <- integrated[,intersect]

#combine seurat objects
subrna[['ATAC']] <- subatac[['ATAC']]

save(subrna, file = "Integrate_RNA_ATAC/RNA_04_subRNA_Intersect.RData")
save(subatac, file = "Integrate_RNA_ATAC/ATAC_02_subATAC_Intersect.RData")

#### Integrate signac atac seurat RNA

#recluster ATAC
DefaultAssay(subrna) <- "ATAC"
subrna <- RunTFIDF(subrna)
subrna <- FindTopFeatures(subrna)
subrna <- RunSVD(subrna)
subrna <- RunUMAP(subrna, dims = 2:50, reduction = 'lsi',
                  reduction.name = "umap.atac", reduction.key = "atacUMAP_")

# build a joint neighbor graph using both assays
subrna <- FindMultiModalNeighbors(
  object = subrna,
  reduction.list = list("pca", "lsi"),
  dims.list = list(1:50, 2:50),
  modality.weight.name = "RNA.weight",
  verbose = TRUE
)

subrna <- RunUMAP(
  object = subrna,
  #reduction.name = "int_umap",
  reduction.name = "wnn.umap", reduction.key = "wnnUMAP_",
  nn.name = "weighted.nn",
  assay = "RNA",
  verbose = TRUE
)

subrna <- FindClusters(
  object = subrna, 
  resolution = c(0.3,0.4,0.5,0.6,0.7,0.8), 
  graph.name = "wsnn",
  algorithm = 3,
  verbose = TRUE
)


## Transfer ATAC metadata to rna object
metaATAC <- as.data.frame(subatac@meta.data)
subrna <- AddMetaData(object = subrna, 
                      metadata = metaATAC)

save(subrna, file = "Integrate_RNA_ATAC/RNA_ATAC_01_Integrated.RData")


pdf("Integrate_RNA_ATAC/RNA_ATAC_01_UMAP_IntegratedAssay.pdf", width=6,height=5)
DimPlot(subrna, label = TRUE, group.by = "wsnn_res.0.5",repel = TRUE, reduction = "wnn.umap") + NoLegend()
dev.off()

p1 <- DimPlot(subrna, reduction = "umap",group.by = "integrated_snn_res.0.5", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(subrna, reduction = "umap.atac",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(subrna, reduction = "wnn.umap",group.by = "wsnn_res.0.5", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

pdf("Integrate_RNA_ATAC/RNA_ATAC_02_UMAP_Allassay.pdf", width=20,height=5)
plot_grid(p1,p2,p3,ncol=3)
dev.off()