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
  library(scDblFinder)
  library(DropletQC)
})

load("Integrate_RNA_ATAC/RNA_ATAC_01_Integrated.RData")

#### Predict droplet QC cell status
df_nf <- subrna
df_nf <- as.data.frame(df_nf@meta.data) %>% dplyr::select(NF,nUMI)
df_nf <- identify_empty_drops(nf_umi = df_nf)
subrna@meta.data$nf_cell_status <- df_nf$cell_status

table(subrna$nf_cell_status)
#cell empty_droplet 
#63122          1910

pdf("Integrate_RNA_ATAC/RNA_ATAC_04a_QC_UMAP_NF_CellStatus.pdf", width = 6, height = 5)
DimPlot(subrna,
        reduction = "umap",raster=FALSE,
        group.by= "nf_cell_status")
dev.off()

DefaultAssay(subrna) <- "ATAC"

pdf("Integrate_RNA_ATAC/RNA_ATAC_03_DensityPlot_TssEnrichment.pdf", width=5,height=5)
DensityScatter(subrna, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

pdf("Integrate_RNA_ATAC/RNA_ATAC_04_QC_plot.pdf", width=18,height=10)
VlnPlot(
  object = subrna,
  features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment", "nucleosome_signal","NF"),
  ncol = 3, group.by = "Genotype",
  pt.size = 0)
dev.off()

# filter out low quality cells
subrna_filt <- subset(
  x = subrna,
  subset = nUMI < 25000 & pMito < 5 & nGene > 250 &
    nCount_ATAC < 180000 & nucleosome_signal < 2 & TSS.enrichment > 1 &
    NF > 0.5 & nf_cell_status == "cell")

pdf("Integrate_RNA_ATAC/RNA_ATAC_05_QC_plot_Filtered.pdf", width=18,height=10)
VlnPlot(
  object = subrna_filt,
  features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment", "nucleosome_signal","NF"),
  ncol = 3, group.by = "Genotype",
  pt.size = 0)
dev.off()

pdf("Integrate_RNA_ATAC/RNA_ATAC_06_UMAP_QC_Filtered.pdf", width=6,height=5)
DimPlot(subrna_filt, label = TRUE, group.by = "wsnn_res.0.5",repel = TRUE, reduction = "wnn.umap") + NoLegend()
dev.off()

### find doublets
#Convert into Single Cell Experiment object
sce <- as.SingleCellExperiment(subrna_filt)

#Find Doublets by genotype
sce0 <- scDblFinder(sce,
                    samples="Genotype",
                    #BPPARAM=MulticoreParam(3),
                    nfeatures = 3000,
                    dims = 30,
                    dbr.sd = 0,
                    multiSampleMode="split")

# add in the new column "Doublets" to make if it doublet or singlet
subrna_filt@meta.data$Doublets <- sce0$scDblFinder.class

save(subrna_filt, file = "Integrate_RNA_ATAC/RNA_ATAC_02_Integrated_QCFiltered.RData")

subrna_filt_nodoub <- subset(subrna_filt, subset = Doublets == "singlet")

#RunUMAP again on filtered data
source("Utils.R")

DefaultAssay(subrna_filt_nodoub) <- "RNA"

subrna_filt_nodoub <- processing_seurat_sctransform(subrna_filt_nodoub, 
                             vars_to_regress = c("nCount_RNA","pMito","pRibo",
                                      "Age","Sex","Hemis","EpDur","SeqBatch"), 
                              npcs = 50, 
                              res = c(0.4,0.5,0.6))

DefaultAssay(subrna_filt_nodoub) <- "RNA"

subrna_filt_nodoub <- NormalizeData(object = subrna_filt_nodoub, 
                                        normalization.method = "LogNormalize", 
                                        scale.factor = 10000)

pdf("Integrate_RNA_ATAC/RNA_ATAC_07_UMAP_Reclustered_Res05.pdf", width = 10, height = 6)
p1 <- DimPlot(object = subrna_filt_nodoub, group.by = "SCT_snn_res.0.5",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = subrna_filt_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

library(harmony)
Haley_MEA <- RunHarmony(subrna_filt_nodoub, assay.use="RNA", 
                       group.by.vars = c("Batch","Genotype"))
Haley_MEA <- RunUMAP(Haley_MEA, reduction = "harmony", dims = 1:50)
Haley_MEA <- FindNeighbors(Haley_MEA, reduction = "harmony", dims = 1:50) 
Haley_MEA <- FindClusters(Haley_MEA, resolution = c(0.3,0.4,0.5,0.6,0.7,0.8), 
                          algorithm = 1, n.iter = 1000,graph.name = "SCT_snn")

pdf("Integrate_RNA_ATAC/RNA_ATAC_08_UMAP_NoDlbt_Reclustd_withHarmony_Res0.3.pdf", width = 10, height = 6)
p1 <- DimPlot(object = Haley_MEA, group.by = "SCT_snn_res.0.3",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = Haley_MEA, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Integrate_RNA_ATAC/RNA_ATAC_08_UMAP_NoDlbt_Reclustd_withHarmony_Res0.4.pdf", width = 10, height = 6)
p1 <- DimPlot(object = Haley_MEA, group.by = "SCT_snn_res.0.4",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = Haley_MEA, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Integrate_RNA_ATAC/RNA_ATAC_08_UMAP_NoDlbt_Reclustd_withHarmony_Res0.5.pdf", width = 10, height = 6)
p1 <- DimPlot(object = Haley_MEA, group.by = "SCT_snn_res.0.5",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = Haley_MEA, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Integrate_RNA_ATAC/RNA_ATAC_08_UMAP_NoDlbt_Reclustd_withHarmony_Res0.6.pdf", width = 10, height = 6)
p1 <- DimPlot(object = Haley_MEA, group.by = "SCT_snn_res.0.6",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = Haley_MEA, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Integrate_RNA_ATAC/RNA_ATAC_08_UMAP_NoDlbt_Reclustd_withHarmony_Res0.7.pdf", width = 10, height = 6)
p1 <- DimPlot(object = Haley_MEA, group.by = "SCT_snn_res.0.7",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = Haley_MEA, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Integrate_RNA_ATAC/RNA_ATAC_08_UMAP_NoDlbt_Reclustd_withHarmony_Res0.8.pdf", width = 10, height = 6)
p1 <- DimPlot(object = Haley_MEA, group.by = "SCT_snn_res.0.8",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = Haley_MEA, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

Haley_MEA$seurat_clusters <- Haley_MEA$SCT_snn_res.0.5
Idents(Haley_MEA) <- "seurat_clusters"

pdf("Integrate_RNA_ATAC/RNA_ATAC_09_UMAP_NoDlbt_Reclustd_withHarmony_Final_Res05.pdf", width = 10, height = 6)
p1 <- DimPlot(object = Haley_MEA, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = Haley_MEA, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

Haley_MEA_Int <- FindMultiModalNeighbors(Haley_MEA, 
                                         reduction.list = list("harmony", "lsi"),
                                         dims.list = list(1:50, 2:50),
                                         modality.weight.name = "RNA.weight",
                                         verbose = TRUE)

Haley_MEA_Int <- RunUMAP(Haley_MEA_Int, 
                         nn.name = "weighted.nn", 
                         reduction.name = "wnn.umap", 
                         reduction.key = "wnnUMAP_", 
                         min.dist = 0.27, n.neighbors = 50)

#find clusters
Haley_MEA_Int <- FindClusters(
  object = Haley_MEA_Int, 
  resolution = c(0.3,0.4,0.5,0.6,0.7,0.8), 
  graph.name = "wsnn",
  algorithm = 3,
  verbose = TRUE
)

colnames(Haley_MEA_Int@meta.data)

p1 <- DimPlot(Haley_MEA_Int, reduction = "umap", group.by = "SCT_snn_res.0.5",label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
p2 <- DimPlot(Haley_MEA_Int, reduction = "umap.atac",  label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(Haley_MEA_Int, reduction = "wnn.umap",group.by = "wsnn_res.0.6", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")

p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

pdf("Integrate_RNA_ATAC/RNA_ATAC_10_UMAP_RNA_ATAC_WNN_AllAssay.pdf", width=20,height=5)
plot_grid(p1,p2,p3,ncol=3)
dev.off()

pdf("Integrate_RNA_ATAC/RNA_ATAC_11_QC_plot_byClusters.pdf", width=18,height=15)
VlnPlot(
  object = Haley_MEA_Int,
  features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment","nucleosome_signal","NF"),
  ncol = 3, group.by = "SCT_snn_res.0.5",
  pt.size = 0)
dev.off()

Haley_MEA_Int <- CellCycleScoring(Haley_MEA_Int,
                                  g2m.features = cc.genes$g2m.genes,
                                  s.features = cc.genes$s.genes)

pdf("Integrate_RNA_ATAC/RNA_ATAC_12_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(Haley_MEA_Int,
        reduction = "umap",raster=FALSE,
        group.by= "Phase")
dev.off()

#save(subrna_filt, file = "Integrate_RNA_ATAC/RNA_ATAC_02_Integrated_QCFiltered.RData")
save(subrna_filt_nodoub, file = "Integrate_RNA_ATAC/RNA_ATAC_03_NoDlbt_Reclusterd.RData")
save(Haley_MEA, file = "Integrate_RNA_ATAC/RNA_ATAC_04_NoDlbt_Reclustd_withHarmony_Res05.RData")
save(Haley_MEA_Int, file = "Integrate_RNA_ATAC/RNA_ATAC_05_Final_Integrated.RData")





