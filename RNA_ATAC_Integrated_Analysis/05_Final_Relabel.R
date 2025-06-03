suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(Matrix)
  library(tidyverse)
  library(ggplot2)
  library(ggpubr)
  library(sctransform)
  library(ggrastr)
  library(clustree)
  library(cowplot)
  library(anndata)
})

#Relabel the clusters with Final annotations

load("Integrate_RNA_ATAC/RNA_ATAC_07_Prediction_Labelled.RData")


Idents(Haley_MEA_Int) <- "seurat_clusters"

Labels <- read.table("Integrate_RNA_ATAC/Label_Final.txt",header=T,sep="\t")

#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$Cell)

current.cluster.ids
new.cluster.ids


##Step 3: Rename cluster with new Labels
Haley_MEA_Int@active.ident <- plyr::mapvalues(x = Haley_MEA_Int@active.ident, 
                                              from = current.cluster.ids, 
                                              to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
Haley_MEA_Int@meta.data$Cell <- Haley_MEA_Int@active.ident


Haley_MEA_Int@meta.data <- Haley_MEA_Int@meta.data %>%
  mutate(Cell_Class = case_when(grepl("Oligo", Cell) ~ "Oligo",
                                grepl("OPC", Cell) ~ "OPC",
                                grepl("Astro", Cell) ~ "Astro",
                                grepl("Micro", Cell) ~ "Micro",
                                grepl("Exc", Cell) ~ "Exc_Neurons",
                                grepl("Inh", Cell) ~ "Inh_Neurons",
                                #grepl("Peri", Cell) ~ "Pericytes",
                                #grepl("Epen", Cell) ~ "Ependymal",
                                grepl("Peri|Endo|VLMC|VECA|VECV|ABC|Epen", Cell) ~ "Vascular",
                                .default = "NoClass")) %>%
  mutate(sub_Class = case_when(grepl("Oligo", Cell) ~ "Glia",
                                grepl("OPC", Cell) ~ "Glia",
                                grepl("Astro", Cell) ~ "Glia",
                                grepl("Micro", Cell) ~ "Glia",
                                grepl("Exc", Cell) ~ "Exc_Neurons",
                                grepl("Inh", Cell) ~ "Inh_Neurons",
                                #grepl("Peri", Cell) ~ "Pericytes",
                                #grepl("Epen", Cell) ~ "Ependymal",
                                grepl("Peri|Endo|VLMC|VECA|VECV|ABC|Epen", Cell) ~ "Glia",
                                .default = "NoClass")) %>%
  mutate(Super_Class = case_when(grepl("Oligo", Cell) ~ "NonNeurons",
                               grepl("OPC", Cell) ~ "NonNeurons",
                               grepl("Astro", Cell) ~ "NonNeurons",
                               grepl("Micro", Cell) ~ "NonNeurons",
                               grepl("Exc", Cell) ~ "Neurons",
                               grepl("Inh", Cell) ~ "Neurons",
                               #grepl("Peri", Cell) ~ "Pericytes",
                               #grepl("Epen", Cell) ~ "Ependymal",
                               grepl("Peri|Endo|VLMC|VECA|VECV|ABC|Epen", Cell) ~ "NonNeurons",
                               .default = "NoClass"))

pdf("Integrate_RNA_ATAC/13_UMAP_FinalLabelled_cell.pdf", width = 8, height = 6)
DimPlot(object = Haley_MEA_Int, reduction = "umap", label = TRUE, pt.size = 0.5 ) + xlim(-15,15) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/14_UMAP_FinalLabelled_Cell_Class.pdf", width = 8, height = 6)
DimPlot(object = Haley_MEA_Int, group.by = "Cell_Class",reduction = "umap", label = TRUE, pt.size = 0.5 ) + xlim(-15,15) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/15_UMAP_FinalLabelled_sub_Class.pdf", width = 8, height = 6)
DimPlot(object = Haley_MEA_Int, group.by = "sub_Class",reduction = "umap", label = TRUE, pt.size = 0.5 ) + xlim(-15,15) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/16_UMAP_FinalLabelled_Super_Class.pdf", width = 8, height = 6)
DimPlot(object = Haley_MEA_Int, group.by = "Super_Class",reduction = "umap", label = TRUE, pt.size = 0.5 ) + xlim(-15,15) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/17_VlnPlot_Final_Cell.pdf", width = 15, height = 15)
VlnPlot(Haley_MEA_Int, features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment", "nucleosome_signal","NF"),
        ncol = 3, group.by = "Cell",
        pt.size = 0) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/17_VlnPlot_Final_PredictionCell.pdf", width = 15, height = 15)
VlnPlot(Haley_MEA_Int, features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment", "nucleosome_signal","NF"),
        ncol = 3, group.by = "Prediction_Cell",
        pt.size = 0) + theme(legend.position="none")
dev.off()

order <- c('Astro_1','Micro_1','Micro_2','Oligo_1','OPC_1',
           'Exc_1','Exc_2','Exc_3','Exc_4','Exc_5',
           'Inh_1','Inh_2','Inh_3','Inh_4',
           'Endo_1','Endo_2','VLMC_1')

Idents(Haley_MEA_Int) <- factor(Idents(Haley_MEA_Int), levels = order) # reorder the factor based on cell type
Haley_MEA_Int$Cell <- factor(Haley_MEA_Int$Cell, levels = order)


markers <- c('GJA1','DOCK8','MBP','PCDH15',
             'SLC17A7','GAD1',
             'FLT1','CPED1')

pdf("Integrate_RNA_ATAC/18_DotPlot_Final_Cell.pdf", width = 10, height = 6)
DotPlot(Haley_MEA_Int, features = markers) #+ rotate_x_text(angle = 45)
dev.off()

save(Haley_MEA_Int, file="Integrate_RNA_ATAC/RNA_ATAC_08_Final_Labelled.RData")
