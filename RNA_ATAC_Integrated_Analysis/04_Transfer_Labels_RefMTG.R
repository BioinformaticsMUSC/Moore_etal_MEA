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

##Load reference and add meta data
load("AllenRef/Allen_MultiRegion_SeuObjWithAnno.RData") #seuRef

load("Integrate_RNA_ATAC/RNA_ATAC_05_Final_Integrated.RData") #Haley_MEA_Int

seuQuery <- Haley_MEA_Int

# perform standard preprocessing on each object
seuRef <- NormalizeData(seuRef)
seuRef <- FindVariableFeatures(seuRef)
seuRef <- ScaleData(seuRef)

seuQuery <- NormalizeData(seuQuery)
seuQuery <- FindVariableFeatures(seuQuery)
seuQuery <- ScaleData(seuQuery)

# find anchors
anchors <- FindTransferAnchors(reference = seuRef, query = seuQuery)
# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = seuRef$cell_type_alias_label)

save(anchors,predictions, file = "Integrate_RNA_ATAC/RNA_ATAC_06_LabelTransferAnchor.RData")

seuQuery <- AddMetaData(object = seuQuery, metadata = predictions)
save(seuQuery, file = "Integrate_RNA_ATAC/RNA_ATAC_06_SeuQuery_withPredictions.RData")

Haley_MEA_Int <- AddMetaData(object = Haley_MEA_Int, 
                             metadata = predictions$predicted.id,
                             col.name = 'predicted.id')


Haley_MEA_Int@meta.data <- Haley_MEA_Int@meta.data %>%
  mutate(predicted.class = sapply(X = strsplit(Haley_MEA_Int$predicted.id, split = " "), FUN = "[", 1)) %>%
  mutate(predicted.subclass = sapply(X = strsplit(Haley_MEA_Int$predicted.id, split = " "), FUN = "[", 2)) %>%
  mutate(predicted.marker1 = sapply(X = strsplit(Haley_MEA_Int$predicted.id, split = " "), FUN = "[", 3)) %>%
  mutate(predicted.marker2 = sapply(X = strsplit(Haley_MEA_Int$predicted.id, split = " "), FUN = "[", 4)) %>%
  unite(predicted.label,predicted.class, predicted.subclass,predicted.marker1,sep = "_", remove = FALSE) %>%
  unite(predicted.label2,predicted.class, predicted.subclass,sep = "_", remove = FALSE)

save(Haley_MEA_Int, file="Integrate_RNA_ATAC/RNA_ATAC_06_Final_Integrated_LabelTransfer.RData")

pdf("Integrate_RNA_ATAC/10_UMAP_LAble_Transfer_name.pdf", width = 6, height = 4)
DimPlot(object = Haley_MEA_Int, group.by = "predicted.id",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/10_UMAP_LAble_Transfer_class.pdf", width = 6, height = 4)
DimPlot(object = Haley_MEA_Int, group.by = "predicted.class",reduction = "umap", label = FALSE, pt.size = 0.5) #+ theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/10_UMAP_LAble_Transfer_subclass.pdf", width = 6, height = 4)
DimPlot(object = Haley_MEA_Int, group.by = "predicted.subclass",reduction = "umap", label = FALSE, pt.size = 0.5) #+ theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/10_UMAP_LAble_Transfer_label.pdf", width = 6, height = 4)
DimPlot(object = Haley_MEA_Int, group.by = "predicted.label",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/10_UMAP_LAble_Transfer_label2.pdf", width = 6, height = 4)
DimPlot(object = Haley_MEA_Int, group.by = "predicted.label2",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()


pdf("Integrate_RNA_ATAC/10_UMAP_LAble_Transfer_class_splitby.pdf", width = 20, height = 15)
DimPlot(object = Haley_MEA_Int,split.by= "predicted.class",group.by = "predicted.label",
        reduction = "umap", label = TRUE, ncol = 4, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/10_UMAP_LAble_Transfer_label_splitby.pdf", width = 20, height = 40)
DimPlot(object = Haley_MEA_Int,split.by= "predicted.label2",group.by = "predicted.label",
        reduction = "umap", label = TRUE, ncol = 4, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("Integrate_RNA_ATAC/10_UMAP_LAble_Transfer_cluster_splitby.pdf", width = 30, height = 40)
DimPlot(object = Haley_MEA_Int,split.by= "SCT_snn_res.0.5",group.by = "predicted.label2",
        reduction = "umap", label = TRUE, ncol = 6, pt.size = 0.5) + theme(legend.position="none")
dev.off()

cnt <- table(Haley_MEA_Int$predicted.label,Haley_MEA_Int$SCT_snn_res.0.5)
head(cnt)
openxlsx::write.xlsx(cnt,file = "Integrate_RNA_ATAC/Annotation_Count.xlsx", colNames = TRUE, borders = "columns")


