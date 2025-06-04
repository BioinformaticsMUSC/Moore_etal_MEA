module load r/4.2.3
module load hdf5/1.14.3
module load gsl/2.7.1
module load gcc/13.2.0

R

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

load("output_InhNeu/01_SeuratObj_InhNeruon_Reclust_withHarmony_Res05.RData")
load("/scratch/subramanian/Haley_MEA/AllenRef/Allen_MultiRegion_SeuObjWithAnno.RData")

table(seuRef$subclass_label)
Astrocyte     Endothelial              IT           L4 IT           L5 ET 
1188              70           21884            3725             158 
L5/6 IT Car3         L5/6 NP           L6 CT             L6b           LAMP5 
1055             818            2562            1083            2438 
Microglia Oligodendrocyte             OPC            PAX6        Pericyte 
750            1932             774             325              32 
PVALB             SST             VIP            VLMC 
2805            2361            3538              11

table(seuRef$class_label) #GABAergic - Inhibitory Glutamatergic - Excitatory
#GABAergic Glutamatergic  Non-neuronal 
#11467         31285          4757

Ref_Inc <- subset(seuRef, subset = class_label == "GABAergic")
table(Ref_Inc$cell_type_alias_label)
#Inh L1 ADARB2 ADAM33      Inh L1 ADARB2 DISP2       Inh L1 LAMP5 GGT8P


# perform standard preprocessing on each object
Ref_Inc <- NormalizeData(Ref_Inc)
Ref_Inc <- FindVariableFeatures(Ref_Inc)
Ref_Inc <- ScaleData(Ref_Inc)

save(Ref_Inc, file = "AllenRef/SubsetRef_Inhibitory.RData")

MEA_InhNeu_withHar <- NormalizeData(MEA_InhNeu_withHar)
MEA_InhNeu_withHar <- FindVariableFeatures(MEA_InhNeu_withHar)
MEA_InhNeu_withHar <- ScaleData(MEA_InhNeu_withHar)

# find anchors
anchors <- FindTransferAnchors(reference = Ref_Inc, query = MEA_InhNeu_withHar)
# transfer labels
predictions <- TransferData(anchorset = anchors, refdata = Ref_Inc$cell_type_alias_label)
save(anchors,predictions, file = "output_InhNeu/00_LabelTransferAnchor.RData")

MEA_InhNeu_withHar <- AddMetaData(object = MEA_InhNeu_withHar, 
                                  metadata = predictions$predicted.id,
                                  col.name = 'predicted.inhid')

MEA_InhNeu_withHar
An object of class Seurat 
291559 features across 4252 samples within 4 assays 
Active assay: RNA (20017 features, 2000 variable features)
3 other assays present: SCT, integrated, ATAC
6 dimensional reductions calculated: pca, umap, lsi, umap.atac, wnn.umap, harmony

save(MEA_InhNeu_withHar, file = "output_InhNeu/01_SubsetInhibitory_AllenPrediction.RData")

load("/scratch/subramanian/Haley_MEA/output_InhNeu/01_SubsetInhibitory_AllenPrediction.RData")

setwd("/Users/SuganyaSubramanian/Haley_MEA_Signac_P6")
load("output_InhNeu/01_SubsetInhibitory_AllenPrediction.RData")

library(MatrixGenerics)
tmp <- as.data.frame.matrix(table(MEA_InhNeu_withHar$SCT_snn_res.0.5,MEA_InhNeu_withHar$predicted.inhid))
tmp$Max_score = rowMaxs(as.matrix(tmp[,c(1,2)]))
tmp$cell.type<-colnames(tmp)[apply(tmp,1,which.max)]
tmp <- tibble::rownames_to_column(tmp, "SCT_snn_res.0.4")
tmp
df <- tmp %>% select(cell.type) %>% rownames_to_column('seurat_clusters')
write.table(df,"output_InhNeu/Celltype_Bymax.txt",sep="\t",row.names=FALSE)

cnt <- table(MEA_InhNeu_withHar$predicted.inhid,MEA_InhNeu_withHar$SCT_snn_res.0.5)
head(cnt)
#cnt[cnt == 0] <- NA
#cnt_NA <- na.omit(cnt)
openxlsx::write.xlsx(cnt,file = "output_InhNeu/Annotation_Count.xlsx", colNames = TRUE, borders = "columns")

smarkers <- c('VIP','ZNF322P1','SSTR1','TOX2',########### Allen prediction Markers
              'PVALB','C8orf4','FAM150B','CNTNAP3P2','SCUBE3',
              'SST','MAFB','MTHFD2P6',
              'LAMP5','DUSP4','SFTA3','NDNF',
              'ADARB2','ADAM33',
              'PAX6','CA4')

jmarkers <- c('PAX6','CDH12','TNFAIP8L3', ###### Notfound BAGE2, CASC6, HYPD, NPM1P10, MIR5548F2
              'LAMP5','NMBR','LCP2','DBP','CA1','PVRL2',
              'SST','CHRNA4','BAGE2',#L1-2
              'ADARB2','MC4R','FAM150B',
              'VIP','SYT6','TSPAN12','CHRNA6','ADAMTSL1','PENK',
              'QPCT','HS3ST3A1','PCDH20','SERPINF1','TYR',
              'CHRM2','CBLN1','CCDC184','GGH','LBH','CASC6','SPAG17','OPRM1',
              'NPY','HYPD','B3GAT2','KLHDC8A','NPM1P10',#SST L3-6
              'GXYLT2','STK32A','CALB1','ADGRG6','FRZB','TH','MIR5548F2',#SST L3-6
              'LHX6','GLP1R',
              'PVALB','LGR5','MEPE','WFDC2','SULF1','SCUBE3')

j2markers <- c('CXCL14','CHRNA7','CNR1','LAMP5','SV2C','SULF1',
               'PDE1A','TOX','SYNPR','ADRA2A','SP8','VIP','NR2F2',
               'MAFB','LHX6','SOX6','SATB1','NES','BEAN1','MEPE')
              


pdf("output_InhNeu/02_DotPlot_check.pdf", width = 10, height = 6)
DotPlot(MEA_InhNeu_withHar, features = smarkers,group.by = "SCT_snn_res.0.5") + rotate_x_text(angle = 45)
dev.off()

pdf("output_InhNeu/02_DotPlot_check_Journal.pdf", width = 18, height = 6)
DotPlot(MEA_InhNeu_withHar, features = jmarkers,group.by = "SCT_snn_res.0.5") + rotate_x_text(angle = 45)
dev.off()

pdf("output_InhNeu/02_DotPlot_check_Journal_2.pdf", width = 18, height = 6)
DotPlot(MEA_InhNeu_withHar, features = j2markers,group.by = "SCT_snn_res.0.5") + rotate_x_text(angle = 45)
dev.off()

pdf("output_InhNeu/02b_UMAP_Res05_Predictionid.pdf", width = 15, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell",label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "predicted.inhid" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.5" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2, p3, ncol = 3)
dev.off()

MEA_InhNeu_withHar@meta.data <- MEA_InhNeu_withHar@meta.data %>%
  mutate(Inh1 = sapply(X = strsplit(MEA_InhNeu_withHar$predicted.inhid, split = " "), FUN = "[", 1)) %>%
  mutate(Inh2 = sapply(X = strsplit(MEA_InhNeu_withHar$predicted.inhid, split = " "), FUN = "[", 2)) %>%
  mutate(Inh3 = sapply(X = strsplit(MEA_InhNeu_withHar$predicted.inhid, split = " "), FUN = "[", 3)) %>%
  mutate(Inh4 = sapply(X = strsplit(MEA_InhNeu_withHar$predicted.inhid, split = " "), FUN = "[", 4)) %>%
  unite(Inh123,Inh1,Inh2,Inh3,sep = "_", remove = FALSE) %>%
  unite(Inh1234,Inh1,Inh2,Inh3,Inh4,sep = "_", remove = FALSE)

pdf("output_InhNeu/03_Prediction_InhLable.pdf", width = 10, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Inh123" ,
        label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.5" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2,ncol = 2)
dev.off()


library(MatrixGenerics)
tmp <- as.data.frame.matrix(table(MEA_InhNeu_withHar$SCT_snn_res.0.5,MEA_InhNeu_withHar$Inh123))
tmp$Max_score = rowMaxs(as.matrix(tmp[,c(1,2)]))
tmp$cell.type<-colnames(tmp)[apply(tmp,1,which.max)]
tmp <- tibble::rownames_to_column(tmp, "SCT_snn_res.0.5")
tmp
df <- tmp %>% dplyr::select(cell.type) %>% rownames_to_column('SCT_snn_res.0.5')
write.table(df,"output_InhNeu/Celltype_Bymax_Inh123.txt",sep="\t",row.names=FALSE)


##### Relabelling again

#Idents(MEA_InhNeu_withHar) <- "seurat_clusters"
Idents(MEA_InhNeu_withHar) <- "SCT_snn_res.0.5"
MEA_InhNeu_withHar@meta.data$seurat_clusters <- MEA_InhNeu_withHar@meta.data$SCT_snn_res.0.5

table(Idents(MEA_InhNeu_withHar))

#Labels <- read.table("output_InhNeu/Try1_LAbel.txt",header=T,sep="\t")
Labels <- read.table("output_InhNeu/Try2_LAbel.txt",header=T,sep="\t")

#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$SubInh)

current.cluster.ids
new.cluster.ids

MEA_InhNeu_withHar@active.ident <- plyr::mapvalues(x = MEA_InhNeu_withHar@active.ident, 
                                                   from = current.cluster.ids, 
                                                   to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
MEA_InhNeu_withHar@meta.data$Sub_InhCell <- MEA_InhNeu_withHar@active.ident

table(MEA_InhNeu_withHar@meta.data$Sub_InhCell)
Inh_L1-3_VIP_GGH  Inh_L4-6_PVALB_SULF1    Inh_L1-3_SST_CALB1 
452                   492                   381 
Inh_L1-4_VIP_CHRNA6   Inh_L4-6_SST_B3GAT2    Inh_L1-2_LAMP5_DBP 
342                   543                   283 
Oligo_MBP           Micro_DOCK4  Inh_L1_ADARB2_ADAM33 
272                   233                   223 
Astro_AQP4    Inh_L2-6_LAMP5_CA1 Inh_L2-5_VIP_SERPINF1 
182                   168                   157 
ExcNeu_SATB2     Inh_L1_LAMP5_NMBR   Inh_L5-6_PVALB_NDNF 
149                    72                    72 
Inh_L1_PAX6_CA4     Inh_L1-3_VIP_SYT6 Inh_L1-6_PVALB_SCUBE3 
70                    58                    58 
Inh_L1-4_LAMP5_DUSP4 
45


pdf("output_InhNeu/08_VlnPlot_nUMI.pdf", width = 20, height = 15)
VlnPlot(MEA_InhNeu_withHar, 
        features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment", "nucleosome_signal","NF"),
        ncol = 3, group.by = "Sub_InhCell",
        pt.size = 0) + theme(legend.position="none")
dev.off()

save(MEA_InhNeu_withHar, file = "output_InhNeu/02_Relabelled_InhNeurons_Res05.RData")

pdf("output_InhNeu/09_FinalLabelled.pdf", width = 8, height = 6)
DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Sub_InhCell" ,
        label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
dev.off()

pdf("output_InhNeu/09a_FinalLabelled_clusters.pdf", width = 12, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Sub_InhCell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "SCT_snn_res.0.5" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
plot_grid(p1,p2, ncol = 2)
dev.off()

pdf("output_InhNeu/09b_UMAP_Labelled_Splitby_Condition.pdf", width = 16, height = 6)
DimPlot(object = MEA_InhNeu_withHar, reduction = "umap", label = TRUE,
        pt.size = 0.5,split.by = "Condition") + theme(legend.position="none") + xlim(-15,15)
dev.off()

order <- c('Inh_L1_ADARB2_ADAM33','Inh_L1_PAX6_CA4',
           'Inh_L1-3_VIP_GGH','Inh_L1-3_VIP_SYT6','Inh_L1-4_VIP_CHRNA6','Inh_L2-5_VIP_SERPINF1',
           'Inh_L4-6_PVALB_SULF1','Inh_L5-6_PVALB_NDNF','Inh_L1-6_PVALB_SCUBE3',
           'Inh_L1-3_SST_CALB1','Inh_L4-6_SST_B3GAT2',
           'Inh_L1_LAMP5_NMBR','Inh_L1-2_LAMP5_DBP','Inh_L1-4_LAMP5_DUSP4','Inh_L2-6_LAMP5_CA1',
           'Oligo_MBP','Micro_DOCK4','Astro_AQP4','ExcNeu_SATB2')
           

Idents(MEA_InhNeu_withHar) <- factor(Idents(MEA_InhNeu_withHar), levels = order) # reorder the factor based on cell type
MEA_InhNeu_withHar$Sub_InhCell <- factor(MEA_InhNeu_withHar$Sub_InhCell, levels = order)

markers <- c('ADARB2','ADAM33','PAX6','CA4',
             'VIP','GGH','SYT6','CHRNA6','SERPINF1',
             'PVALB','SULF1','NDNF','SCUBE3',
             'SST','CALB1','B3GAT2',
             'LAMP5','NMBR','DBP','DUSP4','CA1',
             'MBP','DOCK4','AQP4','SATB2')
pdf("output_InhNeu/10_DotPlot_Final.pdf", width = 10, height = 6)
DotPlot(MEA_InhNeu_withHar, features = markers) + rotate_x_text(angle = 45)
dev.off()



#save(MEA_InhNeu_withHar, file = "output_InhNeu/02_Relabelled_ExcNeurons_Res04.RData")

#load("output_InhNeu/02_Relabelled_ExcNeurons_Res04.RData")
pdf("output_InhNeu/11_FinalLabelled_check.pdf", width = 12, height = 6)
p1 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Sub_InhCell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p2 <- DimPlot(object = MEA_InhNeu_withHar, reduction = "umap",group.by = "Cell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
plot_grid(p1,p2, ncol = 2)
dev.off()


MEA_Inh <- subset(MEA_InhNeu_withHar, subset = Sub_InhCell %in% c('Oligo_MBP','Micro_DOCK4','Astro_AQP4','ExcNeu_SATB2'), invert = TRUE)
table(MEA_Inh@meta.data$Sub_InhCell)

save(MEA_Inh, file = "output_InhNeu/03_Relabelled_ExcNeurons_Res05_NoGlia.RData")

pdf("output_InhNeu/13_FinalLabelled_noGlia.pdf", width = 8, height = 6)
DimPlot(object = MEA_Inh, reduction = "umap",group.by = "Sub_InhCell" ,
        label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
dev.off()

pdf("output_InhNeu/13_UMAP_Labelled_Splitby_Condition_noGlia.pdf", width = 16, height = 6)
DimPlot(object = MEA_Inh, reduction = "umap", label = TRUE,
        pt.size = 0.5,split.by = "Condition") + theme(legend.position="none") + xlim(-15,15)
dev.off()


































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




