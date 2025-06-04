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
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  library(Libra)
})

load("output_InhNeu/02_SeuObj_InhNeurons_NoGlia_ReclusteredwithHarmony.RData")
MEA_InhNeu
#An object of class Seurat 
#427264 features across 3416 samples within 6 assays 
#Active assay: macpeaks (116543 features, 29214 variable features)
#2 layers present: counts, data
#5 other assays present: RNA, SCT, integrated, ATAC, GeneActivity
#6 dimensional reductions calculated: pca, umap, lsi, umap.atac, wnn.umap, harmony

table(MEA_InhNeu$SCT_snn_res.0.5)

MEA_Inh <- MEA_InhNeu
MEA_Inh$seurat_clusters <- MEA_Inh$SCT_snn_res.0.5
Idents(MEA_Inh) <- 'seurat_clusters'
table(MEA_Inh$seurat_clusters)

Labels <- read.table("output_InhNeu/SubCellLabels.txt",header=T,sep="\t")
#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$SubInh)

current.cluster.ids
new.cluster.ids

MEA_Inh@active.ident <- plyr::mapvalues(x = MEA_Inh@active.ident, 
                                        from = current.cluster.ids, 
                                        to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
MEA_Inh@meta.data$SubInhNew <- MEA_Inh@active.ident

table(MEA_Inh@meta.data$SubInhNew)

dir.create("output_InhNeu")

DefaultAssay(MEA_Inh) <- 'RNA'


pdf("output_InhNeu/14_FinalLabelled_NewLabel.pdf", width = 8, height = 6)
DimPlot(object = MEA_Inh, reduction = "umap",group.by = "SubInhNew" ,
        label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
dev.off()

pdf("output_InhNeu/14_UMAP_Labelled_Splitby_Condition_NewLabel.pdf", width = 16, height = 6)
DimPlot(object = MEA_Inh, reduction = "umap", label = TRUE,
        pt.size = 0.5,split.by = "Condition") + theme(legend.position="none") + xlim(-15,15)
dev.off()

order <- c('Pax6_1','Sncg_4',
           'Vip_2','Vip_7','Vip_9','Vip_13','Vip_18',
           'Sst_1','Sst_11','Sst_12','Sst_24',
           'Pvalb_3','Pvalb_4','Pvalb_5',
           'Lamp5_2','Lamp5_4','Lamp5_Lhx6_1')

Idents(MEA_Inh) <- factor(Idents(MEA_Inh), levels = order) # reorder the factor based on cell type
MEA_Inh$SubInhNew <- factor(MEA_Inh$SubInhNew, levels = order)


markers <- c('PAX6','SNCG','VIP','SST','PVALB','LAMP5')

markers2 <- c('ADARB2','ADAM33','PAX6','CA4',
             'VIP','GGH','SYT6','CHRNA6','SERPINF1',
             'PVALB','SULF1','NDNF','SCUBE3',
             'SST','CALB1','B3GAT2',
             'LAMP5','NMBR','DBP','DUSP4','CA1')

pdf("output_InhNeu/15_DotPlot_Final_NewLabel.pdf", width = 6, height = 6)
DotPlot(MEA_Inh, features = markers) + rotate_x_text(angle = 45)
dev.off()

pdf("output_InhNeu/15_DotPlot2_Final_NewLabel.pdf", width = 8, height = 6)
DotPlot(MEA_Inh, features = markers2) + rotate_x_text(angle = 45)
dev.off()

pdf("output_InhNeu/15_VlnPlot_Final_NewLabel.pdf", width = 15, height = 15)
VlnPlot(MEA_Inh, features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment", "nucleosome_signal","NF"),
        ncol = 3, group.by = "SubInhNew",
        pt.size = 0) + theme(legend.position="none")
dev.off()

save(MEA_Inh,file = 'output_InhNeu/03_InhNeurons_Res05_NewLabel.RData')

MEA_Inh

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
  dplyr::arrange(SubInhNew) %>%
  dplyr::group_by(SubInhNew,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(SubInhNew) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(SubInhNew) %>%
  ggbarplot("SubInhNew", "percent",
            fill = "Genotype", color = "Genotype", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_InhNeu/16_SubInhNewProportion_By_Genotype_FinalNoglia.pdf", plot = prop_Cell, width = 18, height = 15, units = "in", dpi = 150)

all_markers_clustID <- presto::wilcoxauc(MEA_Inh, 'SubInhNew', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_InhNeu/SubInhNew-Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_InhNeu/SubInhNew-PrestoTop10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")