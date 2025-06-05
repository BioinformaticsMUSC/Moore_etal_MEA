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


load("output_ExcNeu/02_SeuObj_ExcNeurons_NoGlia_ReclusteredwithHarmony.RData")

### Removing low umi clusters before labelling
MEA_Exc <- subset(MEA_ExcNeu_withHar_NewLabel, subset = seurat_clusters %in% c('5','9'),invert = TRUE)

MEA_Exc$seurat_clusters <- MEA_Exc$seurat_clusters %>% droplevels()

Idents(MEA_Exc) <- 'seurat_clusters'

Labels <- read.table("output_ExcNeu/SubCellLabels.txt",header=T,sep="\t")
#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$SubExc)

current.cluster.ids
new.cluster.ids

MEA_Exc@active.ident <- plyr::mapvalues(x = MEA_Exc@active.ident, 
                                        from = current.cluster.ids, 
                                        to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
MEA_Exc@meta.data$SubExcNew <- MEA_Exc@active.ident

pdf("output_ExcNeu/14a_FinalLabelled_clusters_NewLabel.pdf", width = 18, height = 6)
p1 <- DimPlot(object = MEA_Exc, reduction = "umap",group.by = "Sub_ExcCell_Final" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p2 <- DimPlot(object = MEA_Exc, reduction = "umap",group.by = "SCT_snn_res.0.4" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p3 <- DimPlot(object = MEA_Exc, reduction = "umap",group.by = "SubExcNew" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
plot_grid(p1,p2,p3, ncol = 3)
dev.off()

pdf("output_ExcNeu/14_FinalLabelled_NewLabel.pdf", width = 8, height = 6)
DimPlot(object = MEA_Exc, reduction = "umap",group.by = "SubExcNew" ,
        label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
dev.off()

pdf("output_ExcNeu/14_UMAP_Labelled_Splitby_Condition_NewLabel.pdf", width = 16, height = 6)
DimPlot(object = MEA_Exc, reduction = "umap", label = TRUE,
        pt.size = 0.5,split.by = "Condition") + theme(legend.position="none") + xlim(-15,15)
dev.off()

order <- c('L6_IT_1','L6_IT_Car3','L6_CT_1','L6b_3',
           'L5/6_NP_2','L5_IT_1','L5_IT_2',
           'L4_IT_1','L4_IT_2',
           'L2/3_IT_1a','L2/3_IT_1b','L2/3_IT_4')

Idents(MEA_Exc) <- factor(Idents(MEA_Exc), levels = order) # reorder the factor based on cell type
MEA_Exc$SubExcNew <- factor(MEA_Exc$SubExcNew, levels = order)

markers <- c('LAMP5','LTK','RORB','DAPK2','TWIST2','FILIP1L','C1R',
             'FEZF2','MYBPHL','ABO',
             'THEMIS','C1QL3','FGF10','HS3ST4')

pdf("output_ExcNeu/15_DotPlot_Final_NewLabel.pdf", width = 10, height = 6)
DotPlot(MEA_Exc, features = markers) + rotate_x_text(angle = 45)
dev.off()

pdf("output_ExcNeu/15_VlnPlot_Final_NewLabel.pdf", width = 15, height = 15)
VlnPlot(MEA_Exc, features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment", "nucleosome_signal","NF"),
        ncol = 3, group.by = "SubExcNew",
        pt.size = 0) + theme(legend.position="none")
dev.off()

save(MEA_Exc,file = 'output_ExcNeu/03_ExcNeurons_Res04_NewLabel.RData')

color <- c("#92a8d1","#0e3386",
           "#c071d0","#5c2893",
           "#db821d","#890c0a",
           "#a6ea7c","#28935c",
           "#b0b0b2","#4c3c31",
           "#e2fd07","#7f6000",
           "#fddad2","#f6471c",
           "#cdf5ff","#05cefe")


prop_Cell <-  MEA_Exc@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(SubExcNew) %>%
  dplyr::group_by(SubExcNew,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(SubExcNew) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(SubExcNew) %>%
  ggbarplot("SubExcNew", "percent",
            fill = "Genotype", color = "Genotype", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_ExcNeu/16_SubExcNewProportion_By_Genotype_FinalNoglia.pdf", plot = prop_Cell, width = 18, height = 15, units = "in", dpi = 150)

all_markers_clustID <- presto::wilcoxauc(MEA_Exc, 'SubExcNew', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_ExcNeu/SubExcNew-Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_ExcNeu/SubExcNew-PrestoTop10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")