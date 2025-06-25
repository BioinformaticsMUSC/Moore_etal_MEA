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

load("output_Glia/02_SeuObj_Glia_NoUnDef_ReclusteredwithHarmony.RData")

MEA_Glia <- MEA_Glia_withHar_NewLabel

Idents(MEA_Glia) <- 'seurat_clusters'

Labels <- read.table("output_Glia/SubCellLabels.txt",header=T,sep="\t")
#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$SubExc)

current.cluster.ids
new.cluster.ids

MEA_Glia@active.ident <- plyr::mapvalues(x = MEA_Glia@active.ident, 
                                        from = current.cluster.ids, 
                                        to = new.cluster.ids)

### Add New column - Cell with corresponding Labels
MEA_Glia@meta.data$SubGliaNew <- MEA_Glia@active.ident

table(MEA_Glia@meta.data$SubGliaNew)
Oligo_1 Micro/PVM_1a      Astro_4        OPC_1      Astro_2 
14108         4476         3593         4164         1988 
Oligo_2 Micro/PVM_1b       Endo_1      VLMC_1a      VLMC_1b 
983          950          721          224          192 
Micro/PVM_1c       TCells 
97           78

dir.create("output_Glia_New")

DefaultAssay(MEA_Glia) <- 'RNA'
pdf("output_Glia_New/14a_FinalLabelled_clusters_NewLabel.pdf", width = 18, height = 6)
p1 <- DimPlot(object = MEA_Glia, reduction = "umap",group.by = "Sub_GliaCell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p2 <- DimPlot(object = MEA_Glia, reduction = "umap",group.by = "SCT_snn_res.0.3" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p3 <- DimPlot(object = MEA_Glia, reduction = "umap",group.by = "SubGliaNew" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
plot_grid(p1,p2,p3, ncol = 3)
dev.off()

pdf("output_Glia_New/14_FinalLabelled_NewLabel.pdf", width = 8, height = 6)
DimPlot(object = MEA_Glia, reduction = "umap",group.by = "SubGliaNew" ,
        label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
dev.off()

pdf("output_Glia_New/14_UMAP_Labelled_Splitby_Condition_NewLabel.pdf", width = 16, height = 6)
DimPlot(object = MEA_Glia, reduction = "umap", label = TRUE,
        pt.size = 0.5,split.by = "Condition") + theme(legend.position="none") + xlim(-15,15)
dev.off()
                      
order <- c('Astro_2','Astro_4','Micro/PVM_1a','Micro/PVM_1b','Micro/PVM_1c',
           'Oligo_1','Oligo_2','OPC_1','Endo_1','VLMC_1a','VLMC_1b','TCells')

Idents(MEA_Glia) <- factor(Idents(MEA_Glia), levels = order) # reorder the factor based on cell type
MEA_Glia$SubGliaNew <- factor(MEA_Glia$SubGliaNew, levels = order)

DefaultAssay(MEA_Glia) <- 'RNA'

markers <- c('GFAP','GJA1','FGFR3','ETNPPL',
             'MEF2C','P2RY12','TYROBP','CD74','C1QC',
             'MBP','MOBP','OPALIN',
             'MYT1','PCDH15',
             'FLT1','CLDN5',
             'SLC7A11','CPED1','DCN',
             'CD2','CD3D','PTPRC')

pdf("output_Glia_New/15_DotPlot_Final_NewLabel.pdf", width = 8, height = 6)
DotPlot(MEA_Glia, features = markers) + rotate_x_text(angle = 45)
dev.off()

pdf("output_Glia_New/15_VlnPlot_Final_NewLabel.pdf", width = 15, height = 15)
VlnPlot(MEA_Glia, features = c("nUMI","nGene","pMito","nCount_ATAC","TSS.enrichment", "nucleosome_signal","NF"),
        ncol = 3, group.by = "SubGliaNew",
        pt.size = 0) + theme(legend.position="none")
dev.off()

save(MEA_Glia,file = 'output_Glia_New/03_Glia_Res03_NewLabel.RData')

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


prop_Cell <-  MEA_Glia@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(SubGliaNew) %>%
  dplyr::group_by(SubGliaNew,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(SubGliaNew) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(SubGliaNew) %>%
  ggbarplot("SubGliaNew", "percent",
            fill = "Genotype", color = "Genotype", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_Glia_New/16_SubGliaNewProportion_By_Genotype.pdf", plot = prop_Cell, width = 18, height = 15, units = "in", dpi = 150)

all_markers_clustID <- presto::wilcoxauc(MEA_Glia, 'SubGliaNew', assay = 'data')
#all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")
all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_Glia_New/SubGliaNew-Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_Glia_New/SubGliaNew-PrestoTop10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")