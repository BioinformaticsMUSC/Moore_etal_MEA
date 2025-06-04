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
  library(CellChat)
  library(BiocSingular)
  library(ComplexHeatmap)
  library(Augur)
})

# Cell type prioritization in high-dimensional single-cell data
load("output_InhNeu/03_InhNeurons_Res05_NewLabel.RData")#MEA_Inh

dir.create("output_InhNeu/Augur")

augur1 = calculate_auc(MEA_Inh,
                       cell_type_col = "SubInhNew", 
                       label_col = "Condition",
                       classifier = "rf", #logistic reg rather than random forest (rf)
                       n_threads = 10)

save(augur1,file = "output_InhNeu/Augur/Inh_Augur.RData")

pdf("output_InhNeu/Augur/Inh_UMAP_Augur1.pdf", width = 3, height = 3)
plot_umap(augur1, MEA_Inh, mode = "rank",
          palette = c("white","red"),cell_type_col = "SubInhNew")
dev.off()

pdf("output_InhNeu/Augur/Inh_UMAP_Augur1_lollipop.pdf", width = 3, height = 3)
plot_lollipop(augur1)
dev.off()

p1 <- plot_umap(augur1, MEA_Inh, mode = "rank",
                palette = c("white","red"),cell_type_col = "SubInhNew")
p2 <- DimPlot(object = MEA_Inh, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- plot_lollipop(augur1)

pdf("output_InhNeu/Augur/Inh_Augur_Final_Result.pdf", width = 15, height = 5)
plot_grid(p2, p1, p3, ncol =3)
dev.off()

tmp <- augur1$AUC
write.table(tmp,'output_InhNeu/Augur/AUC_Score.txt',
            sep="\t",quote=F,row.names=F)
