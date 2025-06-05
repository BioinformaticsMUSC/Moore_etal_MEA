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
load("output_ExcNeu/03_ExcNeurons_Res04_NewLabel.RData")#MEA_Exc

dir.create("output_ExcNeu/Augur")

augur1 = calculate_auc(MEA_Exc,
                       cell_type_col = "SubExcNew", 
                       label_col = "Condition",
                       classifier = "rf", #logistic reg rather than random forest (rf)
                       n_threads = 10)

save(augur1,file = "output_ExcNeu/Augur/Exc_Augur.RData")

pdf("output_ExcNeu/Augur/Exc_UMAP_Augur1.pdf", width = 3, height = 3)
plot_umap(augur1, MEA_Exc, mode = "rank",
          palette = c("white","red"),cell_type_col = "SubExcNew")
dev.off()

pdf("output_ExcNeu/Augur/Exc_UMAP_Augur1_lollipop.pdf", width = 3, height = 3)
plot_lollipop(augur1)
dev.off()

p1 <- plot_umap(augur1, MEA_Exc, mode = "rank",
                palette = c("white","red"),cell_type_col = "SubExcNew")
p2 <- DimPlot(object = MEA_Exc, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- plot_lollipop(augur1)

pdf("output_ExcNeu/Augur/Exc_Augur_Final_Result.pdf", width = 15, height = 5)
plot_grid(p2, p1, p3, ncol =3)
dev.off()

tmp <- augur1$AUC
write.table(tmp,'output_ExcNeu/Augur/AUC_Score.txt',
            sep="\t",quote=F,row.names=F)
