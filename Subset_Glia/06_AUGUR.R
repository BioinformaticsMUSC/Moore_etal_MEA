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
load("output_Glia_New/03_Glia_Res03_NewLabel.RData")

dir.create("output_Glia_New/Augur")

augur1 = calculate_auc(MEA_Glia,
                       cell_type_col = "SubGliaNew", 
                       label_col = "Condition",
                       classifier = "rf", #logistic reg rather than random forest (rf)
                       n_threads = 10)

save(augur1,file = "output_Glia_New/Augur/Glia_Augur.RData")

pdf("output_Glia_New/Augur/Glia_UMAP_Augur1.pdf", width = 3, height = 3)
plot_umap(augur1, MEA_Glia, mode = "rank",
          palette = c("white","red"),cell_type_col = "SubGliaNew")
dev.off()

pdf("output_Glia_New/Augur/Glia_UMAP_Augur1_lollipop.pdf", width = 3, height = 3)
plot_lollipop(augur1)
dev.off()

p1 <- plot_umap(augur1, MEA_Glia, mode = "rank",
                palette = c("white","red"),cell_type_col = "SubGliaNew")
p2 <- DimPlot(object = MEA_Glia, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- plot_lollipop(augur1)

pdf("output_Glia_New/Augur/Glia_Augur_Final_Result.pdf", width = 15, height = 5)
plot_grid(p2, p1, p3, ncol =3)
dev.off()
