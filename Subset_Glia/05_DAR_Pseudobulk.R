# Load libraries
suppressPackageStartupMessages({
  library(scater)
  library(Seurat)
  library(Signac)
  library(tidyverse)
  library(cowplot)
  library(Matrix.utils)
  library(edgeR)
  library(dplyr)
  library(magrittr)
  library(Matrix)
  library(purrr)
  library(reshape2)
  library(S4Vectors)
  library(tibble)
  library(SingleCellExperiment)
  library(pheatmap)
  library(png)
  library(DESeq2)
  library(RColorBrewer)
  library(ggpubr)
})


load("output_Glia_New/03_Glia_Res03_NewLabel.RData")

MEA_Glia
#An object of class Seurat 
#519134 features across 31574 samples within 6 assays 
#Active assay: RNA (20017 features, 2000 variable features)
#3 layers present: counts, data, scale.data
#5 other assays present: SCT, integrated, ATAC, macpeaks, GeneActivity
#6 dimensional reductions calculated: pca, umap, lsi, umap.atac, wnn.umap, harmony

MEA_Glia@meta.data$SubGliaNew <- gsub("_","-",MEA_Glia$SubGliaNew)
MEA_Glia@meta.data$SubGliaNew <- gsub("/","-",MEA_Glia$SubGliaNew)

DefaultAssay(MEA_Glia) <- 'ATAC'

cnt <- AggregateExpression(MEA_Glia,
                           group.by = c("SubGliaNew","Genotype"),
                           slot = "counts",
                           assays = "ATAC",
                           return.seurat = FALSE)
cts <- cnt$ATAC

#transpose column to row
cts.t <- t(cts)
#view
cts.t[1:2,1:2]

splitrows <- gsub("_.*","",rownames(cts.t))
cts.split <- split.data.frame(cts.t,f = factor(splitrows))
cts.split$`Astro-2`[1:2,1:2]            .

cts.split.mod <- lapply(cts.split,function(x){
  rownames(x) <- gsub('.*_(.*)','\\1',rownames(x))
  t(x)
})
cts.split.mod$`Astro-2`[1:2,1:2]

#Get sample names for each of the cell type clusters
# prep. data.frame for plotting
get_sample_ids <- function(x){
  cts.split.mod[[x]] %>%
    colnames()
}

kids <- purrr::set_names(levels(as.factor(MEA_Glia$SubGliaNew)))
kids

de_samples <- map(1:length(kids), get_sample_ids) %>%
  unlist()
head(de_samples)

#Get cluster IDs for each of the samples
samples_list <- map(1:length(kids), get_sample_ids)
get_cluster_ids <- function(x){
  rep(names(cts.split.mod)[x], 
      each = length(samples_list[[x]]))
}
de_cluster_ids <- map(1:length(kids), get_cluster_ids) %>%
  unlist()
head(de_cluster_ids)

# Named vector of Genotype names
sids <- purrr::set_names(levels(as.factor(MEA_Glia$Genotype)))

# Total number of Genotypes 
ns <- length(sids)
m <- match(sids,MEA_Glia$Genotype)

## Turn named vector into a numeric vector of number of cells per Genotype
n_cells <- as.numeric(table(MEA_Glia$Genotype))

## Create the Genotype level metadata by combining the reordered metadata with the number of cells corresponding to each Genotype.
ei <- data.frame(MEA_Glia@meta.data[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("Condition","Genotype","n_cells",
                'Batch','Age','Sex','Hemis','EpDur','SeqBatch')
ei

#Finally, let’s create a data frame with the cluster IDs and the corresponding sample IDs. We will merge together the condition information.
#Create a data frame with the sample IDs, cluster IDs and condition
gg_df <- data.frame(SubGliaNew = de_cluster_ids,
                    Genotype = de_samples) #"SubGliaNew", "Genotype"
head(gg_df)

gg_df <- left_join(gg_df, ei[, c("Genotype","Condition","n_cells",
                                 'Batch','Age','Sex','Hemis','EpDur','SeqBatch')]) 

metadata <- gg_df %>%
  dplyr::select(SubGliaNew,Genotype,n_cells,Condition,Batch,Age,Sex,Hemis,EpDur,SeqBatch)

head(metadata)

dir.create('output_Glia_New/DAR_Pseudobulk/')

save(cts.split.mod,metadata, file = "output_Glia_New/DAR_Pseudobulk/PseudobulkMatrix.RData")

# Generate vector of cluster IDs
clusters <- unique(metadata$SubGliaNew)
clusters

pd <- list()
count <- list()
count_filt <- list()
filter <- list()
dds <- list()
FullTabStats <- list()
DGETabStats <- list()

clusters <- unique(metadata$SubGliaNew)
clusters

for (i in 1:length(clusters)){
  pd[[i]] <- metadata %>% 
    filter(SubGliaNew == clusters[[i]]) %>%
    select(-SubGliaNew,-SeqBatch) %>% 
    relocate(n_cells, .after = EpDur) %>% 
    column_to_rownames("Genotype")
  
  
  count[[i]] <- cts.split.mod[[i]]
  #filter[[i]] <- apply(count[[i]], 1, function(x) (all(x[1:6] > 0 ) | all(x[7:12] > 0)))
  #count_filt[[i]] <- count[[i]][filter[[i]],]
  
  dds[[i]] <- DESeqDataSetFromMatrix(countData = round(count[[i]]),colData = pd[[i]],design = ~ Condition + Batch + Age + Sex + Hemis + EpDur + n_cells)
  dds[[i]] <- estimateSizeFactors(dds[[i]])
  dds[[i]] <- dds[[i]][rowSums(counts(dds[[i]], normalized=TRUE) > 0 ) >= ncol(count[[i]])-2,]
  dds[[i]] <- DESeq(dds[[i]], full=design(dds[[i]]), test="Wald",fitType = "glmGamPoi")
  
  
  
  FullTabStats[[i]] <- as.data.frame(results(dds[[i]],contrast=c("Condition","Stim","Ctrl"),alpha = 0.05)) %>%
    rownames_to_column("Gene") %>%
    mutate(CellType = clusters[[i]])
  
  DGETabStats[[i]] <- FullTabStats[[i]] %>%
    mutate(Abs = abs(log2FoldChange)) %>%
    filter(padj < 0.05 & Abs > 0.3) %>%
    arrange(desc(Abs))
  
  
  names(FullTabStats)[[i]] <- clusters[[i]]
  names(DGETabStats)[[i]] <- clusters[[i]]
  
}

# Save all the tables into RDS
purrr::iwalk(FullTabStats, ~saveRDS(.x, file = paste0("output_Glia_New/DAR_Pseudobulk/",.y, "_FullStats.RDS")))
purrr::iwalk(DGETabStats, ~saveRDS(.x, file = paste0("output_Glia_New/DAR_Pseudobulk/",.y, "_DGEStats.RDS")))


purrr::iwalk(FullTabStats, ~openxlsx::write.xlsx(.x, file = paste0("output_Glia_New/DAR_Pseudobulk/", .y, "_FullStats.xlsx"), colNames=T,rowNames=F,boarders="columns"))
purrr::iwalk(DGETabStats, ~openxlsx::write.xlsx(.x, file = paste0("output_Glia_New/DAR_Pseudobulk/", .y, "_DGEStats.xlsx"), colNames=T,rowNames=F,boarders="columns"))

openxlsx::write.xlsx(do.call(rbind, FullTabStats), 
                     file = "output_Glia_New/DAR_Pseudobulk/AllCellCombined_FullStats.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(do.call(rbind, DGETabStats), 
                     file = "output_Glia_New/DAR_Pseudobulk/AllCellCombined_DGEStats.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")


tmp1 <- do.call(rbind,FullTabStats) %>%
  group_by(CellType) %>%
  mutate(Abs = abs(log2FoldChange)) %>%
  #dplyr::filter(padj < 0.05 & Abs > 0.2) %>%
  mutate(Threshold = if_else(padj < 0.05 & Abs > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(log2FoldChange > 0.3 & padj < 0.05 ~ "UpReg", 
                               log2FoldChange < -0.3 & padj < 0.05 ~ "DownReg")) %>%
  dplyr::filter(Threshold == "TRUE")
table(tmp1$CellType)


color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols <- sample(color,15)

pdf("output_Glia_New/DAR_Pseudobulk/01_Barplot_Number_Signac_DARs.pdf",width = 4,height = 4)
df <- as.data.frame(table(tmp1$CellType))
df$Var1 <- gsub("-","_",df$Var1)

ggplot(df,aes(x=Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity",colour="black") + 
  scale_fill_manual(values=cols) +
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  ggtitle("DARs-Pseudobulk") +
  xlab("Cell") +
  ylab("DARs") + #ylim(0,320) +
  theme_classic() + rotate_x_text(45) +
  theme(legend.position="none")
dev.off()


pdf("output_Glia_New/DAR_Pseudobulk/02_Barplot_Number_Signac_DARs_Direction.pdf",width = 6,height = 4)
df2 <- as.data.frame(table(tmp1$CellType,tmp1$Direction))
df2$Var1 <- gsub("-","_",df2$Var1)
ggplot(df2,aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat='identity', position='dodge',colour="black") +
  scale_fill_manual(values=c('#dda0cd','#092277')) +
  #geom_text(aes(label=Freq)) + #, vjust=-0.3, size=3.5
  ggtitle("DARs-by Pseudobulk") +
  xlab("Cell") +
  ylab("DARs") + #ylim(0,250) +
  theme_classic() + rotate_x_text(45) 
dev.off()

### Closest Genes for DARs
load("output_Glia_New/03_Glia_Res03_NewLabel.RData")

tmp1 <- readxl::read_excel("output_Glia_New/DAR_Pseudobulk/AllCellCombined_DGEStats.xlsx",sheet = "Stats")
head(tmp1)
length(tmp1$baseMean) #18877

DefaultAssay(MEA_Glia) <- 'ATAC'

tmp1 <- tmp1 %>% unite("Gene_Cell", c(Gene,CellType),remove = F)

df <- tmp1 %>%distinct(Gene,.keep_all = TRUE)
df <- df %>% column_to_rownames('Gene') %>% arrange(pvalue)
length(df$baseMean)#18222

#Closest Genes for DARs
darStim <- rownames(df[df$log2FoldChange > 0.3, ])
length(darStim)#9517
darCtrl <- rownames(df[df$log2FoldChange < -0.3, ])
length(darCtrl)#8705

closest_genes_darStim <- ClosestFeature(MEA_Glia, regions = darStim)
head(closest_genes_darStim)

closest_genes_darCtrl <- ClosestFeature(MEA_Glia, regions = darCtrl)
head(closest_genes_darCtrl)


write.table(closest_genes_darStim,"output_Glia_New/DAR_Pseudobulk/ClosestGene_Stim.txt",
            sep="\t",quote=F,row.names=F)

write.table(closest_genes_darCtrl,"output_Glia_New/DAR_Pseudobulk/ClosestGene_Ctrl.txt",
            sep="\t",quote=F,row.names=F)
