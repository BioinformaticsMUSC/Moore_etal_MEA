suppressPackageStartupMessages({
  library(tidyverse)
  library(MetBrewer)
  library(SCpubr)
  library(ggsankey)
  library(ggplot2)
  library(ggrepel)
  library(ggpubr)
  library(data.table)
  library(RColorBrewer)
  library(tidyverse)
  library(preprocessCore)
  library(future.apply)
  library(DESeq2)
  library(pheatmap)
  library(sva)
  library(viridis)
  library(limma)
  library(janitor)
  library(UpSetR)
  #library(GGally)
  library(ComplexHeatmap)
  library(igraph)
  library(here)
  library(Seurat)
  library(SeuratExtend)
  library(variancePartition)
  library(scToppR)
})

load("output_ExcNeu/03_ExcNeurons_Res04_NewLabel.RData")#MEA_Exc

MEA_Exc@meta.data$SubExcNew <- gsub("_","-",MEA_Exc$SubExcNew)
MEA_Exc@meta.data$SubExcNew <- gsub("/","-",MEA_Exc$SubExcNew)

DefaultAssay(MEA_Exc) <- 'ATAC'

cnt <- AggregateExpression(MEA_Exc,
                           group.by = c("SubExcNew","Genotype"),
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
cts.split$`Lamp5-2`[1:2,1:2]

cts.split.mod <- lapply(cts.split,function(x){
  rownames(x) <- gsub('.*_(.*)','\\1',rownames(x))
  t(x)
})
cts.split.mod$`Lamp5-2`[1:2,1:2]

#Get sample names for each of the cell type clusters
# prep. data.frame for plotting
get_sample_ids <- function(x){
  cts.split.mod[[x]] %>%
    colnames()
}

kids <- purrr::set_names(levels(as.factor(MEA_Exc$SubExcNew)))
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

# Named vector of Genotype names
sids <- purrr::set_names(levels(as.factor(MEA_Exc$Genotype)))

# Total number of Genotypes 
ns <- length(sids)
m <- match(sids,MEA_Exc$Genotype)

## Turn named vector into a numeric vector of number of cells per Genotype
n_cells <- as.numeric(table(MEA_Exc$Genotype))

## Create the Genotype level metadata by combining the reordered metadata with the number of cells corresponding to each Genotype.
ei <- data.frame(MEA_Exc@meta.data[m, ], 
                 n_cells, row.names = NULL) %>% 
  dplyr::select("Condition","Genotype","n_cells",
                'Batch','Age','Sex','Hemis','EpDur','SeqBatch')
ei

#Finally, let’s create a data frame with the cluster IDs and the corresponding sample IDs. We will merge together the condition information.
#Create a data frame with the sample IDs, cluster IDs and condition
gg_df <- data.frame(SubExcNew = de_cluster_ids,
                    Genotype = de_samples) #"SubExcNew", "Genotype"
head(gg_df)

gg_df <- left_join(gg_df, ei[, c("Genotype","Condition","n_cells",
                                 'Batch','Age','Sex','Hemis','EpDur','SeqBatch')]) 

metadata <- gg_df %>%
  dplyr::select(SubExcNew,Genotype,n_cells,Condition,Batch,Age,Sex,Hemis,EpDur,SeqBatch)

head(metadata)

dir.create('output_ExcNeu/DAR_Pseudobulk_WALD_Signac')

write.table(metadata,'output_ExcNeu/DAR_Pseudobulk_WALD_Signac/Pseudobulk_Meta.txt',
            sep="\t",quote=F,row.names=T)

save(cts.split.mod,metadata, file = "output_ExcNeu/DAR_Pseudobulk_WALD_Signac/PseudobulkMatrix.RData")

pd <- list()
count <- list()
count_filt <- list()
filter <- list()
dds <- list()
FullTabStats <- list()
DGETabStats <- list()

clusters <- unique(metadata$SubExcNew)
clusters

for (i in 1:length(clusters)) 
{
  pd[[i]] <- metadata %>% 
    filter(SubExcNew == clusters[[i]]) %>%
    select(-SubExcNew,-SeqBatch) %>% 
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
purrr::iwalk(FullTabStats, ~saveRDS(.x, file = paste0("output_ExcNeu/DAR_Pseudobulk_WALD_Signac/",.y, "_FullStats.RDS")))
purrr::iwalk(DGETabStats, ~saveRDS(.x, file = paste0("output_ExcNeu/DAR_Pseudobulk_WALD_Signac/",.y, "_DGEStats.RDS")))


purrr::iwalk(FullTabStats, ~openxlsx::write.xlsx(.x, file = paste0("output_ExcNeu/DAR_Pseudobulk_WALD_Signac/", .y, "_FullStats.xlsx"), colNames=T,rowNames=F,boarders="columns"))
purrr::iwalk(DGETabStats, ~openxlsx::write.xlsx(.x, file = paste0("output_ExcNeu/DAR_Pseudobulk_WALD_Signac/", .y, "_DGEStats.xlsx"), colNames=T,rowNames=F,boarders="columns"))

openxlsx::write.xlsx(do.call(rbind, FullTabStats), 
                     file = "output_ExcNeu/DAR_Pseudobulk_WALD_Signac/AllCellCombined_FullStats.xlsx",
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Stats")

openxlsx::write.xlsx(do.call(rbind, DGETabStats), 
                     file = "output_ExcNeu/DAR_Pseudobulk_WALD_Signac/AllCellCombined_DGEStats.xlsx", 
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



color = grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)]
cols <- sample(color,10)

pdf("output_ExcNeu/DAR_Pseudobulk_WALD_Signac/01_Barplot_Number_Signac_DARs.pdf",width = 4,height = 4)
df <- as.data.frame(table(tmp1$CellType))
df$Var1 <- gsub("-","_",df$Var1)

ggplot(df,aes(x=Var1, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity",colour="black") + 
  scale_fill_manual(values=cols) +
  geom_text(aes(label=Freq), vjust=-0.3, size=3.5) +
  ggtitle("DEGs-Pseudobulk") +
  xlab("Cell") +
  ylab("DEGs") + #ylim(0,320) +
  theme_classic() + rotate_x_text(45) +
  theme(legend.position="none")
dev.off()


pdf("output_ExcNeu/DAR_Pseudobulk_WALD_Signac/02_Barplot_Number_Signac_DARs_Direction.pdf",width = 6,height = 4)
df2 <- as.data.frame(table(tmp1$CellType,tmp1$Direction))
df2$Var1 <- gsub("-","_",df2$Var1)
ggplot(df2,aes(x=Var1, y=Freq, fill=Var2)) +
  geom_bar(stat='identity', position='dodge',colour="black") +
  scale_fill_manual(values=c('#dda0cd','#092277')) +
  #geom_text(aes(label=Freq)) + #, vjust=-0.3, size=3.5
  ggtitle("DEGs-by Pseudobulk") +
  xlab("Cell") +
  ylab("DEGs") + #ylim(0,250) +
  theme_classic() + rotate_x_text(45) 
dev.off()

### Closest Genes for DARs

tmp1 <- readxl::read_excel("output_ExcNeu/DAR_Pseudobulk_WALD_Signac/AllCellCombined_DGEStats.xlsx",sheet = "Stats")
head(tmp1)

DefaultAssay(MEA_Exc) <- 'ATAC'

df <- tmp1 %>% column_to_rownames('Gene') %>% arrange(pvalue)
length(df$baseMean)#163
#Closest Genes for DARs
darStim <- rownames(df[df$log2FoldChange > 0.3, ])
length(darStim)#38
darCtrl <- rownames(df[df$log2FoldChange < -0.3, ])
length(darCtrl)#125

closest_genes_darStim <- ClosestFeature(MEA_Exc, regions = darStim)
head(closest_genes_darStim)

closest_genes_darCtrl <- ClosestFeature(MEA_Exc, regions = darCtrl)
head(closest_genes_darCtrl)

write.table(closest_genes_darStim,"output_ExcNeu/DAR_Pseudobulk_WALD_Signac/ClosestGene_Stim.txt",
            sep="\t",quote=F,row.names=F)

write.table(closest_genes_darCtrl,"output_ExcNeu/DAR_Pseudobulk_WALD_Signac/ClosestGene_Ctrl.txt",
            sep="\t",quote=F,row.names=F)







