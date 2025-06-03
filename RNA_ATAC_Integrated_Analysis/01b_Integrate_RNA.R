suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(tidyverse)
  library(ggplot2)
  library(sctransform)
  library(hdf5r) #install.packages("hdf5r")
  library(ggrastr)
  library(clustree)
  library(cowplot)
  library(ggpubr)
})

dir <- "cellbender_out/"
#files <- list.files(dir) 
f8 <- fs::dir_ls(dir, recurse = FALSE, glob = "*_filtered.h5") #list by pattern
files <- f8[-c(1, 2, 5, 6)] ##removing sample ctrl/stim 165 and 173

h5_seurat=list()

for (i in 1:length(files)){
  #counts <- Read10X_h5(paste0(dir,files[i]))
  print(paste("working on ",files[i]))
  counts <- Read10X_h5(files[i])
  t1 <- gsub("cellbender_out/","",files[i])
  #print(t1)
  t1 <- gsub("_filtered.h5","",t1)
  print(t1)
  h5_seurat[[i]] <- CreateSeuratObject(counts = counts,min.features = 100,assay="RNA")
  h5_seurat[[i]] <- RenameCells(h5_seurat[[i]], add.cell.id=t1)
  h5_seurat[[i]]$Genotype <- t1
  print(head(h5_seurat[[i]]@meta.data))
}

##### Merge all ########
seuObject  <- merge(h5_seurat[[1]], y = h5_seurat[2:length(h5_seurat)], 
                    project = "MEA")

########## Mito column ##########
seuObject[["pMito"]] <- PercentageFeatureSet(seuObject, pattern = "^MT-")

seuObject[["pRibo"]] <- PercentageFeatureSet(seuObject, pattern = "^RP[SL]")

seuObject@meta.data <- seuObject@meta.data %>% dplyr::rename(nUMI = nCount_RNA, nGene = nFeature_RNA)

seuObject@meta.data <- seuObject@meta.data %>%
  mutate(Condition = case_when(grepl("Ctrl", Genotype) ~ "Ctrl",
                               grepl("Stim", Genotype) ~ "Stim"))

##Add Meta data to seurat
meta <- read.table("MEA_Meta_scrna.txt",sep=",",header=T)

x <- merge(seuObject@meta.data, meta, by.x="Genotype",by.y="Sample")

seuObject@meta.data$Age <- x$Age
seuObject@meta.data$Sex <- as.factor(x$Sex)
seuObject@meta.data$Hemis <- as.factor(x$Hemis)
seuObject@meta.data$EpDur <- x$EpDur
seuObject@meta.data$SeqBatch <- as.factor(x$SeqBatch)

dir.create("Integrate_RNA")

save(seuObject,file="Integrate_RNA/01_SeuObj_Unfilt.RData")

pdf("Integrate_RNA/01_Quality_Control_plots_Genotype.pdf", width=15,height=4)
feats <- c("nUMI", "nGene", "pMito","NF")
VlnPlot(seuObject, group.by = "Genotype", raster=FALSE,features = feats, pt.size = 0) + 
  NoLegend()
dev.off()

load("/zfs/musc3/SS/HgProteinCodingGenesSara.rda")

#filter for ptn coding genes
seuObject <- seuObject[rownames(seuObject)%in%ptn_genes,]

### Filter Mito ######

mito_filtered <- seuObject@assays$RNA@counts[-grep("^MT-",rownames(seuObject@assays$RNA@counts)),]

# Initialize the Seurat object with the raw (non-normalized data).
seuObject_final <- CreateSeuratObject(counts = mito_filtered, project = "Haley_IVS")

## Add pMito info from meta data for all cells before filtering
metaAll <- as.data.frame(seuObject@meta.data)
seuObject_final <- AddMetaData(object = seuObject_final, metadata = as.data.frame(seuObject@meta.data))
seuObject_final@meta.data$nCount_RNA <- NULL
seuObject_final@meta.data$nFeature_RNA <- NULL

save(seuObject_final,file="Integrate_RNA/02_SeuObj_mitoFilt.RData")

######### STEP 10 Data Integration ###############
seuObject_split <- SplitObject(seuObject_final, split.by = "Condition")

seuObject_split <- seuObject_split[c("Ctrl","Stim")]

for (i in 1:length(seuObject_split)) {
  seuObject_split[[i]] <- SCTransform(seuObject_split[[i]],
                                      vars.to.regress = c("nUMI","pMito","pRibo",
                                                          "Age","Sex","Hemis","EpDur","SeqBatch"),
                                      verbose = FALSE)
}

integ_features <- SelectIntegrationFeatures(object.list = seuObject_split,
                                            nfeatures = 4000)

seuObject_split <- PrepSCTIntegration(object.list = seuObject_split,
                                      anchor.features = integ_features)

integ_anchors <- FindIntegrationAnchors(object.list = seuObject_split,
                                        normalization.method = "SCT",
                                        anchor.features = integ_features)

seuObject_integrated <- IntegrateData(
  anchorset = integ_anchors,
  new.assay.name = "integrated",
  normalization.method = "SCT",
  dims = 1:30,
  k.weight = 100,
  sd.weight = 1,
  eps = 0.5,
  verbose = TRUE
)

DefaultAssay(seuObject_integrated) <- "integrated"

seuObject_integrated <- RunPCA(object = seuObject_integrated,
                               features=NULL,
                               weight.by.var = TRUE,
                               ndims.print = 1:5,
                               nfeatures.print = 30,
                               npcs = 50,
                               reduction.name = "pca")

seuObject_integrated <- FindNeighbors(object = seuObject_integrated,
                                      reduction = "pca",
                                      dims = 1:50,
                                      nn.eps = 0.5)

seuObject_integrated <- FindClusters(object = seuObject_integrated,
                                     resolution = c(0.3,0.4,0.5,0.6,0.7,0.8),
                                     algorithm = 1,
                                     n.iter = 1000)

pdf("Integrate_RNA/02_Clustree-Data_Integrated.pdf", width = 15, height = 6)
clustree(seuObject_integrated@meta.data, prefix = "integrated_snn_res.", node_colour = "sc3_stability")
dev.off()

# Select the RNA counts slot to be the default assay
DefaultAssay(seuObject_integrated) <- "RNA"

Idents(object = seuObject_integrated) <- "integrated_snn_res.0.5"

seuObject_integrated <- RunUMAP(object = seuObject_integrated,
                                reduction = "pca",
                                dims = 1:30)

seuObject_integrated <- NormalizeData(object = seuObject_integrated,
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

pdf("Integrate_RNA/03_UMAP-Data_Integrated_Genotype.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_integrated, raster=FALSE,reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_integrated, raster=FALSE,reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("Integrate_RNA/04_UMAP-Data_Integrated_splitby_Condition.pdf", width = 10, height = 6)
DimPlot(object = seuObject_integrated, raster=FALSE,reduction = "umap", ncol = 3, label = FALSE, pt.size = 0.5, split.by = "Condition")
dev.off()

save(seuObject_integrated, file = "Integrate_RNA/03_seuObject_Integrated.RData")

##Add Nuclear Fraction to metadata
load("DropletQC/DropletQC_Result_Processed_Merged_withbryan.RData")

rna_names<-colnames(seuObject_integrated) #merged seurat object from raw.h5

nf_names<-rownames(new_df)

intersect<-intersect(nf_names, rna_names)
length(intersect) #82283

subdf <- new_df[rownames(new_df) %in% intersect, ]

seuObject_integrated@meta.data$NF <- subdf$nuclear_fraction
seuObject_integrated@meta.data$NF_cell <- subdf$cell

save(seuObject_integrated,file="Integrate_RNA/03_seuObject_Integrated_withNF.RData")

