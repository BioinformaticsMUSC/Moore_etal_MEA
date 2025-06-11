suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
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
  library(future)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(EnsDb.Hsapiens.v86)
})

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotation) <- 'UCSC'
genome(annotation) <- "hg38"

##################################################
###### Create seurat and Integrate ATAC ###########
##################################################

dir <- "CR_OUT_ALLB1B2/"
f8 <- list.files(dir)
files <- f8[-c(1, 2, 5, 6)] ##removing sample ctrl/stim 165 and 173

bc=list()
metric_list=list()
fragment_list=list()
mat_list=list()
assay_list=list()
obj_list=list()
peaks=list()
gr=list()

for ( i in 1:length(files)){
  peaks[[i]] <- read.table(file=paste0(dir,files[i],"/outs/atac_peaks.bed"),col.names=c("chr", "start", "end"))
  gr[[i]] <- makeGRangesFromDataFrame(peaks[[i]])
}

combined.peaks <- reduce(do.call(c,gr))
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]

for (i in 1:length(files)){
  bc[[i]] <- read.table(gzfile(paste0(dir,files[i],"/outs/filtered_feature_bc_matrix/barcodes.tsv.gz")),header=F)
  
  metric_list[[i]] <- read.csv(paste0(dir,files[i],"/outs/per_barcode_metrics.csv"),header=T,row.names=1)
  metric_list[[i]] <- metric_list[[i]][bc[[i]]$V1,]
}

for (i in 1:length(files)){
  fragment_list[[i]] <- CreateFragmentObject(path=paste0(dir,files[i],"/outs/atac_fragments.tsv.gz"),cells = bc[[i]]$V1)
  mat_list[[i]] <- FeatureMatrix(fragments = fragment_list[[i]],features = combined.peaks,cells = bc[[i]]$V1)
  assay_list[[i]] <- CreateChromatinAssay(mat_list[[i]], fragments = fragment_list[[i]])
  obj_list[[i]] <- CreateSeuratObject(assay_list[[i]], assay = "ATAC",meta.data=metric_list[[i]])
}
save(obj_list,file="IntegrateATAC/Step1_Objlist.RData")

dir <- "CR_OUT_ALLB1B2/"
f8 <- list.files(dir)
files <- f8[-c(1, 2, 5, 6)] ##removing sample ctrl/stim 165 and 173


for (i in 1: length(files)){
  Annotation(obj_list[[i]]) <- annotation
  obj_list[[i]] <- FindTopFeatures(obj_list[[i]], min.cutoff = 10)
  obj_list[[i]] <- RunTFIDF(obj_list[[i]])
  obj_list[[i]] <- RunSVD(obj_list[[i]])
  obj_list[[i]]<- NucleosomeSignal(object = obj_list[[i]])
  obj_list[[i]]<- TSSEnrichment(object = obj_list[[i]],fast=F)
  obj_list[[i]]$pct_reads_in_peaks <- obj_list[[i]]$atac_peak_region_fragments / obj_list[[i]]$atac_fragments * 100
  obj_list[[i]]$blacklist_ratio <- obj_list[[i]]$blacklist_region_fragments / obj_list[[i]]$peak_region_fragments
  obj_list[[i]]@meta.data$dataset <- files[i]
  obj_list[[i]]@meta.data$BC <- rownames(obj_list[[i]]@meta.data)
  obj_list[[i]] <- RenameCells(obj_list[[i]],new.names=paste(obj_list[[i]]@meta.data$dataset,obj_list[[i]]@meta.data$BC,sep="_"))
}

save(obj_list,file="IntegrateATAC/Step2_Objlist.RData")


combined <- merge(x = obj_list[[1]],y = obj_list[2:length(files)])
combined <- FindTopFeatures(combined, min.cutoff = 'q80')
combined <- RunTFIDF(combined)
combined <- RunSVD(combined)
combined <- RunUMAP(combined, reduction = "lsi", dims = 2:30)

save(combined,file="IntegrateATAC/Step3_combined.RData")

integration.anchors <- FindIntegrationAnchors(object.list = obj_list,anchor.features = rownames(combined),reduction = "rlsi",dims = 2:30)
integrated <- IntegrateEmbeddings(anchorset = integration.anchors,reductions = combined[["lsi"]],new.reduction.name = "integrated_lsi",dims.to.integrate = 1:30,k.weights=100)
integrated <- RunUMAP(integrated, reduction = "integrated_lsi", dims = 2:50)
#save(combined,file="IntegrateATAC/Step3_combined.RData")
save(integrated,file="IntegrateATAC/Step4_Integrated.RData")

#To make genotype column
table(integrated$dataset)
integrated@meta.data <- integrated@meta.data %>%
  mutate(tmp1 = sapply(X = strsplit(integrated$dataset, split = "_"), FUN = "[", 2)) %>%
  mutate(tmp2 = sapply(X = strsplit(integrated$dataset, split = "_"), FUN = "[", 4))

integrated$tmp2 <- gsub("Control","Ctrl",integrated$tmp2)

integrated@meta.data <- integrated@meta.data %>%
  unite("Genotype", c(tmp2,tmp1),sep = "",remove = F)

integrated$tmp1 <- NULL
integrated$tmp2 <- NULL

integrated <- RenameCells(integrated,new.names=paste(integrated@meta.data$Genotype,integrated@meta.data$BC,sep="_"))

save(integrated, file="IntegrateATAC/ATAC_01_Integrated.RData")
