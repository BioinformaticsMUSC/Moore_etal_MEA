library(Seurat)
library(Signac)
library(stringr)
library(ChIPseeker)
library(ggplot2)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

setwd("~/biocm/projects/Haley_IVS/peaks/MEA_cell_class/")
bedfiles = list.files(pattern="bed")
bedfiles

peakAnnolist <- lapply(bedfiles,
                       annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=TRUE)
names(peakAnnolist) <- str_replace_all(str_replace_all(bedfiles, "^bed_", ""), ".bed", "")
plotAnnoBar(peakAnnolist)
ggsave("MEA_CC_anno_bar.pdf", width = 10, height = 5)
saveRDS(peakAnnolist, "MEA_CC_peakannoList.rds")
