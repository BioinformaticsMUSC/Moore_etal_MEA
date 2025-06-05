suppressPackageStartupMessages({
  library(WGCNA)
  library(cluster)
  library(ggplot2)
  library(reshape2)
  library(RColorBrewer)
  library(dplyr)
  library(magrittr)
  library(tidyverse)
  library(ggpubr)
  library(Seurat)
})

################################################################## 
#Prepare the data for enrichment analysis from presto markers file
##################################################################

all_markers_clustID <- presto::wilcoxauc(MEA_Exc, 'seurat_clusters', assay = 'data')

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_ExcNeu/Presto_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

markers <- readxl::read_excel("output_ExcNeu/Presto_padjLT05_logfcGT0.xlsx")

markers_sign <- markers %>%
  dplyr::filter(padj < 0.05, pct_in > 20, pct_out < 20) %>% 
  dplyr::rename(Gene = feature, Cluster = group)

write.table(markers_sign,"output_ExcNeu/Markers_Significant_Res06.txt",sep="\t",quote=F,row.names=F)

df <- markers_sign %>% dplyr::select(Gene,Cluster)
write.table(df,"output_ExcNeu/Markers_Significant_DGE_Res06.txt",sep="\t",quote=F,row.names=F)

################################################################## 
#Prepare the data for enrichment analysis from Jorstad reference
##################################################################
#trying with both ACC and DLPFC Jorstad markers

jormarkers <- readxl::read_excel("output_ExcNeu/Jorstad_Exc_ACC_DLPFC.xlsx") 

Jordf <- jormarkers %>% dplyr::select(gene,cluster)
Jordf <- Jordf %>% dplyr::rename(Gene = gene, Cluster = cluster) 
head(Jordf)
write.table(Jordf,"output_ExcNeu/FE_Jorstad_Exc_ACC_DLPFC_Markers_Significant.txt",sep="\t",quote=F,row.names=F)

jorlist <- split(Jordf,f = Jordf$Cluster)
save(jorlist, file ="output_ExcNeu/PFC_Exc_ACC_DLPFC_GeneSets_FisherEnrichment.RData")

################################################################## 
#Prepare the data for enrichment analysis from presto markers file
##################################################################
allowWGCNAThreads(nThreads = 8) # in R studio

wgcna = list.files(path = "output_ExcNeu",pattern = 'Markers_Significant_DGE_Res04.txt',full.names = TRUE)
tab=read.table(wgcna,sep="\t",header=T)
colnames(tab)=c("Gene","DEFINITION")
tab$DEFINITION <- as.factor(tab$DEFINITION)
Genes=as.data.frame(table(tab$DEFINITION))

# Loop to make the overlap
# The loop goes along the length of the jorlist lists
load("output_ExcNeu/PFC_Exc_ACC_DLPFC_GeneSets_FisherEnrichment.RData")#trying with ACC and PFC
ln=length(jorlist)
cl=length(Genes$Var1)
TEMP=list()
INT=list()
for (i in 1:ln){
  TEMP[[i]]=tab[tab$Gene %in% jorlist[[i]]$Gene,]
  INT[[i]]=as.data.frame(table(TEMP[[i]]$DEFINITION))
}
names(INT)=names(jorlist)
names(TEMP)=names(jorlist)

# Create the matrix for each GeneSet
NROWS <- sapply(jorlist,nrow)

#
#               GeneSet
#              +       -
#          +-------+-------+
#        + |   a   |   b   |
#  Module  +-------+-------+
#        - |   c   |   d   |
#          +-------+-------+

for (i in 1:length(INT)){
  INT[[i]]$b <- NROWS[[i]]-INT[[i]]$Freq
  INT[[i]]$c <- Genes$Freq-INT[[i]]$Freq 
  INT[[i]]$d <- 19776-Genes$Freq-nrow(jorlist[[i]]) #19776 #15585
}

# sum(Genes$Freq)
RunFisher <- function(row, alt = 'greater', cnf = 0.85) {
  f <- fisher.test(matrix(row, nrow = 2), alternative = alt, conf.level = cnf)
  return(c(row,
           P_val = f$p.value,
           LogP = -log10(f$p.value), 
           OR = f$estimate[[1]],
           OR_Low = f$conf.int[1],
           OR_Up = f$conf.int[2]))
}

# run
FisherMat=list()
for (i in 1:length(INT)){
  FisherMat[[i]] <- t(apply(INT[[i]][,2:5], 1, RunFisher))
  rownames(FisherMat[[i]]) <- INT[[i]]$Var1
  FisherMat[[i]] <- FisherMat[[i]][rownames(FisherMat[[i]]) != "grey",]
}
names(FisherMat)<-names(INT)

# Create matrices of Pval
tmp<-list()
FisherP<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,5]))
  FisherP <- do.call(cbind,tmp)
}
rownames(FisherP) <- rowNames
colnames(FisherP) <- colNames

# Create matrices of OR
tmp<-list()
FisherOR<-matrix()
rowNames <- rownames(FisherMat[[1]])
colNames <- names(FisherMat)
for (i in 1:length(INT))
{
  tmp[[i]] <- cbind(as.numeric(FisherMat[[i]][,7]))
  FisherOR <- do.call(cbind,tmp)
}
rownames(FisherOR) <- rowNames
colnames(FisherOR) <- colNames

# Labels by OR
tmp <- as.data.frame(FisherOR)
tmp$Rows <- rownames(tmp)
tmp <- melt(tmp)

tmp2 <- tmp %>% 
  group_by(Rows) %>%
  dplyr::filter(value == max(value)) %>%
  dplyr::arrange(Rows,variable,value) %>%
  as.data.frame()

write.table(tmp2,"output_ExcNeu/Jorstad_Labels_Clusters_ACC&PFC.txt",sep="\t",quote=F)

#library(reshape2)
p <- FisherP %>% 
  as.data.frame() %>%
  rownames_to_column("Module") %>%
  melt() %>%
  group_by(Module) %>%
  mutate(FDR = p.adjust(value,method="BH")) %>%
  mutate(log = -log10(FDR)) %>%
  dplyr::rename(Pval = value) %>%
  as.data.frame()

OR <- FisherOR %>%
  as.data.frame() %>%
  rownames_to_column("Module") %>%
  melt()%>%
  dplyr::rename(OR = value)

p <- Reduce(dplyr::left_join, list(p, OR))
p$OR[!is.finite(p$OR)] <- max(p$OR[is.finite(p$OR)])
p$log[p$log < 1.3] <- NA
p$OR<- ifelse(is.na(p$log), p$log, p$OR)

#p$variable <- factor(p$variable,levels=c("Micro_L1-3","Endo_L2-6","Astro_L1-2", "Astro_L1-6","OPC_L1-6","Oligo_L1-6","Exc_L2","Exc_L2-3","Exc_L2-4","Exc_L3-4","Exc_L3-5","Exc_L4-5","Exc_L4-6","Exc_L5-6",
#                                          "Exc_L6","Exc_L1","Exc_L1-2","Exc_L1-3","Exc_L1-4","Exc_L2-3","Exc_L2-4", "Exc_L2-5","Exc_L2-6","Exc_L3-5","Exc_L3-6","Exc_L4-5","Exc_L4-6","Exc_L5-6"))

pdf("output_ExcNeu/04_Jorstad_Enrichment_ACC&PFC.pdf",width=18,height=8)
ggscatter(p, 
          x = "variable",
          y = "Module",
          size="OR",
          color="log",
          alpha = 0.8,
          xlab = "",
          ylab = "",) +
  theme_minimal() +
  rotate_x_text(angle = 45)+
  #coord_flip()+
  scale_size(range = c(0.5, 10))+ 
  gradient_color(c("red","darkred"))
dev.off()
