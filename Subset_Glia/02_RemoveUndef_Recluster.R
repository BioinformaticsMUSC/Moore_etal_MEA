MEA_Glia_withHar_noundef <- subset(MEA_Glia_withHar, 
                                   subset = seurat_clusters %in% c('7','13'), 
                                   invert = TRUE)

MEA_Glia <- processing_seurat_sctransform(MEA_Glia_withHar_noundef, 
                        vars_to_regress = c("nUMI","pMito","pRibo",
                                          "Age","Sex","Hemis","EpDur","SeqBatch"),
                        npcs = 50, 
                        res = c(0.3,0.4,0.5,0.6))

DefaultAssay(MEA_Glia) <- "RNA"

MEA_Glia <- NormalizeData(object = MEA_Glia, 
                                   normalization.method = "LogNormalize", 
                                   scale.factor = 10000)

library(harmony)
MEA_Glia_undef_withHar <- MEA_Glia

MEA_Glia_undef_withHar <- RunHarmony(MEA_Glia_undef_withHar, assay.use="RNA", group.by.vars = "Genotype")
Reductions(MEA_Glia_undef_withHar)

MEA_Glia_undef_withHar <- RunUMAP(MEA_Glia_undef_withHar, reduction = "harmony", dims = 1:30)
MEA_Glia_undef_withHar <- FindNeighbors(MEA_Glia_undef_withHar, reduction = "harmony", dims = 1:30) 
MEA_Glia_undef_withHar <- FindClusters(MEA_Glia_undef_withHar, resolution = c(0.1,0.2,0.3,0.4,0.5,0.6), 
                                          algorithm = 1, n.iter = 1000,
                                          graph.name = "SCT_snn")

Idents(MEA_Glia_undef_withHar) <- "SCT_snn_res.0.3"
MEA_Glia_undef_withHar@meta.data$seurat_clusters <- MEA_Glia_undef_withHar@meta.data$SCT_snn_res.0.3

pdf("output_Glia/13_FinalLabelled_Noundef_reclust.pdf", width = 18, height = 6)
p1 <- DimPlot(object = MEA_Glia_undef_withHar, reduction = "umap",group.by = "Sub_GliaCell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p2 <- DimPlot(object = MEA_Glia_undef_withHar, reduction = "umap",group.by = "SCT_snn_res.0.3" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
p3 <- DimPlot(object = MEA_Glia_withHar_noundef, reduction = "umap",group.by = "Sub_GliaCell" ,label = TRUE, pt.size = 0.5) + theme(legend.position="none") + xlim(-15,15)
plot_grid(p1,p2,p3, ncol = 3)
dev.off()

pdf("output_Glia/13_FinalLabelled_Noundef_reclust_splitby.pdf", width = 15, height = 15)
DimPlot(object = MEA_Glia_undef_withHar, reduction = "umap",split.by = "SCT_snn_res.0.3",group.by = "Sub_GliaCell",label = TRUE, pt.size = 0.5,ncol=4) + theme(legend.position="none") + xlim(-15,15)
dev.off()

save(MEA_Glia_undef_withHar, file = "output_Glia/02_SeuObj_Glia_NoUnDef_ReclusteredwithHarmony.RData")