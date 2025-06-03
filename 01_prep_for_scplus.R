#prep signac object for scenicplus

library(Seurat)
library(Signac)
library(Matrix)
library(stringr)

seurat_rdata_path = "path"
load(seurat_rdata_path)

DefaultAssay(seurat_obj) <- "ATAC"

#export fragment matrix
frag_mtx <- LayerData(seurat_obj, assay = "ATAC", layer = "counts")
writeMM(frag_mtx, file = "atac_frag_mtx.mtx")

#export fragment names
frag_names <- rownames(seurat_obj)
write.table(frag_names, file = "atac_frag_names.txt",
		quote=F, row.names=F, col.names=F)

#export cell names
cell_names <- colnames(seurat_obj)
write.table(cell_names, file = "cell_names.txt",
		quote=F, row.names=F, col.names=F)

#export cell metadata
metadata <- seurat_obj@meta.data
write.csv(metadata, file = "cell_metadata.csv")


############ if need to update fragment file paths
test_path = seurat_obj@assays$ATAC@fragments[[1]]@path
if (str_detect(test_path, "/zfs/musc3")) {
	for (n in 13:length(seurat_obj@assays$ATAC@fragments)) {
	print(seurat_obj@assays$ATAC@fragments[[n]]@path)
	new_path = stringr::str_replace(seurat_obj@assays$ATAC@fragments[[n]]@path,
		pattern = "/zfs/musc3",
		replacement = "/project/stefanoberto/musc")
	print(new_path)
	seurat_obj@assays$ATAC@fragments[[n]] <- UpdatePath(seurat_obj@assays$ATAC@fragments[[n]],
		new.path = new_path,
		verbose = T)
	}
}

###save frag paths
fp = list()
for (n in 1:length(seurat_obj@assays$ATAC@fragments)) {
	fp <- append(fp, seurat_obj@assays$ATAC@fragments[[n]]@path)
	}
write.table(unlist(fp), "frag_paths.txt", sep = "\t", quote = F, row.names=F, col.names=F)
