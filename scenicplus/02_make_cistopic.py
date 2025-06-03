#02a_scenicplus_QC_python
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import os
import scanpy as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import anndata
import pycisTopic
from scipy import io, sparse
import pickle
from pycisTopic.cistopic_class import create_cistopic_object

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['pdf.use14corefonts'] = True

adata_dir = "/project/stefanoberto/musc/Haley_IVS/NEW/exc/adata/"
scplus_dir = "/project/stefanoberto/musc/Haley_IVS/NEW/exc/scenicplus/"
os.chdir(scplus_dir)

#get matrix from processed seurat obj

X = io.mmread(os.path.join(adata_dir, "atac_frag_mtx.mtx"))

X = sparse.csr_matrix(X)


#get fragment file paths
with open(os.path.join(adata_dir, "frag_paths.txt"), "r") as infile:
    fp = infile.read().splitlines()

fragment_path = {}
for f in fp:
    dataset = f.split("/")[7]
    fragment_path[dataset] = f

# for row in df.itertuples():
#     #if s in cell_data['dataset'].unique():
#     fragment_file = row.frag_path
#     if os.path.exists(fragment_file):
#         print(f"{fragment_file} exists")
#         fragment_path[row.sample] = fragment_file

#save fragment paths to pkl file
pickle.dump(fragment_path,
            open(os.path.join(scplus_dir, 'fragment_paths.pkl'), 'wb'))

##if frag_path exists, open via pickle
# with open(os.path.join(scplus_dir, 'fragment_paths.pkl'), 'rb') as infile:
#     fragment_path = pickle.load(infile)

#get cell names
with open(os.path.join(adata_dir, "cell_names.txt"), "r") as cell_infile:
    cell_names = [c.strip() for c in cell_infile.readlines()]

#function to make sure barcodes are the same in RNA and ATAC
def change_sample(x):
    bc = x.split("_")[-1]
    if x.startswith("Ctrl"):
        new_x = f"HM_{x[4:7]}_IVS_Control" + "_" + bc
    else:
        new_x = f"HM_{x[4:7]}_IVS_Stim" + "_" + bc
    return new_x

new_cell_names = [change_sample(t) for t in cell_names]

with open(os.path.join(adata_dir, "atac_frag_names.txt"), "r") as frag_infile:
    region_names = [f.strip() for f in frag_infile.readlines()]

#correct region names to format: "chr??:123-456"
def region_fix(x):
    chrom, start, end = x.split("-")
    return f"{chrom}:{start}-{end}"

new_region_names = [region_fix(r) for r in region_names]

#create cistopic object
cistopic_object = create_cistopic_object(
    fragment_matrix=X,
    cell_names=new_cell_names,
    region_names=new_region_names,
    path_to_blacklist="/project/stefanoberto/musc/reference/exclusion_files/hg38_ENCFF356LFX.bed",
    path_to_fragments=fragment_path)

#add metadata to cistopic object, fix barcodes if needed
cell_data = pd.read_csv(os.path.join(adata_dir, "cell_metadata.csv"), index_col = [0])
cell_data = cell_data.reset_index(names = "barcode")
cell_data['new_barcode'] = cell_data['barcode'].apply(change_sample)
cell_data = cell_data.set_index("new_barcode")

cistopic_object.add_cell_data(cell_data)

##save cistopic object
pickle.dump(
    cistopic_object,
    open(os.path.join(scplus_dir, "cistopic_obj.pkl"), "wb")
)


#### MODELS - may take a while - can be run as a batch script
from pycisTopic.lda_models import run_cgs_models
# Run models
models=run_cgs_models(
    cistopic_object,
    n_topics=[10,20,30],
    n_cpu=12,
    n_iter=100,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    save_path="/project/stefanoberto/musc/Haley_IVS/NEW/exc/scenicplus",
)


pickle.dump(
    models,
    open("models.pkl", "wb")
)
