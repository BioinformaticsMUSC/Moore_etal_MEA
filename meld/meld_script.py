#meld
import os
import sys
import pandas as pd
import numpy as np 
import graphtools as gt
import phate
import magic
import scprep
import meld
import sklearn
import scipy
import seaborn as sns
import diffxpy.api as de
import scanpy as sc

# setting defaults for matplotlib font sizes
import matplotlib.pyplot as plt
plt.rc('font', size=14)

# making sure plots & clusters are reproducible
np.random.seed(42)

adata_file_path = sys.argv[1]
celltype_col = sys.argv[2]
experimental_col = sys.arv[3]
stim_value = sys.argv[4]

#adata = sc.read_h5ad("/project/stefanoberto/musc/Haley_MEA/adata/final/MEA_exc_adata.h5ad")
adata = sc.read_h5ad(adata_file_path)

CELLTYPE = celltype_col
EXP_COL = experimental_col
STIM_VALUE = stim_value
LIKELIHOOD_COL = STIM_VALUE + "_likelihood"

assert CELLTYPE in adata.obs.columns
assert EXP_COL in adata.obs.columns
assert STIM_VALUE in adata.obs[EXP_COL].values
meld_path = "/project/stefanoberto/musc/Haley_MEA/inh/NEW/meld"
if not os.path.exists(meld_path):
	os.mkdir(meld_path)

os.chdir(meld_path)

data = pd.DataFrame(adata.X.todense(), index=adata.obs_names, columns=adata.var_names)
metadata = adata.obs

# Library size normalize
data = scprep.normalize.library_size_normalize(data)
data = scprep.transform.sqrt(data)

#run PHATE to vizualize data
phateop = phate.PHATE()
data_phate = phateop.fit_transform(data)

#Run MELD to calculate relative likelihood of stimulation
meld_op = meld.MELD()
sample_densities = meld_op.fit_transform(data, metadata[EXP_COL])
sample_densities.index = data.index
metadata[LIKELIHOOD_COL] = meld.utils.normalize_densities(sample_densities)[STIM_VALUE]



scprep.plot.scatter2d(data_phate, c=metadata[LIKELIHOOD_COL],
                      ticks=False, figsize=(6,5),
                     vmin=0, vmax=1, cmap=meld.utils.get_meld_cmap())
plt.tight_layout()
plt.savefig(f"{LIKELIHOOD_COL}_umap.pdf")
adata.obs[LIKELIHOOD_COL] = metadata[LIKELIHOOD_COL]
sc.pl.umap(adata, color=[CELLTYPE, LIKELIHOOD_COL], save="joint_umap_cell.pdf", ncols=1)
sc.pl.violin(adata, keys=LIKELIHOOD_COL, groupby=CELLTYPE, save="stim_cell.pdf", rotation=90)

#thresholding
import sklearn.mixture

mixture_model = sklearn.mixture.GaussianMixture(n_components=3)
classes = mixture_model.fit_predict(adata.obs[LIKELIHOOD_COL].values.reshape(-1,1))
#classes = scprep.utils.sort_clusters_by_values(classes, metadata['stimulated_likelihood'])
adata.obs['classes'] = pd.Categorical(classes)
plt.clf()
g = sns.histplot(data=adata.obs, x=LIKELIHOOD_COL, hue="classes")
plt.tight_layout()
plt.savefig("class1_hist.pdf")
plt.clf()

t = adata.obs.groupby("classes")[LIKELIHOOD_COL].median().sort_values(ascending=True)
mid_class = t.index[1]

def fix_classes(cl, stim_likelihood, mid_class):
	if cl == mid_class:
		return "mid"
	else:
		if stim_likelihood > 0.5:
			return "high"
		else:
			return "low"

def manual_classes(cl, low, high):
	if cl < low:
		return "low"
	elif cl > high:
		return "high"
	return "mid"


adata.obs['new_class'] = adata.obs.apply(lambda x: fix_classes(x['classes'], x[LIKELIHOOD_COL], mid_class), axis=1)

adata.obs['manual_class'] = adata.obs[LIKELIHOOD_COL].apply(lambda x: manual_classes(x, 0.4, 0.57))
adata.obs['new_class'] = pd.Categorical(adata.obs['new_class'], categories = ["high", "low",  "mid"])
adata.obs['manual_class'] = pd.Categorical(adata.obs['manual_class'], categories = ["high", "low",  "mid"])
adata.uns['new_class_colors'] = ["tab:red", "tab:blue", "gainsboro"]
adata.uns['manual_class_colors'] = ["tab:red", "tab:blue", "gainsboro"]
sc.pl.umap(adata, color=[CELLTYPE, "new_class"], ncols=1, save="new_class_cell.pdf", size=5)
###save adata metadata
adata.obs.to_csv("metadata_meld.csv")


##scatterplot
# stim_cmap = {k:v for k,v in zip(adata.obs['Genotype'].unique(), ["red" if x.startswith("St") else "blue" for x in adata.obs['Genotype'].unique()])}
condition_cmap = {"Stim": "tab:red", 
					"Ctrl": "tab:blue"}
plt.clf()
fig, ax = plt.subplots(1, figsize=(10,10))
ax.axhline(y=0.5, linestyle="--", color="gray")
ax.set_ylim(0,1)
scprep.plot.jitter(adata.obs[CELLTYPE], adata.obs[], c=adata.obs[EXP_COL], 
                   cmap=condition_cmap, discrete=True, legend=True, plot_means=True, xlabel=False, ylabel='Mean AD likelihood',
                   ax=ax)
plt.xticks(rotation=90)
#fig.suptitle(dataset.replace("_", " "))
fig.tight_layout()
fig.savefig(f"jitter_celltype1.pdf")


plt.clf()
plt.figure(figsize=(20,5))
#adata.obs['group'] = pd.Categorical(adata.obs['group'], categories=["WT", "AD"], ordered=True)
g = sns.boxplot(data=adata.obs, x=CELLTYPE, y='stim_likelihood', hue=EXP_COL,
                whis=(0,100), palette = condition_cmap)  

plt.xticks(rotation=45, ha="right")
plt.legend(loc="upper right", bbox_to_anchor=(0.55, 0.25))
plt.ylabel("Stim likelihood")
plt.tight_layout()
plt.savefig("boxplot_split2.pdf")

###save adata with meld scores
adata.write(adata_file_path)
