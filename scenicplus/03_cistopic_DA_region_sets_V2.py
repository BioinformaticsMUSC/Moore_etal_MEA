import os
import sys
import pickle
import pycisTopic
from pycisTopic.lda_models import evaluate_models
from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['pdf.use14corefonts'] = True

####TAKEN FROM TUTORIAL HERE:
# https://pycistopic.readthedocs.io/en/latest/notebooks/human_cerebellum.html#Model-selection


work_dir = "/project/stefanoberto/musc/Haley_IVS/NEW/exc/scenicplus"
model_path  = "/project/stefanoberto/musc/Haley_IVS/NEW/exc/scenicplus/models.pkl"
os.chdir(work_dir)

CELLTYPE = "SubExcNew"

###Make region sets directories

os.makedirs(os.path.join(work_dir, "region_sets"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "region_sets", "Topics_otsu"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "region_sets", "Topics_top_3k"), exist_ok = True)
os.makedirs(os.path.join(work_dir, "region_sets", "DARs_cell_type"), exist_ok = True)


#### ADD MODELS FROM PREVIOUS

with open(os.path.join(work_dir, "cistopic_obj.pkl"), "rb") as infile:
    cistopic_obj = pickle.load(infile)

with open(model_path,"rb") as infile:
    models = pickle.load(infile)

model = evaluate_models(models)
plt.tight_layout()
plt.savefig("model_eval.pdf")

###need to check the plot to determine the best model
selected_model = evaluate_models(models,
    select_model = 20,
    return_model = True)

cistopic_obj.add_LDA_model(selected_model)

with open(os.path.join(work_dir, "cistopic_obj_model.pkl"), "wb") as outfile:
    pickle.dump(cistopic_obj, outfile)

###VIZ and cluster

from pycisTopic.clust_vis import (
    find_clusters,
    run_umap,
    run_tsne,
    plot_metadata,
    plot_topic,
    cell_topic_heatmap
)


find_clusters(
    cistopic_obj,
    target  = 'cell',
    k = 10,
    res = [0.6, 1.2, 3],
    prefix = 'pycisTopic_',
    scale = True,
)

run_umap(
    cistopic_obj,
    target  = 'cell', scale=True)

plot_topic(
    cistopic_obj,
    reduction_name = 'UMAP',
    target = 'cell',
    num_columns=5
)
plt.savefig("02_pycisTopic_topics.pdf")

plot_metadata(
    cistopic_obj,
    reduction_name='UMAP',
    variables=[CELLTYPE, "Condition"],
    target='cell', num_columns=2,
    text_size=10,
    dot_size=2)
plt.savefig("03_pycisTopic_umap_new.pdf")

####### BINARIZE
from pycisTopic.topic_binarization import binarize_topics
region_bin_topics_top_3k = binarize_topics(
    cistopic_obj, method='ntop', ntop = 3_000,
    plot=True, num_columns=5
)
plt.savefig("04_binarize_topics_top_3k.pdf")

region_bin_topics_otsu = binarize_topics(
    cistopic_obj, method='otsu',
    plot=True, num_columns=5
)
plt.savefig("05_binarize_topics_otsu.pdf")

binarized_cell_topic = binarize_topics(
    cistopic_obj,
    target='cell',
    method='li',
    plot=True,
    num_columns=5, nbins=100)
plt.savefig("06_binarized_cell_topics.pdf")

### TOPIC QC
from pycisTopic.topic_qc import compute_topic_metrics, plot_topic_qc, topic_annotation
import matplotlib.pyplot as plt
from pycisTopic.utils import fig2img

topic_qc_metrics = compute_topic_metrics(cistopic_obj)

fig_dict={}
fig_dict['CoherenceVSAssignments']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Log10_Assignments', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['AssignmentsVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Log10_Assignments', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSCells_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Cells_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSRegions_in_bin']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Regions_in_binarized_topic', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSMarginal_dist']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Marginal_topic_dist', var_color='Gini_index', plot=False, return_fig=True)
fig_dict['CoherenceVSGini_index']=plot_topic_qc(topic_qc_metrics, var_x='Coherence', var_y='Gini_index', var_color='Gini_index', plot=False, return_fig=True)

# Plot topic stats in one figure
fig=plt.figure(figsize=(40, 43))
i = 1
for fig_ in fig_dict.keys():
    plt.subplot(2, 3, i)
    img = fig2img(fig_dict[fig_]) #To convert figures to png to plot together, see .utils.py. This converts the figure to png.
    plt.imshow(img)
    plt.axis('off')
    i += 1


plt.subplots_adjust(wspace=0, hspace=-0.70)
plt.savefig("07_TOPIC_QC.pdf")


###ANNOTATE TOPICS BY CELLTYPE
topic_annot = topic_annotation(
    cistopic_obj,
    annot_var=CELLTYPE,
    binarized_cell_topic=binarized_cell_topic,
    general_topic_thr = 0.2
)
topic_annot.to_csv("annotated_topics.csv")

from pycisTopic.utils import region_names_to_coordinates
for topic in region_bin_topics_otsu:
    region_names_to_coordinates(
        region_bin_topics_otsu[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "region_sets", "Topics_otsu", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )


for topic in region_bin_topics_top_3k:
    region_names_to_coordinates(
        region_bin_topics_top_3k[topic].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "region_sets", "Topics_top_3k", f"{topic}.bed"),
        sep = "\t",
        header = False, index = False
    )



###SAVE CISTOPIC
with open(os.path.join(work_dir, "cistopic_obj_model_preDAR.pkl"), "wb") as outfile:
    pickle.dump(cistopic_obj, outfile)


with open(os.path.join(work_dir, "cistopic_obj_model_preDAR.pkl"), "rb") as infile:
    cistopic_obj = pickle.load(infile)
### FIND DARs

from pycisTopic.diff_features import (
    impute_accessibility,
    normalize_scores,
    find_highly_variable_features,
    find_diff_features
)

imputed_acc_obj = impute_accessibility(
    cistopic_obj,
    selected_cells=None,
    selected_regions=None,
    scale_factor=10**6
)

normalized_imputed_acc_obj = normalize_scores(imputed_acc_obj, scale_factor=10**4)

variable_regions = find_highly_variable_features(
    normalized_imputed_acc_obj,
    min_disp = 0.05,
    min_mean = 0.0125,
    max_mean = 3,
    max_disp = np.inf,
    n_bins=20,
    n_top_features=None,
    plot=True,
    save="HVG_imputed_acc.pdf"
)
len(variable_regions)

markers_dict= find_diff_features(
    cistopic_obj,
    imputed_acc_obj,
    variable=CELLTYPE,
    var_features=variable_regions,
    contrasts=None,
    adjpval_thr=0.05,
    log2fc_thr=np.log2(1.5),
    n_cpu=5,
    _temp_dir='/scratch/bryangranger/tmp',
    split_pattern = '-'
)
from pycisTopic.clust_vis import plot_imputed_features

plot_imputed_features(
    cistopic_obj,
    reduction_name='UMAP',
    imputed_data=imputed_acc_obj,
    features=[markers_dict[x].index.tolist()[0] for x in markers_dict.keys()],
    scale=False,
    num_columns=4,
    save="plot_imputed_features.pdf"
)


#if need to change the celltypes to remove /s
new_markers_dict = {k.replace("/", "_"):v for k,v in markers_dict.items()}

####create region sets from earlier code


from pycisTopic.utils import region_names_to_coordinates
for cell_type in new_markers_dict:
    region_names_to_coordinates(
        new_markers_dict[cell_type].index
    ).sort_values(
        ["Chromosome", "Start", "End"]
    ).to_csv(
        os.path.join(work_dir, "region_sets", "DARs_cell_type", f"{cell_type}.bed"),
        sep = "\t",
        header = False, index = False
    )


with open(os.path.join(work_dir, "cistopic_obj_final.pkl"), "wb") as outfile:
    pickle.dump(cistopic_obj, outfile)

