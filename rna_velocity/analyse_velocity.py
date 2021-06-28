import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import re
import os



# Parameters ----------------------------------------------------------------

scv.settings.n_jobs = 4
scv.settings.verbosity = 3

# folders when running scvelo locally
# LOOM_PATH = "../data_raw/spliced_unspliced_counts"
# METADATA_FILE = "../data_generated/metadata.csv"
# SAMPLEDATA_FILE = "../metadata/sample_groups.csv"
# NB_SUBCLUSTER_FILE = "../metadata/nb_subclusters.csv"
# NB_SUBUMAP_FILE = "../metadata/nb_subumap.csv"
# OUT_FILE_ALL = "../data_generated/rna_velocity_all.h5ad"
# OUT_FILE_TUMOR = "../data_generated/rna_velocity_tumor.h5ad"

# folders when running scvelo on the cluster
LOOM_PATH = "/media/AGFORTELNY/data/spliced_unspliced_counts"
METADATA_FILE = "/media/AGFORTELNY/analysis/wolfgang/data_generated/metadata.csv"
SAMPLEDATA_FILE = "/media/AGFORTELNY/analysis/wolfgang/metadata/sample_groups.csv"
NB_SUBCLUSTER_FILE = "/media/AGFORTELNY/analysis/wolfgang/metadata/nb_subclusters.csv"
NB_SUBUMAP_FILE = "/media/AGFORTELNY/analysis/wolfgang/metadata/nb_subumap.csv"
OUT_FILE_ALL = "/media/AGFORTELNY/analysis/wolfgang/data_generated/rna_velocity_all.h5ad"
OUT_FILE_TUMOR = "/media/AGFORTELNY/analysis/wolfgang/data_generated/rna_velocity_tumor.h5ad"



# Load data, filter for used barcodes ---------------------------------------

nb_metadata = pd.read_csv(METADATA_FILE).set_index("cell")
sample_metadata = pd.read_csv(SAMPLEDATA_FILE, comment="#").set_index("bsf_id")

samples = []

for sample_path in os.listdir(LOOM_PATH):
    print("Reading", sample_path)
    bsf_id = re.findall("output_(.*)", sample_path)[0]
    sample_name = sample_metadata.loc[bsf_id, "sample"]
    bsf_order = sample_metadata.loc[bsf_id, "bsf_order"]

    sample = scv.read_loom(
      os.path.join(LOOM_PATH, sample_path, "possorted_genome_bam.loom")
    )
    sample.var_names_make_unique()
    sample.obs.index = [
        "{}-{}".format(re.findall(":([ACGT]+)", s)[0], str(bsf_order))
        for s in sample.obs.index
    ]
    sample = sample[np.isin(sample.obs.index, nb_metadata.index)]

    samples.append(sample)


nb_all = anndata.concat(samples)



# Add UMAP coordinates and clusters -----------------------------------------

## Overall dataset ----

nb_metadata = nb_metadata.loc[nb_all.obs.index]
nb_metadata["cellont_cluster"] = (nb_metadata["cellont_cluster"]
                                  .astype("category"))

cluster_colors = (nb_metadata["cellont_abbr"]
                  .replace(dict(B="#33a02c", T="#1f78b4", NK="#a6cee3",
                                M="#ff7f00", NB="#e31a1c", E="#b15928",
                                SC="#6a3d9a", pDC="#fdbf6f", other="black"))
                  .values)

nb_all.obs = nb_metadata
nb_all.obsm["X_umap"] = (nb_metadata[["umap_1_monocle", "umap_2_monocle"]]
                         .values)
nb_all.uns["cluster_colors"] = cluster_colors


## Tumor cells ----

nb_tumor = nb_all[nb_all.obs["cellont_cluster"] == "NB (8)"]

nb_tumor.obsm["X_umap"] = (
    pd.read_csv(NB_SUBUMAP_FILE)
    .set_index("cell")
    .loc[nb_tumor.obs.index]
    [["umap_1_nb", "umap_2_nb"]]
    .values
)

nb_tumor.obs["tumor_subcluster"] = (
    pd.read_csv(NB_SUBCLUSTER_FILE)
    .set_index("cell")
    .loc[nb_tumor.obs.index]
    ["tumor_subcluster"]
    .astype("str")

)

nb_tumor.uns["cluster_colors"] = (
    nb_tumor.obs["tumor_subcluster"]
    .replace({"1": "#FF0000", "2": "#80FF00", "3": "#00FFFF", "4": "#8000FF"})
    .values
)



# Run scvelo ----------------------------------------------------------------

## Overall dataset ----

scv.pp.filter_and_normalize(nb_all, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(nb_all, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(nb_all)
scv.tl.velocity(nb_all, mode="dynamical")
scv.tl.velocity_graph(nb_all)

nb_all.write(OUT_FILE_ALL)


## Tumor cells ----

scv.pp.filter_and_normalize(nb_tumor, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(nb_tumor, n_pcs=30, n_neighbors=30)
scv.tl.recover_dynamics(nb_tumor)
scv.tl.velocity(nb_tumor, mode="dynamical")
scv.tl.velocity_graph(nb_tumor)

nb_tumor.write(OUT_FILE_TUMOR)
