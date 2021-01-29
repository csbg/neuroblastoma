# Analysis of scRNA-seq and scATAC-seq data from bone marrow cells of ultra-high-risk neuroblastoma patients

## Folders

* `data_generated`: output files generated by the scripts in this repository. The following datasets are available:

  - `3_datasets`: Three datasets integrated with SCTransform; small datasets that can be used for testing purposes
  - `all_datasets_v1`: Initial attempt at integrating all datasets with SCTransform; one sample was duplicated (with different read lengths)
  - `all_datasets_v2`: All (13) datasets integrated with SCTransform, but no QC
  - `all_datasets_v3`: All (13) datasets with QC based on feature count and mitochondrial genes
  - `all_datasets_current`: Adds four additional datasets (now 17)

* `data_raw`: Cell Ranger output downloaded from [here](https://biomedical-sequencing.at/projects/BSA_0407_STM_Neuroblastoma_2ba0210fb73d412397728e8a97a3e423). Subfolders ending in `_transcriptome` and `_ATAC` contain scRNA-seq and ATAC-seq data, respectively. `metadata` contains various additional information (as provided by our collaborators). Notably, `sample_groups.csv` contains mappings from samples to neuroblastoma groups.

* `doc`: project documentation

* `literature`: relevant publications

* `plots`: generated plots


## Main workflow

* [integrate_samples.R](integrate_samples.R) - integrate selected scRNA-seq samples; invoke via (sge_job_integrate.sh)[sge_job_integrate.sh]
* [extract_seurat_data.R] - run basic analyses (dimensional reduction, clustering ...) and extract resulting data from a single (large) Seurat object to (smaller) CSV files and Seurat objects that only contain a single assay
* [cell_types.R] - perform cell type classification via singler; invoke via [sge_job_singler.sh]
* [detect_doublets.R] - doublet detection via scds; invoke via [sge_job_doublet.sh]
* [make_subclusters.R] - perform subclustering
* [assemble_metadata.R] - generate one CSV and RDS file with all metadata
* [conserved_markers.R] - find conserved markers; invoke via [sge_job_consmark.sh]
* [analysis_integration.R] - analysis of integrated data
* [plot_markers.R] - plot canonical cell type and NB markers
* [diff_exp_seurat.R] - differential expression via Seurat; invoke via [sge_job_de_seurat.sh]

## Other analyses

* [barcodes_seurat_cellranger.R] - create a table that maps cellranger and Seurat barcodes
* [compare_reads.Rmd] - compare two samples that differ in read length (50 vs 75 bp)
* [download_data.R] - download raw data from the BSF
* [merge_samples.R] - merge selected samples without actually integrating them; invoke via [sge_job_merge.sh]
* [plot_qc_metrics.R] - plot QC metrics for cell filtering (i.e., % of mitochondrial genes etc.)
* [summary_cellranger_metrics.R] - generate summary plots of Cell Ranger QC metrics
* [summary_cellranger_tsne.R] - generate overview of tSNEs calculated by Cell Ranger
