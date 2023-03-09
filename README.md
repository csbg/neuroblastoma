# Dissecting the cellular architecture of neuroblastoma bone marrow metastasis using single-cell transcriptomics and epigenomics unravels the role of monocytes at the metastatic niche

## Folders

(Not all of these folders appear in the git repository.)

- `data_generated`: output files generated by the scripts in this repository
- `data_raw`: raw input data
- `doc`: project documentation
- `literature`: relevant publications
- `metadata`: additional required data
- `misc`: miscellaneous scripts
- `plots`: generated plots
- `renv`: R environment data
- `tables`: exported tables



## Download data

Create a folder `data_raw` that will contain raw data in the following subfolders:

- `adrmed`:
  - `adrenal_medulla_Seurat.RDS`: reference expression data for adrenal medullary cells; download from https://adrenal.kitz-heidelberg.de/developmental_programs_NB_viz/ (Download data -> Download Adrenal medulla data -> Seurat object (RDS))
- `rna_seq`: Download `GSE216155_RAW.tar` from GEO Series GSE216155 and extract all files.
- `atac_seq`: Download all files from GEO Series GSE216175 (`GSE216175_barcodes.tsv.gz`, `GSE216175_barcodes_samples.csv.gz`, `GSE216175_filtered_peak_bc_matrix.h5`, `GSE216175_matrix.mtx.gz`, `GSE216175_peaks.bed.gz`, and `GSE216175_RAW.tar`). Extract all files from the tarball.
- `GSE137804`: download the following files from GEO series GSE137804:
  - `GSE137804_tumor_dataset_annotation.csv.gz`
  - `GSE137804_RAW.tar`, from which the following eleven files must be extracted:
    - `GSM4088774_T10_gene_cell_exprs_table.xls.gz`
    - `GSM4088776_T27_gene_cell_exprs_table.xls.gz`
    - `GSM4088777_T34_gene_cell_exprs_table.xls.gz`
    - `GSM4088780_T69_gene_cell_exprs_table.xls.gz`
    - `GSM4088781_T71_gene_cell_exprs_table.xls.gz`
    - `GSM4088782_T75_gene_cell_exprs_table.xls.gz`
    - `GSM4088783_T92_gene_cell_exprs_table.xls.gz`
    - `GSM4654669_T162_gene_cell_exprs_table.xls.gz`
    - `GSM4654672_T200_gene_cell_exprs_table.xls.gz`
    - `GSM4654673_T214_gene_cell_exprs_table.xls.gz`
    - `GSM4654674_T230_gene_cell_exprs_table.xls.gz`
- `snp_array`: Extract the contents of `snp_array.tgz` provided in Zenodo repository https://doi.org/10.5281/zenodo.7707614

Optionally, obtain intermediary data: Extract the contents of `R_data_generated.tgz` from Zenodo repository https://doi.org/10.5281/zenodo.7707614 to folder `data_generated`.



## Main workflow

Run these R scripts in the given order to generate all files
required by figures and tables.

- [perform_qc.R](perform_qc.R) -
  perform QC filtering, ensure unique cell names
- [integrate_rna.R](integrate_rna.R) -
  integrate scRNA-seq samples with monocle, perform basic analyses
  and extract resulting metadata
- [correct_ambiance.R](correct_ambiance.R) –
  remove cell-free RNA contamination
- [classify_cell_types.R](classify_cell_types.R) -
  perform cell type classification via SingleR
- [assemble_metadata.R](assemble_metadata.R) -
  generate one CSV and RDS file with all metadata
- [analyse_dge.R](analyse_dge.R) -
  analyse differential gene expression using mixed models
- [analyse_cnv.R](analyse_cnv.R) -
  analyse copy number variations
- [analyse_ccc.R](analyse_ccc.R) -
  analyse cell-cell communication
- [analyse_myeloid.R](analyse_myeloid.R) -
  analyse myeloid subpopulation
- [prepare_data_dong.R](prepare_data_dong.R) –
  prepare dataset by Dong et al for analysis
- [classify_as_adrenal.R](classify_as_adrenal.R) –
  classify tumor cells as adrenal medullary cell types
- [compare_tumors_pseudobulk.R](compare_tumors_pseudobulk.R) –
  comparison of tumor samples via pseudobulk correlation



## Plotting functions

Run these R scripts in arbitrary order to generate publication figures and tables:

- [plot_figure_1_S1.R](plot_figure_1_S1.R) - includes Tables S2, S3, and S10
- [plot_figure_2_S2.R](plot_figure_2_S2.R)
- [plot_figure_3_S3.R](plot_figure_3_S3.R) - includes Table S4
- [plot_figure_4_S4_S7b.R](plot_figure_4_S4_S7b.R) - includes Tables S5 and S6
- [plot_figure_S5_S7c.R](plot_figure_S5_S7c.R) - includes Table S7 and S8
- [plot_figures_revision.R](plot_figures_revision.R) – plots for the reply to reviewers



## Other scripts

- [common_functions.R](common_functions.R) -
  functions used throughout the project
- [plot_dependencies.R](plot_dependencies.R) –
  plot the dependency graph
- [styling.R](styling.R) –
  functions for generating publication-quality figures and tables
