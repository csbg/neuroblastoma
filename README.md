# Analysis of scRNA-seq data from bone marrow samples of ultra-high-risk neuroblastoma patients

## Folders

(Not all of these folders appear in the git repository.)

- `data_generated`: output files generated by the scripts in this repository.
- `data_raw`: various raw input data
- `data_wip`: data from work in progress
- `doc`: project documentation
- `literature`: relevant publications
- `metadata`: additional information; e.g., `sample_groups.csv` contains mappings from samples to neuroblastoma groups.
- `plots`: generated plots
- `renv`: R environment data
- `tables`: exported tables



## Main workflow

![workflow](plots/dependency_graph.png)

- [perform_qc.R](perform_qc.R) -
  perform QC filtering, ensure unique cell names
- [integrate_rna.R](integrate_rna.R) -
  integrate scRNA-seq samples with monocle, perform basic analyses
  and extract resulting metadata
- [correct_ambiance.R](correct_ambiance.R) –
  remove cell-free RNA contamination
- [classify_cell_types.R](classify_cell_types.R) -
  perform cell type classification via SingleR;
  invoke via [sge_classify_cell_types.sh](sge_classify_cell_types.sh)
- [assemble_metadata.R](assemble_metadata.R) -
  generate one CSV and RDS file with all metadata
- [analyse_rna.R](analyse_rna.R) -
  analysis of integrated data
- [plot_markers.R](plot_markers.R) -
  plot manually selected canonical cell type markers
- [analyse_dge_pb.R](analyse_dge_pb.R) -
  analyse differential gene expression using pseudobulks
- [plot_dge_pb.R](plot_dge_pb.R) -
  plot pseudobulk DGEA results
- [analyse_dge_mm.R](analyse_dge_mm.R) -
  analyse differential gene expression using mixed models
- [plot_dge_mm.R](plot_dge_mm.R) -
  plot mixed model DGEA results
- [analyse_cnv.R](analyse_cnv.R) -
  analyse copy number variations
- [plot_cnv.R](plot_cnv.R) - 
  plot copy number variations



## Other scripts

- [common_functions.R](common_functions.R) -
  functions used throughout the project
- [plot_dependencies.R](plot_dependencies.R) –
  plot the dependency graph, which is shown above
- [styling.R](styling.R) –
  functions for generating publication-quality figures and tables
- R scripts starting with `wip` represent work in progress
  that has not yet been integrated into the main workflow
