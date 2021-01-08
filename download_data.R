# Download files specified in the variables rna_files and atac_files for each
# sample in parent_url, and save the files in data_dir

library(tidyverse)
library(fs)
library(RCurl)
library(XML)

parent_url <- "https://biomedical-sequencing.at/projects/BSA_0407_STM_Neuroblastoma_2ba0210fb73d412397728e8a97a3e423/COUNT/"
count_dir <- getURL(parent_url)
data_dir <- "data_raw"

rna_files <- c(
  # "raw_feature_bc_matrix.h5",
  "filtered_feature_bc_matrix.h5",
  "metrics_summary.csv",
  "analysis/tsne/2_components/projection.csv",
  "analysis/clustering/graphclust/clusters.csv"
)
rna_folder_re <- "[^L]_transcriptome/"  # exclude aggregated samples

atac_files <- c(
  # "raw_peak_bc_matrix.h5",
  # "filtered_peak_bc_matrix.h5",
  # "filtered_tf_bc_matrix.h5",
  "summary.csv",
  "analysis/tsne/2_components/projection.csv",
  "analysis/clustering/graphclust/clusters.csv"
)
atac_folder_re <- "ATAC/"

download_cellranger_files <- function(folder_re, files, overwrite = FALSE) {
  folders <- 
    readHTMLTable(count_dir)[[1]] %>% 
    as_tibble(.name_repair = "unique") %>%
    filter(str_detect(Name, folder_re)) %>%
    pull(Name)
  
  walk(
    folders,
    function(folder) {
      walk(
        files,
        function(file) {
          url <- str_c(parent_url, folder, file)
          dest <- path_join(c(data_dir, folder, file))
          dir_create(path_dir(dest))
          if (file_exists(dest) && !overwrite)
            message("Skipping ", dest, " since it already exists")
          else {
            message("Downloading ", dest, " from ", url)
            download.file(url, dest)
          }
        }
      )
    }
  )
}

download_cellranger_files(rna_folder_re, rna_files)
download_cellranger_files(atac_folder_re, atac_files)
