# Download files specified in the variables rna_files and atac_files for each
# sample in parent_url, and save the files in data_dir
#
# @DEPO raw data

library(tidyverse)
library(fs)
library(RCurl)
library(XML)



# Parameters --------------------------------------------------------------

# URL that points to the parent directory
parent_url <- "https://biomedical-sequencing.at/projects/BSA_0407_STM_Neuroblastoma_2ba0210fb73d412397728e8a97a3e423/"

# parent directory where files should be saved on the disk
data_dir <- "data_raw"

# if TRUE, only download files in selected_files
only_selected_files <- FALSE

selected_files <- c(
  # scRNA-seq
  "filtered_feature_bc_matrix.h5",
  "metrics_summary.csv",
  "analysis/tsne/2_components/projection.csv",
  "analysis/clustering/graphclust/clusters.csv",
  
  # scATAC-seq
  "summary.csv",
  "analysis/tsne/2_components/projection.csv",
  "analysis/clustering/graphclust/clusters.csv"
)



# Crawl HTML files --------------------------------------------------------

get_links <- function(url) {
  message("In ", url)
  links <- 
    url %>% 
    getURL() %>% 
    readHTMLTable(skip.rows = 1:2) %>%
    pluck(1) %>%
    as_tibble(.name_repair = "unique") %>%
    filter(!is.na(Name)) %>%
    pull(Name) %>% 
    {str_glue("{url}{.}")}
  
  files <- str_subset(links, "/$", negate = TRUE)
 
  subfiles <- map(
    str_subset(links, "/$"),
    get_links
  ) %>% flatten_chr()
  
  c(files, subfiles)
}

all_files <- get_links(parent_url)

download_data <-
  all_files %>% 
  enframe(value = "url") %>% 
  select(!name) %>% 
  mutate(
    dest =
      url %>%
      str_sub(start = str_length(parent_url)) %>% 
      {str_c(data_dir, .)}
  )

if (only_selected_files) {
  download_data <- 
    download_data %>% 
    rowwise() %>%
    mutate(is_selected = any(str_detect(url, coll(selected_files)))) %>% 
    filter(is_selected)
}

download_data %>% 
  write_delim("bsf_files.txt", col_names = FALSE)


# Download files ----------------------------------------------------------

# download.data() often fails for large files. Thus, it is recommended that you
# use the following bash script:

# while read -r url filename; do
#   curl "$url" --create-dirs -C - -o "$filename"
# done < bsf_files.txt

# pwalk(
#   download_data,
#   function(url, dest, ...) {
#     dir_create(path_dir(dest))
#     if (file_exists(dest) && !overwrite) {
#       message("Skipping ", dest, " since it already exists")
#     } else {
#       message("Downloading ", dest, " from ", url)
#       download.file(url, dest)
#     }
#   }
# )
