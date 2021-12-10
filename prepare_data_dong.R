# Generate a Seurat object with all Dong datasets.
#
# @DEPI published datasets
# @DEPO tumor_data_dong.rds


library(Seurat)
library(tidyverse)
library(fs)



files_dong <- read_csv("metadata/samples_dong.csv", comment = "#")

metadata_dong <-
  read_csv("data_raw/GSE137804/GSE137804_tumor_dataset_annotation.csv") %>%
  separate(cellname, into = c("sample", "cell"), sep = "_")

data_dong <- pmap(
  files_dong %>% filter(high_risk),
  function(file, tumor_id, bad_header, ...) {
    if (bad_header) {
      sce <- read_tsv(
        path_join(c("data_raw", "GSE137804", file)),
        skip = 1,
        col_names =
          read_tsv(path_join(c("data_raw", "GSE137804", file)), n_max = 0) %>%
          colnames() %>%
          prepend("Symbol")
      )
    } else {
      sce <-
        read_tsv(
          path_join(c("data_raw", "GSE137804", file)),
          col_types = cols(
            Gene_ID = "c",
            Symbol = "c",
            .default = col_integer()
          )
        ) %>%
        select(!Gene_ID)
    }

    sce <-
      sce %>%
      mutate(Symbol = make.unique(Symbol)) %>%
      column_to_rownames("Symbol") %>%
      as.matrix() %>%
      CreateSeuratObject(project = tumor_id)

    sce@meta.data <-
      sce@meta.data %>%
      as_tibble(rownames = "cell") %>%
      rename(sample = orig.ident) %>%
      left_join(metadata_dong, by = c("cell", "sample")) %>%
      column_to_rownames("cell")

    sce %>%
      subset(subset = celltype == "tumor")
  }
)

data_dong %>%
  {merge(.[[1]], .[-1])} %>%
  saveRDS("data_generated/tumor_data_dong.rds")

