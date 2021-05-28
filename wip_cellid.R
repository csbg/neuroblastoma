# cell type classification via CelliD

library(CellID)
library(scuttle)
library(Seurat)
library(tidyverse)
source("common_functions.R")



# Load data ---------------------------------------------------------------

nb_list <-
  readRDS("data_generated/rna_qcpassed.rds") %>% 
  map(
    SCTransform,
    vars.to.regress = "percent_mito",
    method = "qpoisson",
    verbose = FALSE
  )

markers <-
  read_csv("metadata/microenvironment_markers.csv", comment = "#") %>% 
  as.list() %>%
  map(discard, is.na)



# Analysis ----------------------------------------------------------------

cell_types <- map(
  nb_list,
  function(sample) {
    info("Classifying sample {Project(sample)}")
    sample %>%
      RunMCA() %>%
      RunCellHGT(markers, minSize = 6)
  }
)

cell_type_predictions <-
  map_dfr(
    cell_types,
    function(sample) {
      imap_dfr(
        array_branch(sample, 2),
        ~tibble(
          cell = .y,
          cell_type = names(markers)[which.max(.x)],
          log_q = max(.x)
        )
      )   
    }
  )



# Save data ---------------------------------------------------------------

cell_types %>% saveRDS("data_wip/cellid_celltypes.rds")

cell_type_predictions %>% write_csv("data_wip/cellid_predictions.csv")
