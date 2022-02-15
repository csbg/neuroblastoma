# Creates the following files:
# * metadata.rds - combined metadata
# * celltype_details.rds - processed SingleR details
#
# @DEPI metadata_monocle.csv
# @DEPI cell_types_[ref]_[labels].csv
# @DEPO metadata.rds
# @DEPO celltype_details.rds

library(ontoProc)
library(igraph)
library(celldex)
library(tidyverse)
library(fs)
source("common_functions.R")



# Functions ---------------------------------------------------------------

#' Load cell type data from SingleR.
#'
#' @param file SingleR results file (`cell_types_singler_[ref]_[labels].csv`).
#'
#' @return A dataframe with one row per cell and 12 columns:
#'   * `cell`: barcode
#'   * `first_labels`,
#'     `tuning_scores_first`,
#'     `tuning_scores_second`, and
#'     `labels`: reused from the output of `SingleR::SingleR()`
#'   * `median_score`: median of cell type-specific correlation scores
#'   * `diff_next`: difference between the best score and next best score
#'   * `label_score`: score for the assigned label
#'   * `delta_score`: difference between `label_score` and `median_score`
#'   * `median_delta_score`: median of `delta_score` across all cells
#'                           of the assigned type
#'   * `mad_delta_score`: mean absolute deviation
#'   * `z_score`: `delta_score` normalized via median and MAD
load_singler_data <- function(file) {
  info("Loading {file}")
  df <- 
    read_csv(file) %>%
    rowwise() %>% 
    mutate(median_score = median(c_across(starts_with("score_")))) %>% 
    ungroup() %>% 
    mutate(diff_next = tuning_scores_first - tuning_scores_second) %>%
    pivot_longer(
      starts_with("score_"),
      names_to = "scores",
      names_prefix = "score_",
      values_to = "label_score"
    ) %>%
    filter(scores == labels) %>%
    mutate(delta_score = label_score - median_score) %>%
    select(!c(pruned_labels, scores))
  
  left_join(
    df,
    df %>%
      group_by(labels) %>%
      summarise(
        median_delta_score = median(delta_score),
        mad_delta_score = mad(delta_score)
      ),
    by = "labels"
  ) %>%
    mutate(z_score = (delta_score - median_delta_score) / mad_delta_score)
}


#' Load cell metadata from several CSV files.
#' 
#' @param files CSV files with metadata.
#'
#' @return A dataframe with metadata combined from several files.
load_cell_metadata <- function(files) {
  sample_order <-
    read_csv("metadata/sample_groups.csv", comment = "#") %>%
    arrange(group, sample) %>% 
    pull(sample) %>% 
    unique()
  
  map(files, read_csv) %>% 
    reduce(left_join, by = "cell") %>%
    mutate(
      group =
        as_factor(group) %>%
        fct_relevel("I", "II", "III", "IV"),
      sample =
        as_factor(sample) %>% 
        fct_relevel(sample_order),
      across(
        matches("cluster|partition"),
        ~as_factor(.x) %>% fct_inseq()
      )
    )
}


#' Combine general metadata and cell types.
#'
#' @param df_metadata Dataframe returned by `load_cell_metadata()`.
#' @param singler_list Named list of dataframes as returned by
#'   `load_singler_data()`.
#' @param min_z_score Minimum z score, ...
#' @param min_delta_score minimum delta score, ...
#' @param min_diff_next ... and minimum `diff_next` for a cell type label to be
#'   retained. Other labels are set to `NA`. This filtering essentially
#'   reproduces the functionality of `SingleR::pruneScores()`.
#'
#' @return The dataframe provided by `df_metadata`, with two additional columns
#'   `cell_type_fine` and `cell_type_broad`.
add_cell_types <- function(df_metadata,
                           singler_list,
                           min_z_score = -3,
                           min_delta_score = -Inf,
                           min_diff_next = 0) {
  df_cell_types <- 
    imap(
      singler_list,
      function(df, file) {
        ref <- str_match(file, "cell_types_(.*)_")[, 2]
        label <- ifelse(str_detect(file, "main"), "broad", "fine")
        type_name <- str_glue("cell_type_{ref}_{label}")
        
        df %>% 
          transmute(
            cell = cell,
            {{type_name}} :=
              case_when(
                z_score >= min_z_score &
                  delta_score >= min_delta_score &
                  diff_next > min_diff_next
                ~ labels,
                TRUE
                ~ NA_character_
              ) %>% 
              as_factor() %>% 
              fct_infreq()
          )
      }
    ) %>% 
    reduce(left_join, by = "cell")

  left_join(df_metadata, df_cell_types, by = "cell")
}



#' Label all cells in a cluster with the same cell type. To this end, perform a
#' majority vote on cell ontology IDs associated with all cell type labels in
#' the respective cluster. Use a set of top-level CO IDs to which all fine cell
#' type labels are assigned.
#'
#' @param df_metadata Dataframe returned by `load_cell_metadata()`.
#' @param clusters Column with cluster IDs.
#' @param ancestors Character vector with cell ontology IDs.
#' @param abbrev Character vector with cell type abbreviations
#'
#' @return The dataframe provided by `df_metadata`, with four additional columns
#'   `cellont_name`, `cellont_abbr`, `cellont_id`, and `cellont_cluster`.
add_unified_labels <- function(df_metadata, clusters, ancestors, abbrev) {
  # get cell ontology IDs for cell type reference datasets
  reference_cell_types <-
    list(
      hpca = HumanPrimaryCellAtlasData(),
      blueprint = BlueprintEncodeData(),
      dice = DatabaseImmuneCellExpressionData(),
      dmap = NovershternHematopoieticData(),
      monaco = MonacoImmuneData()
    ) %>% 
    map_dfr(~as_tibble(colData(.)), .id = "ref") %>% 
    select(ref, cell_type = label.fine, coid = label.ont) %>% 
    distinct()
  
  # add CO IDs of assigned cell types to each cell
  data_cellont <-
    df_metadata %>% 
    select(cell, cluster_col = {{clusters}}, ends_with("fine")) %>%
    pivot_longer(
      !c(cell, cluster_col),
      names_to = "ref",
      values_to = "cell_type",
      names_pattern = "type_(.*)_fine"
    ) %>%
    left_join(reference_cell_types, by = c("ref", "cell_type"))
  
  # make DAG from cell ontology and assemble all childs of selected ancestors
  cell_ontology <- getCellOnto()
  
  cell_ontology_graph <- 
    cell_ontology$children %>% 
    imap_dfr(~tibble(from = .y, to = .x)) %>% 
    graph_from_data_frame()
  
  common_ancestors <-
    tibble(ancestor = ancestors) %>% 
    rowwise() %>% 
    mutate(
      cellont_name = cell_ontology$name[ancestor],
      child =
        cell_ontology_graph %>%
        subcomponent(ancestor, mode = "out") %>% 
        magrittr::use_series("name") %>% 
        list()
    ) %>% 
    unnest_longer(child)
  
  # unify labels by cluster-wise majority vote
  unified_labels <- 
    data_cellont %>% 
    left_join(common_ancestors, by = c(coid = "child")) %>% 
    replace_na(list(cellont_name = "other")) %>%
    group_by(cluster_col, cellont_name) %>%
    summarise(cellont_id = first(ancestor), n = n()) %>%
    slice_max(n, n = 1) %>% 
    select(!n) %>% 
    ungroup() %>% 
    mutate(
      cellont_abbr =
        cellont_name %>%
        as_factor() %>% 
        fct_recode(!!!cell_type_abbreviations) %>%
        fct_relevel(names(cell_type_abbreviations)),
    ) %>%
    arrange(cellont_abbr, cluster_col) %>% 
    mutate(
      cellont_cluster =
        str_glue("{cellont_abbr} ({cluster_col})") %>%
        as_factor()
    )
  join_by <- set_names("cluster_col", rlang::as_string(rlang::enexpr(clusters)))
  
  df_metadata %>%
    left_join(unified_labels, by = join_by)
}



#' Split cluster 8 into two clusters, denoted 8 and 21.
#' 
#' This occurs at higher resolution (cluster_20), where cluster 8 splits into
#' 5+32 (NB cells, only groups II-IV) and 36 (other cells, also group I).
#'
#' @param df_metadata Dataframe returned by `load_cell_metadata()`.
#'
#' @return Modified dataframe.
modify_clusters <- function(df_metadata) {
  df_metadata %>% 
    mutate(
      cluster_50 =
        case_when(
          cluster_20 == "36" ~ "21",
          TRUE ~ as.character(cluster_50)
        ) %>% 
        as_factor() %>% 
        fct_inseq()
    )
}



# Load data ---------------------------------------------------------------

folder <- "data_generated"
metadata_files <- "metadata_monocle.csv"

ancestors <- c(
  "CL:0000576",
  "CL:0000084",
  "CL:0000236",
  "CL:0000784",
  "CL:0000764",
  "CL:0000623",
  "CL:0000540",
  "CL:0008001"
)

cell_type_abbreviations <- c(
  "T" = "T cell",
  NK  = "natural killer cell",
  B   = "B cell",
  M   = "monocyte",
  pDC = "plasmacytoid dendritic cell",
  E   = "erythroid lineage cell",
  SC  = "hematopoietic precursor cell",
  NB  = "neuron"
)

singler_data <- map(
  dir_ls(folder, regex = "cell_types"),
  load_singler_data
)

nb_data <-
  load_cell_metadata(str_glue("{folder}/{metadata_files}")) %>% 
  modify_clusters() %>% 
  add_cell_types(singler_data) %>%
  add_unified_labels(cluster_50, ancestors, cell_type_abbreviations)



# Export data -------------------------------------------------------------

nb_data %>% saveRDS(str_glue("{folder}/metadata.rds"))
singler_data %>% saveRDS(str_glue("{folder}/celltype_details.rds"))
