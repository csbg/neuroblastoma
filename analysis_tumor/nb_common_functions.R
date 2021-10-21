# Different helper functions and variables
################################################################################

#Pairwise differential expression testing:
# a data.table version from NF
pairwiseDE_old <- function(cds){
  cds@colData$cluster <- cds@clusters$UMAP$clusters
  (n = length(unique(cds@colData$cluster)))
  DEG2 <- data.table()
  i1 <- 1  # Comparison pair (start=1)
  i2 <- 2  # Comparison pair (start=2)
  for (i1 in 1:(n - 1)) {
    for (i2 in (i1 + 1):n) {
      message(paste("Comparing", i1, "to", i2))  # added by WES
      vcds <- cds[,as.numeric(cds@colData$cluster) %in% c(i1,i2)]
      DEG <- data.table(i1 = i1, i2 = i2, data.table(
        coefficient_table(fit_models(
          vcds, model_formula_str = "~cluster + sample + phase", 
          cores = 1, reduction_method = "UMAP"))))
      DEG2 <- rbind(
        DEG2,
        DEG[,c("i1", "i2", "gene_short_name", "term", "estimate", "q_value")]
        )
      DEG <- data.table()
      closeAllConnections()
    }
    closeAllConnections()
  }
  DEG2 %>% filter(grepl(DEG2$term,pattern = "^cluster"))
}

#Pairwise differential expression testing:
#WES A tidyverse version
pairwiseDE <- function(cds) {
  cds@colData$cluster <- clusters(cds)
  list(i1 = levels(clusters(cds)), i2 = levels(clusters(cds))) %>%
    cross_df(`<=`) %>%
    pmap_dfr(
      function(i1, i2) {
        message(paste("Comparing", i1, "to", i2))
        de_genes <-
          cds[, cds@colData$cluster %in% c(i1, i2)] %>% 
          fit_models(model_formula_str = "~cluster + sample + phase") %>% 
          coefficient_table() %>% 
          mutate(i1 = i1, i2 = i2) %>% 
          select(i1, i2, gene_short_name, term, estimate, q_value)
        closeAllConnections()
        de_genes
      }
    ) %>% 
    filter(str_starts(term, "cluster"))
}

################################################################################
#Gene Quality control
# calculate relative global expression of gene 
#   (how many cells express this gene in a given cluster)
# used to cutoff genes which fit the criteria that at least 
#     20% of cells in every clusters express it


#WES Let's make it shorter using a mapping function from purrr.
#WES According to the tidyverse style guide, function names should be verbs.
#WES I also like to document my functions in a standardized way, since I tend
#WES to forget what their arguments and output are. Roxygen comments, even if
#WES they actually were developed for documenting functions in packages,
#WES are also suitable for documenting any function. Insert such a docstring
#WES by moving the mouse cursor into you function code and choosing
#WES Code -> Insert roxygen skeleton.

#' Calculate relative global expression of each gene in a cds object.
#' 
#' Determine how many cells express a gene in a given cluster. This value may be
#' used to cutoff genes which fit the criteria that at least 20% of cells in
#' every clusters express it.
#'
#' @param cds A monocle cell data set.
#'
#' @return A dataframe with columns `gene`, `n_cells_expr` (number of cells
#'   expressing this gene in the respective cluster), `rel_expr` (expression
#'   relative to the number of cells in the cluster), and `cluster`.
detect_genes_clusterwise <- function(cds) {
  map_dfr(
    levels(clusters(cds)),
    function(cluster_id) {
      # subset of counts matrix for current cluster
      cluster_counts <- counts(cds)[, clusters(cds) == cluster_id]
      
      # number of cells that express the respective gene
      n_cells_expr <- rowSums(cluster_counts > 0)
      
      # number of cells in the current cluster
      n_cells <- sum(clusters(cds) == cluster_id)
      
      n_cells_expr %>% 
        enframe(name = "gene", value = "n_cells_expr") %>% 
        mutate(
          rel_expr = n_cells_expr / n_cells,
          cluster = cluster_id
        )
    }
  )  
}

################################################################################
# Enrichment function adapted from W.E.S
enrich_genes <- function(gene_list, databases) {
  Sys.sleep(1)
  enrichr(gene_list, databases) %>% 
    discard(~nrow(.) == 0) %>%
    bind_rows(.id = "db") %>%
    as_tibble()
}

################################################################################
# helper function adapted from common functions @W.E.S github
ggsave_default <- function(filename,
                           type = "png",
                           plot = NULL,
                           width = 297,
                           height = 210,
                           crop = TRUE,
                           ...) {
  if (is.null(filename))
    return()
  
  filename <- fs::path(path_plot, filename, ext = type)
  filename %>%
    fs::path_dir() %>%
    fs::dir_create()
  
  if (is.null(plot)) {
    ggplot2::ggsave(filename, dpi = 300, units = "mm", limitsize = FALSE,
                    width = width, height = height, ...)  
  } else {
    rlang::exec(type, filename, res = 300, units = "mm",
                width = width, height = height,  ...)
    print(plot)
    dev.off()
  }
  
  if (crop)
    knitr::plot_crop(filename)
  
  invisible(filename)
}

################################################################################

#' Get significant genes from differential expression analysis.
#'
#' @param de_results Results from `pairwiseDE()`.
#' @param cluster Name of the cluster for which results should be summarized.
#' @param min_occurrence Minimum number of comparisons where a gene should
#'   appear a significant in order to be classified as marker.
#'
#' @return A data frame wit four columns `gene`, `cluster`, `estimate`, and
#'   `q_value`.
get_significant_genes <- function(de_results, cluster, min_occurrence = 2) {
  i1_genes <- 
    de_results %>% 
    filter(i1 == cluster)
  
  i2_genes <- 
    de_results %>% 
    filter(i2 == cluster) %>% 
    mutate(estimate = -estimate)
  
  bind_rows(i1_genes, i2_genes) %>% 
    filter(q_value < 0.05, estimate > 0) %>% 
    group_by(gene_short_name) %>%
    filter(n() >= min_occurrence) %>%
    summarise(
      cluster = {{cluster}},
      estimate = mean(estimate),
      q_value = mean(q_value)
    ) 
}


#Enrichr databases character vector

dbs <- c("GO_Molecular_Function_2018",                          #1
         "GO_Biological_Process_2018",                          #2
         "GO_Cellular_Component_2018",                          #3
         "WikiPathways_2019_Human",                             #4 Pathway
         "Panther_2016",                                        #5 Pathway
         "KEGG_2019_Human",                                     #6 Pathway
         "NCI-Nature_2016",                                     #7 Pathway       
         "TRRUST_Transcription_Factors_2019",                   #8 TranscriptionFactor 
         "ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X",           #9 TranscriptionFactor 
         "TRANSFAC_and_JASPAR_PWMs",                            #10 TF  
         "MSigDB_Hallmark_2020",                                #11 
         "CORUM",                                               #12 Complex
         "NCI-60_Cancer_Cell_Lines",                            #13 Cancer cell lines
         "Chromosome_Location",                                 #14 Chromosomal location enrichment
         "ClinVar_2019",                                        #15 genomic variation in relationship to human phenotype
         "Human_Gene_Atlas",                                    #16
         "Jensen_TISSUES",                                      #17
         "Jensen_COMPARTMENTS",                                 #18
         "Jensen_DISEASES",                                     #19
         "Reactome_2016",                                       #20 Pathway
         "ChEA_2016",                                           #21 TF
         "ENCODE_TF_ChIP-seq_2015",                             #22 TF
         "Elsevier_Pathway_Collection",                         #23 Pathway
         "MSigDB_Oncogenic_Signatures",                         #24
         "Tissue_Protein_Expression_from_ProteomicsDB",         #25
         "Cancer_Cell_Line_Encyclopedia",                       #26
         "PPI_Hub_Proteins",
         "RNA-Seq_Disease_Gene_and_Drug_Signatures_from_GEO",
         "ARCHS4_Tissues",
         "GWAS_Catalog_2019",
         "DisGeNET",
         "Enrichr_Libraries_Most_Popular_Genes"
         )