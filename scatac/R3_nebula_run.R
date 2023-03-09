# Set working directory
setwd('/media/AGFORTELNY/people/rohit/')

# Required packages
pacman::p_load(Seurat, 
               SeuratWrappers,
               Signac, 
               SingleCellExperiment, 
               tidyverse, 
               scater,
               nebula)
# File paths
file_path = "projects/nblast_geneActivity.rds" 
anno_path = 'projects/neuroblastoma/idents.txt'

annotation <- function(file_path, anno_path) {
  temp_obj  <- readRDS(file = file_path)
  
  # Add patient group information
  groups <- list(c('2', '3', '6', '15', '16'),
                 c('1', '4', '5', '11'),
                 c('7','12'),
                 c('8', '9', '10','13', '14'))
  temp_obj$group <- temp_obj@meta.data %>% 
    dplyr::mutate(
      group = dplyr::case_when(
        orig.ident %in% groups[[1]] ~"I",
        orig.ident %in% groups[[2]] ~"II", 
        orig.ident %in% groups[[3]] ~"III",
        orig.ident %in% groups[[4]] ~"IV")) %>% 
    pull(group)
  temp_obj$activ.ident <- temp_obj@active.ident
  
  # Change Seurat identities
  annotate <- read.table(anno_path,
                         sep = ',',
                         header = FALSE, 
                         strip.white = TRUE)
  annot <- as.character(annotate$V2)
  names(annot) <-annotate$V1
  temp_obj <- Seurat::RenameIdents(temp_obj, annot)
  temp_obj$activ.ident <- Idents(temp_obj)
  
  #
  return(temp_obj)
}

nb <- annotation(file_path, anno_path)

run_nebula <- function(cell_type) {
  cell_typ <- subset(x = nb, 
                     subset = activ.ident == cell_type)
  path = '/media/AGFORTELNY/people/rohit/projects/clusterObj'
  myfile <- file.path(path, paste0(cell_type,'.csv'))
  sce <- SeuratWrappers::as.cell_data_set(cell_typ)
  sce <- logNormCounts(sce) 
  model_mat <- model.matrix(~group , data = colData(sce))
  model_mat <- model_mat[, colSums(model_mat != 0) > 0] 
  sce_grouped <- group_cell(counts(sce), 
                            id = colData(sce)$orig.ident, 
                            pred = model_mat, 
                            offset = NULL)
  result <- nebula(sce_grouped$count, 
                   sce_grouped$id, 
                   pred=sce_grouped$pred)
  
  result$summary <- result$summary %>%
    as_tibble() %>% 
    mutate(
      algorithm = result$algorithm,
      convergence = result$convergence,
      overdispersion_subject = result$overdispersion$Subject,
      overdispersion_cell = result$overdispersion$cell
    )
  closest_genes <- ClosestFeature(cell_typ, 
    regions = result$summary$gene) %>% 
    mutate(gene = query_region)

  results_with_closest_gene <- merge(result$summary, 
                                     closest_genes, by = 'gene')
  write.table(results_with_closest_gene, file = myfile)
  #return(results_with_closest_gene)
}

set.seed(1)
dge_mm_nebula <- map_dfr(
  nb$activ.ident %>%
    unique() %>%
    setdiff("NB") %>%
    set_names(),
  run_nebula,
  .id = "cell_type"
)

datalist = list()

nebula_results_df <- function(filename) {
  cell_type <- str_replace(
    basename(filename), ".csv", "")
  results_with_closest_gene <- read.table(filename, header= TRUE, sep = ' ')
  results_with_closest_gene <- results_with_closest_gene %>% 
    mutate(q_groupII = qvalue::qvalue(results_with_closest_gene$p_groupII)$qvalues) %>% 
    mutate(q_groupIII = qvalue::qvalue(results_with_closest_gene$p_groupIII)$qvalues) %>%
    mutate(q_groupIV = qvalue::qvalue(results_with_closest_gene$p_groupIV)$qvalues) %>% 
    mutate(logFC_groupII = logFC_groupII/log(2),
           logFC_groupIII = logFC_groupIII/log(2),
           logFC_groupIV = logFC_groupIV/log(2)) %>% 
    filter(!-log10(q_groupII) > 10) %>% 
    filter(!-log10(q_groupIII) > 10) %>%
    filter(!-log10(q_groupIV) > 10) %>%
    mutate(significanceII = case_when(
      (logFC_groupII < -1 & q_groupII < 0.05 )~ 'Down',
      (logFC_groupII > 1 & q_groupII < 0.05 )~ 'Up', TRUE ~ 'No change')) %>%
    mutate(significanceIII = case_when(
      (logFC_groupIII < -1 & q_groupIII < 0.05 )~ 'Down',
      (logFC_groupIII > 1 & q_groupIII < 0.05 )~ 'Up', TRUE ~ 'No change')) %>%
    mutate(significanceIV = case_when(
      (logFC_groupIV < -1 & q_groupIV < 0.05 )~ 'Down',
      (logFC_groupIV > 1 & q_groupIV < 0.05 )~ 'Up', TRUE ~ 'No change')) %>%
    mutate(delabel = case_when(significanceIII != 'No change' ~ gene_name,
                               significanceII != 'No change' ~ gene_name,
                               significanceIV != 'No change' ~ gene_name)) %>% 
    mutate(cell_type = paste0('Cluster:',cell_type))
  
  return(results_with_closest_gene)
}
# Test case ---------------------------------------------------------------

testcase <- nebula_results_df(filepath)
testcase %>% ggplot(aes(x = logFC_groupIII, 
                        y = -log10(q_groupIII), color = significanceIII)) +
  geom_point(size = 0.5) + #geom_text(aes(label=delabel), size = 2) + 
  labs(title = "B-cell", x = 'LogFC_ATRX_MYCN', y = '-Log10(q_value)') +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-2, 2), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  ggpubr::theme_pubclean() +
  theme(legend.title = element_blank()) 


