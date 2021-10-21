source("analysis_tumor/packages.R") 
source("analysis_tumor/common_functions.R") 

# import data
cds <- readRDS("analysis_tumor/data_generated/cds.rds")

# detect genes with 0 cells expressing it in subset
cds <- detect_genes(cds)

# remove genes expressed in 0 cells from CDS
cutoff <- 0
sub1 <- rowData(cds)$num_cells_expressed > cutoff
cds <- cds[sub1, ]

# remove ribosomal genes
sub2 <- !str_detect(rownames(cds), "^RP[LS]\\d")
cds <- cds[sub2, ] 

# remove genes below 20% relative global expression
# i.e., genes expressed in less than 20% of cells per cluster
sub3 <- 
  detect_genes_clusterwise(cds) %>% 
  group_by(gene) %>%
  summarise(rel_expr = max(rel_expr)) %>% 
  filter(rel_expr > 0.2) %>% 
  pull(gene)
cds <- cds[sub3, ]

# collect all removed genes
sub <- list(
  "zero expression" = sub1, 
  "ribosomal genes" = sub2, 
  "low global expression" = sub3
)

# detect differentially expressed genes
DE <- pairwiseDE(cds)

# select significant genes for each cluster
DEG <-
  levels(clusters(cds)) %>%
  map_dfr(~get_significant_genes(DE, .))

# export data
saveRDS(DEG, "analysis_tumor/data_generated/deg.rds")
saveRDS(sub, "analysis_tumor/data_generated/removed_genes.rds")
