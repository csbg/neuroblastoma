################################################################################
# Gene Quality controlSubset cds to get rid of "noisy"/low expression 
# genes for DE analysis 
################################################################################
# load helper functions
source("nb_Packages.R") 
source("nb_common_functions.R") 

# import data
cds <- readRDS(paste0(path,"/nb_CDS.rds"))

#detect genes with 0 cells expressing it in subset
cds <- detect_genes(cds)

#WES simpler approach: get column as vector, use purrr::keep() to remove zeros,
#WES use names of this vector for subsetting 
cutoff <- 0

sub1 <-
  rowData(cds)$num_cells_expressed %>% 
  keep(~. > cutoff) %>% 
  names()

#remove genes expressed in 0 cells from CDS
cds <- cds[sub1, ]

# get ribosomal genes
sub2 <- !str_detect(pattern="^RP[LS]\\d", rownames(cds))

#subset CDS to data w.o. Ribosomal genes
cds <- cds[sub2, ] 

#remove genes below 20% relative global expression for Regression analysis 
#   = Genes expressed in less than 20% of cells/cluster --> remove
#WES see detect_genes_clusterwise() in nb_HelperFunction.R,
#WES which replaces geneQC()

sub3 <- 
  detect_genes_clusterwise(cds) %>% 
  group_by(gene) %>%
  summarise(rel_expr = max(rel_expr)) %>% 
  filter(rel_expr > 0.2) %>% 
  pull(gene)

cds <- cds[sub3, ]

# collect all removed genes
sub <- list("zero expression"=sub1, 
            "ribosomal genes"=sub2, 
            "low global expression"=sub3
            )

# detect differentially expressed genes using helper function
DE <- pairwiseDE(cds)

# select significant genes for each cluster
DEG <-
  levels(clusters(cds)) %>%
  map_dfr(~get_significant_genes(DE, .))

# export data
saveRDS(DEG,file = paste0(path,"/nb_DEG.rds"))
saveRDS(sub,file = paste0(path,"/nb_removedGenes.rds"))
