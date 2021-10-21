source("nb_Packages.R") 
################################################################################
# Detect cell-cycle scores (bootstrap or import of pre-saved file)

# import data
# generated CDS
cds <- readRDS(
  paste0(path,"/nb_CDS_raw.rds")
  ) #change link

# Cell-Cycle DB from scran 
hs.pairs <- readRDS(
  system.file("exdata", 
              "human_cycle_markers.rds", 
              package = "scran"
              ) #change link?
  )

set.seed(100)

assignments <- scran::cyclone(
  cds,
  iter = 2500,
  verbose = TRUE,
  pairs = hs.pairs,
  gene.names = rowData(cds)$ENSEMBL,
  BPPARAM = BiocParallel::MulticoreParam()
)


#WES I used
#assignments <- readRDS("data/nb_cellcycle_assignments.rds")

# assign cell-cycle scores to cds metadata
cds@colData$norm.scores_G1 <- assignments$normalized.scores$G1
cds@colData$norm.scores_G2M <- assignments$normalized.scores$G2M
cds@colData$norm.scores_S <- assignments$normalized.scores$S
cds@colData$phase <- assignments$phases

# export data 
saveRDS(cds, 
        file = paste0(path,"/nb_CDS_raw.rds")
        ) #change link

saveRDS(assignments, 
        file = paste0(path,"/nb_cellcycle_assignments.rds")
        ) #change link
