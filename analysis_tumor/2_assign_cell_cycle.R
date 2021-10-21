source("analysis_tumor/packages.R") 

cds <- readRDS("analysis_tumor/data_generated/cds_raw.rds")

# cell cycle markers from scran 
hs_pairs <- readRDS(
  system.file("exdata", "human_cycle_markers.rds", package = "scran")
)

set.seed(100)
assignments <- cyclone(
  cds,
  pairs = hs_pairs,
  gene.names = rowData(cds)$ENSEMBL,
  iter = 2500,
  verbose = TRUE
)

# assignments <- readRDS("analysis_tumor/data_generated/cell_cycle_assignments.rds")

# assign cell cycle scores to cds metadata
cds@colData$norm.scores_G1 <- assignments$normalized.scores$G1
cds@colData$norm.scores_G2M <- assignments$normalized.scores$G2M
cds@colData$norm.scores_S <- assignments$normalized.scores$S
cds@colData$phase <- assignments$phases

# export data 
saveRDS(cds, "analysis_tumor/data_generated/cds_raw.rds")
saveRDS(assignments, "analysis_tumor/data_generated/cell_cycle_assignments.rds")
