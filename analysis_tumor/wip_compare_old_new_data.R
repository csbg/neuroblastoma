source("analysis_tumor/packages.R")

data_old <- readRDS("analysis_tumor/data_generated/NB_CDS_trajectory.rds")
data_new <- cds_trajectory

data_old
data_new

colData(data_old)
colData(data_new)

plot_cells(data_old) + coord_fixed()
plot_cells(data_new) + coord_fixed()




waldo::compare(data_old, data_new, tolerance = 0.00000001)

