# library(devtools)
# library(reticulate)
# library(rgdal) 
# library(Matrix)
# library(scater)
library(scran)
# library(Seurat)
# library(hdf5r)
library(enrichR)
# library(scales)
library(ComplexHeatmap)
library(latex2exp)
# library(log4r)
# library(magick)
# library(colorspace)
library(RColorBrewer)
# library(grid)
# library(ggpubr)
# library(rstatix)
library(data.table)
library(monocle3) 
library(CytoTRACE)
library(readxl)
library(fs)
library(tidyverse)



# # set environment paths (for whole environment)
# paths <- read.csv("analysis_tumor/0_paths.csv", header = FALSE)
# 
# path <- paths$V2 #ubuntu data
# path <- paths$`3` #windows data
# path <- paths$V7 # CAME data 
# 
# path_plot <- paths$V3 #ubuntu plots
# path_plot <- paths$`4` #windows data
# path_plot <- paths$V8 # CAME plots
# 
# path_DATA <- paths$`5` #win GFS W.E.S data
# path_DATA <- paths$V4 #ubuntu GFS W.E.S data
# path_DATA <- paths$V6 # CAME WES GFS 