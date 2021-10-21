#### Packages:
##################################################
#Set environment paths (set here for whole environment)
paths <- read.csv("nb_Paths.csv", 
                  header = F
)

path <- paths$V2 #ubuntu data
path <- paths$`3` #windows data
path <- paths$V7 # CAME data 

path_plot <- paths$V3 #ubuntu plots
path_plot <- paths$`4` #windows data
path_plot <- paths$V8 # CAME plots

path_DATA <- paths$`5` #win GFS W.E.S data
path_DATA <- paths$V4 #ubuntu GFS W.E.S data
path_DATA <- paths$V6 # CAME WES GFS 


library(devtools)
library(reticulate)
library(readr)
library(stringr) 
library(rgdal) 
library(Matrix)
library(scater)
library(scran)
library(Seurat)
library(hdf5r)
library(enrichR) 
library(ggplot2)
library(scales)
library(ComplexHeatmap)
library(WriteXLS)
library(latex2exp)
library(ggplot2)
library(log4r)
library(magick)
library(fs)
library(colorspace)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(rstatix)
library(data.table)
library(fs)
library(tidyverse)

# library(org.Hs.eg.db)
# library(MyPlotTrajectoryPackage)
################################################################################
# Installing monocle3:

# renv::install('bioc::BiocGenerics', 'bioc::DelayedArray', 'bioc::DelayedMatrixStats', 'bioc::limma', 'bioc::S4Vectors', 'bioc::SingleCellExperiment', 'bioc::SummarizedExperiment', 'bioc::batchelor', 'bioc::Matrix.utils')
# renv::install('cole-trapnell-lab/leidenbase')
# renv::install('cole-trapnell-lab/monocle3')
################################################################################
#Installing cytoTRACE:
# dependencies:

# renv::install("bioc::genefilter")
# renv::install("bioc::sva")

# Cellular (Cyto) Trajectory Reconstruction Analysis using gene Counts and Expression 
# is a computational method that predicts the differentiation state of cells from 
# single-cell RNA-sequencing data

# renv::install("/home/ubuntu/mnt/agfortelny/PROJECTS/Neuroblastoma/analysis/samuel/Data/CytoTRACE_0.3.3.tar.gz")

# for multiple dataset comparison use rminiconda, reticulate to use python packages scanoramaCT & numpy 
py <- rminiconda::find_miniconda_python("my_python")
#rminiconda::rminiconda_pip_install("scanoramaCT","my_python")
reticulate::use_python(py, required = TRUE)

library(CytoTRACE)
library(monocle3) 

#renv::install("efheuston/MyPlotTrajectoryPackage")

#renv::install("bioc::org.Ce.eg.db","bioc::org.Hs.eg.db", "bioc::org.Mm.eg.db")
#renv::install("wjawaid/enrichR")
# renv::install("bioc::org.Hs.eg.db")
 