#!/usr/bin/env Rscript

# Load libraries silently
suppressPackageStartupMessages({
    library(devtools)
    library(remotes)
    library(Cairo)
})


# Install Bioconductor core packages
#BiocManager::install(version = "3.20") #done
#install.packages("MASS", repos = "https://cloud.r-project.org") #done
#install.packages("igraph", repos = "https://cloud.r-project.org") #done



#BiocManager::install(c(
#    "BiocGenerics", #done-no error
#    "DelayedArray" #done-no error
#    "DelayedMatrixStats" #done
#    "limma"  #done
#    "lme4" #done 
#    "S4Vectors" #done-was already installed
#   "SingleCellExperiment" #done
#   "SummarizedExperiment" #done-was already installed
#   "batchelor" #done
#    "HDF5Array" #done
#   "ggrastr" #done
#   "Gviz" #done
#    "GenomicRanges" #done
#   "rtracklayer" #done
##))

# # Install GitHub packages
#remotes::install_github("bnprks/BPCells/r") #done
#remotes::install_github("cole-trapnell-lab/monocle3") #done problmes with units and s2
remotes::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
