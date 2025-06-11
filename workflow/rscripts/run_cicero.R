#!/usr/bin/env Rscript

# ---------------- Check and install Cicero from GitHub if not already installed ---------------



# Check and install monocle3
if (!requireNamespace("monocle3", quietly = TRUE)) {
  # Bioconductor installation
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(version = "3.21")
  BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                          'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                          'SummarizedExperiment', 'batchelor', 'HDF5Array',
                          'ggrastr'))
  
  # devtools and BPCells
  install.packages("devtools")
  remotes::install_github("bnprks/BPCells/r")
  
  # monocle3 installation
  install_dir <- .libPaths()[1]
  devtools::install_github("cole-trapnell-lab/monocle3", lib = install_dir)
}
library(monocle3)

# Check and install cicero
if (!requireNamespace("cicero", quietly = TRUE)) {
  # Bioconductor installation for Cicero
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
  
  # devtools for Cicero
  install.packages("devtools")
  
  # cicero installation
  devtools::install_github("cole-trapnell-lab/cicero-release", ref = "monocle3")
}
library(cicero)


# --------------- Command Line Arguments ---------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: Rscript run_cicero.R <matrix_file> <barcodes_file> <peaks_file> <chrom_sizes_file> <output_folder>")
}

matrix_file <- args[1]
barcodes_file <- args[2]
peaks_file <- args[3]
chrom_sizes_file <- args[4]
output_folder <- args[5]

# --------------- Read in matrix data ---------------
cat("Reading matrix data...\n")
indata <- Matrix::readMM(matrix_file)
indata@x[indata@x > 0] <- 1  # Binarize

# --------------- Read cell and peak info ---------------
cat("Reading cell and peak information...\n")
cellinfo <- read.table(barcodes_file, header = FALSE)
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

peakinfo <- read.table(peaks_file, header = FALSE)
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep = "_")
row.names(peakinfo) <- peakinfo$site_name

# --------------- Format matrix ---------------
cat("Formatting matrix data...\n")
row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# --------------- Create CDS object ---------------
cat("Creating CellDataSet...\n")
input_cds <- suppressWarnings(new_cell_data_set(
  indata,
  cell_metadata = cellinfo,
  gene_metadata = peakinfo
))

# --------------- Preprocessing ---------------
set.seed(2017)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# --------------- Dimensional Reduction (UMAP) ---------------
input_cds <- reduce_dimension(input_cds, reduction_method = "UMAP", preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

# --------------- Create Cicero CDS ---------------
cat("Creating Cicero CDS...\n")
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# --------------- Load chromosome lengths ---------------
cat("Reading chromosome lengths from provided file...\n")
chromosome_length_raw <- read.table(chrom_sizes_file, header = FALSE, stringsAsFactors = FALSE)
chromosome_length <- chromosome_length_raw[, 1:2]
colnames(chromosome_length) <- c("chr", "length")

# --------------- Run Cicero ---------------
cat("Running Cicero...\n")
conns <- run_cicero(cicero_cds, chromosome_length)

# --------------- Save Results ---------------
cat("Saving results...\n")
all_peaks <- row.names(exprs(input_cds))
write.csv(all_peaks, file = file.path(output_folder, "all_peaks.csv"), row.names = FALSE)
write.csv(conns, file = file.path(output_folder, "cicero_connections.csv"), row.names = FALSE)

cat("Cicero analysis completed successfully!\n")
