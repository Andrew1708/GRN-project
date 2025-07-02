#!/usr/bin/env Rscript

# ---------------- Check and install Cicero from GitHub if not already installed ---------------

library(monocle3)
library(cicero)

log <- function(msg) message("[INFO] ", msg)


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
log("Reading matrix data...\n")
indata <- Matrix::readMM(matrix_file)
indata@x[indata@x > 0] <- 1  # Binarize

# --------------- Read cell and peak info ---------------
log("Reading cell and peak information...\n")
cellinfo <- read.table(barcodes_file, header = FALSE)
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# Read the first line
first_line <- readLines(peaks_file, n = 1)
# Split the line and test if column 2 is an integer (no header if it's numeric)
tokens <- strsplit(first_line, "\t")[[1]]
has_header <- !grepl("^\\d+$", tokens[2])

# Read the file accordingly
if (has_header) {
  log("Header detected in peaks file. Reading with header=TRUE\n")
  peakinfo <- read.table(peaks_file, header = TRUE)
} else {
  log("No header detected in peaks file. Reading with header=FALSE\n")
  peakinfo <- read.table(peaks_file, header = FALSE)
}

names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep = "_")
row.names(peakinfo) <- peakinfo$site_name

# --------------- Format matrix ---------------
log("Formatting matrix data...\n")
row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# --------------- Create CDS object ---------------
log("Creating CellDataSet...\n")
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
log("Creating Cicero CDS...\n")
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# --------------- Load chromosome lengths ---------------
log("Reading chromosome lengths from provided file...\n")
chromosome_length_raw <- read.table(chrom_sizes_file, header = FALSE, stringsAsFactors = FALSE)
chromosome_length <- chromosome_length_raw[, 1:2]
colnames(chromosome_length) <- c("chr", "length")

# --------------- Run Cicero ---------------
log("Running Cicero...\n")

conns <- run_cicero(cicero_cds, chromosome_length)

# --------------- Save Results ---------------
log("Saving results...\n")
all_peaks <- row.names(exprs(input_cds))
write.csv(all_peaks, file = file.path(output_folder, "all_peaks.csv"), row.names = FALSE)
write.csv(conns, file = file.path(output_folder, "cicero_connections.csv"), row.names = FALSE)

log("Cicero analysis completed successfully!\n")
