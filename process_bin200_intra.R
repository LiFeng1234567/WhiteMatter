#!/usr/bin/env Rscript

# ============================================================================
# Script: process_bin200_intra.R
# Purpose: Generate a Seurat object (bin200 level) from Stereo-seq intra‑cellular
#          data restricted to corpus callosum or other regions.
# Usage:   Rscript process_bin200_intra.R <sample_id>
# Example: Rscript process_bin200_intra.R 001
# ============================================================================

# 0. Parse command line argument ---------------------------------------------
args <- commandArgs(TRUE)
if (length(args) < 1) {
  stop("Please provide a sample ID, e.g.: Rscript process_bin200_intra.R 001")
}
sample <- args[1]

# 1. Load required packages --------------------------------------------------
library(dplyr)
library(reshape2)
library(Seurat)

# 2. Define file paths (adjust base directories as needed) -------------------
base_dir      <- "xx/intra"
data_dir      <- "xx/data"
position_dir  <- "xx/relative_position_corpus_callosum_with_annot"

# Input files
counts_file   <- file.path(data_dir, paste0("total_gene_T", sample, "_marmoset_f001_2D_20230130-marmoset-cortex-v2.txt"))
position_file <- file.path(position_dir, paste0(sample, "_relative_position_with_annot.csv"))

# Output file
out_rds       <- file.path(base_dir, paste0(sample, "_intra_bin200.rds"))

# 3. Read and preprocess the spatial transcriptomics data --------------------
cat("Reading counts file:", counts_file, "\n")
data <- read.table(counts_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Create bin200 identifier (200x200 pixel bins, equivalent to ~100 µm)
data$loci <- paste(floor(data$x / 200), floor(data$y / 200), sep = "_")

# Keep only intra‑cellular transcripts (cell_label != 0)
intra_data <- data[data$cell_label != 0, ]

# 4. Read corpus callosum position annotation --------------------------------
cc_pos <- read.csv(position_file, stringsAsFactors = FALSE)
cc_pos$loci <- paste(floor(cc_pos$dim1 / 200), floor(cc_pos$dim2 / 200), sep = "_")

# Check coverage: all CC loci present in our intra‑cellular data?
cc_loci <- unique(cc_pos$loci)
data_loci <- unique(intra_data$loci)
missing <- setdiff(cc_loci, data_loci)
if (length(missing) > 0) {
  warning("Some CC loci have no intra‑cellular data. Missing: ", paste(missing, collapse = ", "))
}
cat("CC loci coverage: ", length(intersect(cc_loci, data_loci)), " / ", length(cc_loci), "\n")

# Subset intra‑cellular data to CC loci only
intra_cc <- intra_data %>% filter(loci %in% cc_loci)

# 5. Aggregate UMI counts per gene per bin -----------------------------------
gene_count <- intra_cc %>%
  group_by(loci, gene) %>%
  summarise(umi_sum = sum(umi_count), .groups = "drop")

# Convert to wide matrix (bins as columns, genes as rows)
wide_mat <- gene_count %>%
  dcast(loci ~ gene, value.var = "umi_sum", fill = 0)

# Set bin names as rownames and remove the 'loci' column
rownames(wide_mat) <- wide_mat$loci
wide_mat <- wide_mat[, -1]   # now a genes × bins matrix

# Transpose to bins × genes (Seurat expects features as rows, cells as columns)
# But CreateSeuratObject counts matrix should be features × cells.
# Our wide_mat has genes as rows, bins as columns – that is correct.
# So we do NOT transpose here.

# 6. Prepare metadata for the Seurat object ----------------------------------
# Average relative position per bin
cc_relative <- cc_pos %>%
  group_by(loci) %>%
  summarise(relative_position = mean(relative_position, na.rm = TRUE), .groups = "drop")

# Majority area per bin (most frequent cc_area)
cc_area <- cc_pos %>%
  group_by(loci, cc_area) %>%
  summarise(Freq = n(), .groups = "drop_last") %>%
  slice_max(order_by = Freq, n = 1, with_ties = FALSE) %>%
  select(loci, area = cc_area)

# Merge metadata
metadata <- cc_relative %>%
  left_join(cc_area, by = "loci") %>%
  distinct(loci, .keep_all = TRUE)

# Ensure bins are in the same order as columns of the count matrix
bin_order <- colnames(wide_mat)
metadata <- metadata %>% filter(loci %in% bin_order) %>% arrange(match(loci, bin_order))
rownames(metadata) <- metadata$loci

# 7. Create Seurat object ----------------------------------------------------
seurat_obj <- CreateSeuratObject(
  counts = wide_mat,               # genes × bins
  project = sample,
  min.features = 200,              # bins with <200 genes removed
  min.cells = 3,                   # genes detected in <3 bins removed
  meta.data = metadata
)

# 8. Save the Seurat object --------------------------------------------------
saveRDS(seurat_obj, file = out_rds)
cat("Seurat object saved to:", out_rds, "\n")
