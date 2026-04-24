#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tximport)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key), call. = FALSE)
    }
    if (i == length(args)) {
      stop(sprintf("Missing value for %s", key), call. = FALSE)
    }
    out[[substring(key, 3)]] <- args[[i + 1]]
    i <- i + 2
  }
  out
}

opt <- parse_args(args)
required <- c("loci_csv", "tx2gene", "rsem_dir", "out_rds")
missing <- required[!required %in% names(opt)]
if (length(missing) > 0) {
  stop(sprintf("Missing required args: %s", paste(missing, collapse = ", ")), call. = FALSE)
}

loci_df <- read.csv(opt$loci_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (!"locus" %in% colnames(loci_df)) {
  stop(sprintf("%s must contain a locus column", opt$loci_csv), call. = FALSE)
}
loci <- unique(as.character(loci_df$locus[!is.na(loci_df$locus)]))
if (length(loci) == 0) {
  stop(sprintf("No loci found in %s", opt$loci_csv), call. = FALSE)
}

tx2gene <- read.delim(
  opt$tx2gene,
  header = FALSE,
  stringsAsFactors = FALSE,
  col.names = c("TXNAME", "GENEID")
)
if (nrow(tx2gene) == 0) {
  stop(sprintf("%s is empty", opt$tx2gene), call. = FALSE)
}

sample_dirs <- list.dirs(opt$rsem_dir, full.names = TRUE, recursive = FALSE)
sample_dirs <- sort(sample_dirs[file.info(sample_dirs)$isdir])
files <- character()
for (sample_dir in sample_dirs) {
  sample <- basename(sample_dir)
  isoform_results <- file.path(sample_dir, sprintf("%s.isoforms.results", sample))
  if (file.exists(isoform_results)) {
    files[[sample]] <- isoform_results
  }
}

if (length(files) == 0) {
  stop(sprintf("No *.isoforms.results files found under %s", opt$rsem_dir), call. = FALSE)
}

txi <- tximport(
  files = files,
  type = "rsem",
  tx2gene = tx2gene,
  ignoreTxVersion = FALSE,
  countsFromAbundance = "no"
)

reindex_matrix <- function(mat, fill_value) {
  if (is.null(mat)) {
    return(NULL)
  }
  out <- matrix(
    fill_value,
    nrow = length(loci),
    ncol = ncol(mat),
    dimnames = list(loci, colnames(mat))
  )
  common <- intersect(rownames(mat), loci)
  if (length(common) > 0) {
    out[common, ] <- as.matrix(mat[common, , drop = FALSE])
  }
  out
}

txi$counts <- reindex_matrix(txi$counts, 0)
txi$abundance <- reindex_matrix(txi$abundance, 0)
txi$length <- reindex_matrix(txi$length, NA_real_)

out_dir <- dirname(opt$out_rds)
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
}

saveRDS(
  list(
    txi = txi,
    samples = names(files),
    files = files,
    loci_csv = opt$loci_csv,
    tx2gene = opt$tx2gene
  ),
  file = opt$out_rds
)
