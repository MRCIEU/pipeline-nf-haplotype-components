#!/usr/bin/env Rscript

# Computes haplotype components (HCs) from the combined chunk length matrix written by combine.
#
# The matrix (combine's <prefix>.txt.gz) has three columns, "index1 index2 value", where the indices
# are 1-based positions in the samples file (the sample order of the BCF/VCFs that were painted).
# HCs are the top --nhcs columns of U * sqrt(D), from the SVD of log10(matrix + 1).

suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option("--matrix", type = "character",
              help = "combined chunk length matrix from combine (<prefix>.txt.gz) [required]"),
  make_option("--samples", type = "character",
              help = "sample IDs, one per line, in the order of the painted BCF/VCFs [required]"),
  make_option("--nhcs", type = "integer",
              help = "number of HCs to compute [required]"),
  make_option("--out", type = "character",
              help = "output file for the HCs (CSV: sample, HC1, HC2, ...) [required]"),
  make_option("--sums", type = "character", default = NA,
              help = paste("combine's row sums file (<prefix>.sums: \"index1 sum\"), calculated before",
                           "combine --restrict dropped entries. Needed for --reweight; also used in the diagnostics")),
  make_option("--remove_pairs", type = "character", default = NA,
              help = paste("optional file of sample ID pairs (e.g. relatives) whose entries are removed from the matrix,",
                           "in both directions: the first two whitespace-separated columns, other columns ignored.",
                           "Pairs with an ID not in --samples (or a header line) are skipped")),
  make_option("--reweight", action = "store_true", default = FALSE,
              help = paste("rescale each row's values to add up to its sum in --sums (its total before",
                           "combine --restrict and --remove_pairs removed entries) before the SVD")),
  make_option("--diagnostics", type = "character", default = NA,
              help = "optional output file for per-sample diagnostics (tab-separated)")
)
parser <- OptionParser(
  usage = "calculate_haplotype_components --matrix FILE --samples FILE --nhcs N --out FILE [options]",
  option_list = option_list
)
opt <- parse_args(parser)

fail <- function(...) {
  cat("Error:", ..., "\n", file = stderr())
  quit(save = "no", status = 1)
}

for (required in c("matrix", "samples", "nhcs", "out")) {
  if (is.null(opt[[required]])) {
    print_help(parser)
    fail(sprintf("--%s is required", required))
  }
}
if (opt$reweight && is.na(opt$sums)) fail("--reweight needs --sums")
if (!is.na(opt$remove_pairs) && !opt$reweight) {
  cat(paste("Warning: --remove_pairs without --reweight. Rows that lose entries (e.g. a relative's) keep only",
            "the rest of their total, which can shrink their HCs towards 0; --reweight restores their totals.\n"))
}

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(sparsesvd)
})

# --- inputs ---
cat("Reading samples\n")
samples <- readLines(opt$samples)
samples <- samples[nzchar(samples)]
nind <- length(samples)
if (nind == 0) fail("no samples in", opt$samples)
if (anyDuplicated(samples)) fail("duplicate sample IDs in", opt$samples)
if (opt$nhcs >= nind) fail(sprintf("--nhcs (%d) must be less than the number of samples (%d)", opt$nhcs, nind))
cat(sprintf("  %d samples\n", nind))

cat("Reading matrix\n")
cl <- fread(opt$matrix, header = FALSE, col.names = c("i", "j", "v"))
cat(sprintf("  %d entries\n", nrow(cl)))
if (nrow(cl) > 0 && max(cl$i, cl$j) > nind) {
  fail(sprintf("the matrix has sample index %d, but there are only %d samples in %s",
               max(cl$i, cl$j), nind, opt$samples))
}

# per-row totals before the cap (from combine), and of the entries that are left
rows <- data.table(i = seq_len(nind), sample = samples)
if (!is.na(opt$sums)) {
  sums <- fread(opt$sums, header = FALSE, col.names = c("i", "total_precap"))
  if (nrow(sums) > 0 && max(sums$i) > nind) fail("the sums file has rows that are not in", opt$samples)
  rows[sums, total_precap := i.total_precap, on = "i"]
} else {
  rows[, total_precap := NA_real_]
}
kept <- cl[, .(n_entries = .N, total_kept = sum(v), largest = max(v)), by = i]
rows[kept, `:=`(n_entries = i.n_entries, total_kept = i.total_kept, largest = i.largest), on = "i"]
rows[is.na(n_entries), `:=`(n_entries = 0L, total_kept = 0, largest = NA_real_)]

# --- removing pairs (e.g. relatives) ---
rows[, total_removed := 0]
if (!is.na(opt$remove_pairs)) {
  cat("Removing pairs\n")
  pairs <- fread(opt$remove_pairs, header = FALSE, select = 1:2, col.names = c("id1", "id2"), colClasses = "character")
  pairs[, `:=`(i = match(id1, samples), j = match(id2, samples))]
  nskip <- pairs[is.na(i) | is.na(j), .N]
  pairs <- unique(pairs[!is.na(i) & !is.na(j) & i != j])
  both <- unique(rbind(pairs[, .(i, j)], pairs[, .(i = j, j = i)]))
  removed <- cl[both, on = c("i", "j"), nomatch = 0]
  cl <- cl[!both, on = c("i", "j")]
  rows[removed[, .(r = sum(v)), by = i], total_removed := i.r, on = "i"]
  cat(sprintf("  %d pairs (%d skipped: an ID is not in %s); removed %d entries from %d rows\n",
              nrow(pairs), nskip, opt$samples, nrow(removed), uniqueN(removed$i)))
}
used <- cl[, .(n_used = .N, total_used = sum(v)), by = i]
rows[used, `:=`(n_entries_used = i.n_used, total_used = i.total_used), on = "i"]
rows[is.na(n_entries_used), `:=`(n_entries_used = 0L, total_used = 0)]

# --- reweighting ---
rows[, reweight_factor := 1]
if (opt$reweight) {
  cat("Reweighting rows to their sums before the cap and removals\n")
  if (rows[n_entries_used > 0, anyNA(total_precap)]) fail("rows in the matrix are missing from", opt$sums)
  # every row with entries has values > 0, so total_used > 0
  rows[n_entries_used > 0, reweight_factor := total_precap / total_used]
  cl[rows, v := v * i.reweight_factor, on = "i"]
  cat(sprintf("  scale factors range from %.6f to %.6f\n",
              rows[n_entries_used > 0, min(reweight_factor)], rows[n_entries_used > 0, max(reweight_factor)]))
}

nempty <- rows[n_entries_used == 0, .N]
if (nempty > 0) cat(sprintf("  %d samples have no entries in the matrix; their HCs will be 0\n", nempty))

# --- SVD ---
cat("Making sparse matrix\n")
A <- sparseMatrix(cl$i, cl$j, x = log10(cl$v + 1), dims = c(nind, nind))
rm(cl); invisible(gc())

cat("Calculating SVD\n")
res <- sparsesvd(A, rank = opt$nhcs)
if (length(res$d) < opt$nhcs) {
  cat(sprintf("  Warning: only %d non-zero singular values; writing %d HCs\n", length(res$d), length(res$d)))
}

cat("Calculating HCs\n")
HCs <- res$u %*% diag(sqrt(res$d), nrow = length(res$d))
colnames(HCs) <- paste0("HC", seq_len(ncol(HCs)))

cat("Writing HCs\n")
fwrite(data.table(sample = samples, HCs), opt$out, sep = ",")

# --- diagnostics ---
if (!is.na(opt$diagnostics)) {
  cat("Writing diagnostics\n")
  rows[, `:=`(
    frac_kept = total_kept / total_precap,      # share of the row's pre-cap total left after the cap
    largest_share = largest / total_kept,       # share of the row (after the cap) taken by its largest entry
    frac_removed = total_removed / total_kept,  # share of the row (after the cap) removed by --remove_pairs
    status = "fitted"
  )]
  fwrite(rows[, .(sample, n_entries, total_precap, total_kept, frac_kept, largest_share,
                  frac_removed, n_entries_used, reweight_factor, status)],
         opt$diagnostics, sep = "\t", na = "NA", quote = FALSE)
}
