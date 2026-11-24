#!/usr/bin/env Rscript
# Sanity checks for the HCs produced by `nextflow run main.nf -profile test`.
#
# Usage: Rscript testdata/check_hcs.R <HCs.csv> [nHCs (default 10)] [plot.png]
#
# Checks: one row per sample in testdata/samples.txt, nHCs columns, all values finite, and
# that the leading HCs separate the 1000G superpopulations (rows are in BCF sample order,
# the same as samples_superpop.txt). Exits non-zero if a check fails.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: check_hcs.R <HCs.csv> [nHCs] [plot.png]")
hcs_file <- args[1]
nHCs <- if (length(args) >= 2) as.integer(args[2]) else 10L
png_file <- if (length(args) >= 3) args[3] else NA

testdata_dir <- local({
  f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(f) == 1) dirname(normalizePath(f)) else "testdata"
})

pops <- read.table(file.path(testdata_dir, "samples_superpop.txt"),
                   col.names = c("sample", "superpop"), stringsAsFactors = FALSE)
hcs <- as.matrix(read.csv(hcs_file))

fail <- function(...) { cat("FAIL:", ..., "\n"); quit(save = "no", status = 1) }
ok <- function(...) cat("ok:", ..., "\n")

# --- shape and values ---
if (nrow(hcs) != nrow(pops)) fail("expected", nrow(pops), "rows (samples), got", nrow(hcs))
ok(nrow(hcs), "rows (one per sample)")
if (ncol(hcs) != nHCs) fail("expected", nHCs, "columns, got", ncol(hcs))
ok(ncol(hcs), "columns")
if (!all(is.finite(hcs))) fail("non-finite values in HCs")
ok("all values finite")

# --- population structure ---
# R^2 of each HC regressed on superpopulation: the leading HCs should be largely explained
# by superpopulation in a 5-superpopulation 1000G sample.
r2 <- sapply(seq_len(min(5, ncol(hcs))), function(i)
  summary(lm(hcs[, i] ~ factor(pops$superpop)))$r.squared)
cat("R^2 of HC1..HC", length(r2), " on superpopulation: ",
    paste(sprintf("%.2f", r2), collapse = " "), "\n", sep = "")
if (max(r2) < 0.5) fail("no leading HC is strongly associated with superpopulation (max R^2 < 0.5)")
ok("leading HCs separate superpopulations")

cat("Mean HC1/HC2 by superpopulation:\n")
print(aggregate(hcs[, 1:2], list(superpop = pops$superpop), mean), digits = 3)

if (!is.na(png_file)) {
  cols <- c(AFR = "#E69F00", AMR = "#CC79A7", EAS = "#009E73", EUR = "#0072B2", SAS = "#D55E00")
  png(png_file, width = 1200, height = 600)
  par(mfrow = c(1, 2))
  for (p in list(c(1, 2), c(3, 4))) {
    plot(hcs[, p[1]], hcs[, p[2]], col = cols[pops$superpop], pch = 19,
         xlab = colnames(hcs)[p[1]], ylab = colnames(hcs)[p[2]])
    legend("topright", legend = names(cols), col = cols, pch = 19)
  }
  dev.off()
  ok("plot written to", png_file)
}
