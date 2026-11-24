#!/usr/bin/env Rscript
# Tests of the HC step (calculate_haplotype_components) on the HC-step test data (testdata/hc),
# without painting. Run from the repository root, in the calculate_HCs image (which has the R
# packages), e.g.:
#
#   docker run --rm -v "$PWD":/repo -w /repo mrcieu/pipeline-nf-haplotype-components:calculate_HCs \
#       Rscript tests/test_hc_step.R
#
# It runs the R script in the repository (dependencies/calculate_HCs/calculate_haplotype_components.R),
# not the copy installed in the image, so it tests the working copy. Exits non-zero if a test fails.

suppressPackageStartupMessages(library(data.table))

script <- "dependencies/calculate_HCs/calculate_haplotype_components.R"
hc <- "testdata/hc"
if (!file.exists(script) || !dir.exists(hc)) stop("run from the repository root")

tmp <- tempfile("hc_step_tests_")
dir.create(tmp)
base <- c("--matrix", file.path(hc, "combined.chunklengths.txt.gz"), "--samples", file.path(hc, "samples.txt"),
          "--sums", file.path(hc, "combined.chunklengths.sums"), "--nhcs", "20")

failures <- 0
check <- function(ok, what) {
  cat(if (isTRUE(ok)) "  ok:  " else "  FAIL:", what, "\n")
  if (!isTRUE(ok)) failures <<- failures + 1
}

# run the script; returns the exit status and the log
run <- function(name, ...) {
  log <- file.path(tmp, paste0(name, ".log"))
  status <- system2("Rscript", c(script, ..., "--out", file.path(tmp, paste0(name, ".csv")),
                                 "--diagnostics", file.path(tmp, paste0(name, ".tsv"))),
                    stdout = log, stderr = log)
  list(status = status, log = readLines(log))
}
hcs <- function(name) fread(file.path(tmp, paste0(name, ".csv")))
diag <- function(name) fread(file.path(tmp, paste0(name, ".tsv")))
maxdiff <- function(a, b) max(abs(as.matrix(a[, -1]) - as.matrix(b[, -1])))

samples <- readLines(file.path(hc, "samples.txt"))
spt <- fread(file.path(hc, "samples_superpop.txt"), header = FALSE, col.names = c("sample", "sp"))
pairs <- fread(file.path(hc, "related_pairs.txt"), header = FALSE)
related <- samples %in% c(pairs$V1, pairs$V2)
writeLines(samples, file.path(tmp, "all.txt"))
writeLines(rev(samples), file.path(tmp, "all_reversed.txt"))
writeLines(samples[!related], file.path(tmp, "unrelated.txt"))
writeLines(c(spt[sp == "AFR", sample], "NOTREAL1"), file.path(tmp, "afr.txt"))
writeLines(c("ID1 ID2 Kinship", paste(pairs$V1, pairs$V2, 0.25), "NOTREAL1 NOTREAL2 0.25"), file.path(tmp, "kinship.txt"))

# HC2-HC10 size of related samples relative to unrelated samples of the same superpopulation (median over superpops)
size_ratio <- function(x) {
  x <- x[match(samples, x$sample)]
  n <- sqrt(rowSums(as.matrix(x[, HC2:HC10])^2))
  median(sapply(split(seq_along(related), spt$sp), function(w) median(n[w][related[w]]) / median(n[w][!related[w]])))
}
superpop_r2 <- function(x, hc) {
  x <- x[match(spt$sample, x$sample)]
  summary(lm(x[[hc]] ~ factor(spt$sp)))$r.squared
}

cat("Default run\n")
r <- run("default", base)
check(r$status == 0, "exit status 0")
d0 <- hcs("default")
check(identical(names(d0), c("sample", paste0("HC", 1:20))), "columns are sample, HC1..HC20")
check(identical(d0$sample, samples), "one row per sample, in sample order")
check(all(is.finite(as.matrix(d0[, -1]))), "all HCs finite")
check(superpop_r2(d0, "HC2") > 0.8, "HC2 separates superpopulations (R^2 > 0.8)")
check(size_ratio(d0) < 0.6, sprintf("relatives' HCs are shrunk without any handling (size ratio %.2f < 0.6)", size_ratio(d0)))
g0 <- diag("default")
check(identical(g0$sample, samples) && all(g0$status == "fitted"), "diagnostics: one row per sample, all fitted")
check(median(g0$largest_share[related]) > 0.5 && median(g0$largest_share[!related]) < 0.1,
      "diagnostics: relatives' rows are dominated by their largest entry")

cat("--reweight with nothing removed\n")
r <- run("reweight", base, "--reweight")
check(r$status == 0 && max(abs(diag("reweight")$reweight_factor - 1)) < 1e-4, "scale factors are 1")
check(maxdiff(hcs("reweight"), d0) < 1e-3, "HCs unchanged")

cat("--fit_samples with every sample\n")
for (f in c("all", "all_reversed")) {
  r <- run(paste0("fit_", f), base, "--fit_samples", file.path(tmp, paste0(f, ".txt")))
  check(r$status == 0 && identical(hcs(paste0("fit_", f)), d0), sprintf("identical to the default (%s)", f))
}

cat("--remove_pairs\n")
r <- run("pairs", base, "--remove_pairs", file.path(hc, "related_pairs.txt"))
check(r$status == 0 && any(grepl("Warning: --remove_pairs without --reweight", r$log)), "warns without --reweight")
r <- run("pairs_rw", base, "--remove_pairs", file.path(hc, "related_pairs.txt"), "--reweight")
check(r$status == 0, "exit status 0 with --reweight")
check(any(grepl(sprintf("%d pairs \\(0 skipped", nrow(pairs)), r$log)), "all pairs found")
g <- diag("pairs_rw")
check(all(g$frac_removed[!related] == 0) && all(g$frac_removed[related] > 0), "only relatives' rows lose entries")
check(size_ratio(hcs("pairs_rw")) > 0.75, sprintf("relatives' HCs are restored (size ratio %.2f > 0.75)", size_ratio(hcs("pairs_rw"))))
r <- run("kinship", base, "--remove_pairs", file.path(tmp, "kinship.txt"), "--reweight")
check(r$status == 0 && any(grepl(sprintf("%d pairs \\(2 skipped", nrow(pairs)), r$log)),
      "kinship-style file: header and unknown pair skipped, extra column ignored")
check(identical(hcs("kinship"), hcs("pairs_rw")), "kinship-style file gives the same HCs")

cat("--fit_samples with a group\n")
r <- run("afr", base, "--fit_samples", file.path(tmp, "afr.txt"), "--reweight")
check(r$status == 0 && any(grepl("1 IDs in .* skipped", r$log)), "unknown ID skipped")
a <- hcs("afr")
check(identical(a$sample, spt[sp == "AFR", sample]), "only the group's samples get HCs, in sample order")
g <- diag("afr")
check(sum(g$status == "fitted") == nrow(a) && all(g$status[!samples %in% a$sample] == "not fitted"), "diagnostics status")

cat("--fit_samples unrelated + --project\n")
r <- run("unrel", base, "--fit_samples", file.path(tmp, "unrelated.txt"), "--reweight")
r <- run("proj", base, "--fit_samples", file.path(tmp, "unrelated.txt"), "--project", "--reweight")
check(r$status == 0, "exit status 0")
p <- hcs("proj")
check(identical(p$sample, samples), "every sample gets HCs, in sample order")
u <- hcs("unrel")
check(maxdiff(u, p[match(u$sample, p$sample)]) < 1e-10, "--project doesn't change the fitted samples' HCs")
g <- diag("proj")
check(all(g$status[related] == "projected") && all(g$status[!related] == "fitted"), "diagnostics: relatives projected")
check(size_ratio(p) > 0.75, sprintf("projected relatives' HCs are not shrunk (size ratio %.2f > 0.75)", size_ratio(p)))
pm <- as.matrix(p[match(samples, p$sample), HC1:HC10])
cent <- apply(pm[!related, ], 2, function(v) tapply(v, spt$sp[!related], mean))
nearest <- apply(pm[related, ], 1, function(h) rownames(cent)[which.min(colSums((t(cent) - h)^2))])
check(mean(nearest == spt$sp[related]) >= 0.95,
      sprintf("projected relatives are nearest their own superpopulation (%d/%d)", sum(nearest == spt$sp[related]), sum(related)))
r <- run("proj_all", base, "--fit_samples", file.path(tmp, "all.txt"), "--project")
check(r$status == 0 && identical(hcs("proj_all"), d0), "--project with every sample fitted is identical to the default")

# exact test of the projection formula: a "clone" sample whose row is a copy of a fitted sample's row, and
# which no fitted sample has an entry to, must be projected to exactly that sample's fitted HCs
m <- fread(file.path(hc, "combined.chunklengths.txt.gz"), header = FALSE, col.names = c("i", "j", "v"))
sm <- fread(file.path(hc, "combined.chunklengths.sums"), header = FALSE, col.names = c("i", "sum"))
orig <- 1L
clone <- length(samples) + 1L
fwrite(rbind(m, m[i == orig][, i := clone]), file.path(tmp, "clone.txt"), sep = " ", col.names = FALSE)
fwrite(rbind(sm, sm[i == orig][, i := clone]), file.path(tmp, "clone.sums"), sep = " ", col.names = FALSE)
writeLines(c(samples, "CLONE"), file.path(tmp, "clone_samples.txt"))
clone_base <- c("--matrix", file.path(tmp, "clone.txt"), "--samples", file.path(tmp, "clone_samples.txt"),
                "--sums", file.path(tmp, "clone.sums"), "--nhcs", "20", "--fit_samples", file.path(tmp, "all.txt"), "--project")
for (rw in c(FALSE, TRUE)) {
  name <- if (rw) "clone_rw" else "clone"
  r <- run(name, clone_base, if (rw) "--reweight")
  x <- hcs(name)
  check(r$status == 0 && maxdiff(x[sample == "CLONE"], x[sample == samples[orig]]) < 1e-8,
        sprintf("a projected copy of a fitted row gets exactly that row's HCs%s", if (rw) " (with --reweight)" else ""))
  ref <- if (rw) hcs("reweight") else d0
  check(maxdiff(x[sample != "CLONE"], ref) < 1e-8, sprintf("the projected clone doesn't change the fitted HCs%s", if (rw) " (with --reweight)" else ""))
}

cat("Errors\n")
check(run("e1", base, "--project")$status != 0, "--project without --fit_samples fails")
check(run("e3", c(base[1:6], "--nhcs", "200"), "--fit_samples", file.path(tmp, "afr.txt"))$status != 0,
      "--nhcs >= number of fitted samples fails")
writeLines(samples[1:5], file.path(tmp, "five.txt"))
check(run("e4", "--matrix", file.path(hc, "combined.chunklengths.txt.gz"), "--samples", file.path(tmp, "five.txt"),
          "--nhcs", "2")$status != 0, "matrix indices beyond the sample list fail")
check(run("e5", "--matrix", file.path(hc, "combined.chunklengths.txt.gz"), "--samples", file.path(hc, "samples.txt"),
          "--nhcs", "2", "--reweight")$status != 0, "--reweight without --sums fails")

unlink(tmp, recursive = TRUE)
cat(if (failures == 0) "\nAll HC-step tests passed\n" else sprintf("\n%d HC-step test(s) FAILED\n", failures))
quit(save = "no", status = if (failures == 0) 0 else 1)
