#!/usr/bin/env Rscript

USAGESTRING <- paste0(
  "calculate_haplotype_components <input file> <sample size (number of individuals)> <no. haplotype components required out> <output file> [row sums file]\n",
  "\n",
  "  If a row sums file (combine's <prefix>.sums: \"index1 sum\", one line per row) is given, each row's\n",
  "  values are rescaled to add up to that row's sum before the SVD. The sums are calculated before\n",
  "  combine --restrict drops entries, so this restores each row's total after restricting.\n"
)

args <- commandArgs(trailingOnly = TRUE)

for(argument in args) {
  if(argument == "-h" || argument == "--help") {
    cat(USAGESTRING)
    quit(save = "no", status = 0)
  }
}

if(!(length(args) %in% c(4, 5))) {
  cat("Error: Wrong number of arguments provided.\n")
  cat(USAGESTRING)
  quit(save = "no", status = 1)
}


infile <- args[1]
nind <- as.numeric(args[2])
nHCs <- as.numeric(args[3])
outfile <- args[4]
sumsfile <- if(length(args) == 5) args[5] else NA


library(data.table)
library(Matrix)
library(sparsesvd)


cat("Reading data\n")
cl <- fread(infile)

if(!is.na(sumsfile)) {
  cat("Reweighting rows to their pre-restrict sums\n")
  sums <- fread(sumsfile, col.names = c("index1", "sum"))
  kept <- cl[, .(kept = sum(V3)), by = V1]
  target <- sums$sum[match(kept$V1, sums$index1)]
  if(anyNA(target)) {
    cat("Error: rows in", infile, "are missing from", sumsfile, "\n")
    quit(save = "no", status = 1)
  }
  # every row in cl has at least one entry > 0, so kept > 0
  kept[, factor := target / kept]
  cl[kept, V3 := V3 * i.factor, on = "V1"]
  cat(sprintf("  %d rows; scale factors range from %.6f to %.6f\n",
              nrow(kept), min(kept$factor), max(kept$factor)))
  nempty <- sum(!(sums$index1 %in% kept$V1) & sums$sum > 0)
  if(nempty > 0) {
    cat(sprintf("  %d rows had all their entries dropped, and stay empty\n", nempty))
  }
}

cat("Making sparse matrix\n")
A <- sparseMatrix(cl$V1, cl$V2, x = log10(cl$V3 + 1), dims = c(nind, nind))

cat("Calculating SVD \n")
res <- sparsesvd(A, rank = nHCs)

cat("Calculating HCs\n")
HCs <- res$u %*% diag(sqrt(res$d))

cat("Writing HCs\n")
colnames(HCs) <- paste0("HC", 1:ncol(HCs))
fwrite(HCs, outfile, row.names = FALSE, sep = ",")
