#!/usr/bin/env Rscript

USAGESTRING <- "calculate_haplotype_components <input file> <sample size (number of individuals)> <no. haplotype components required out> <output file>\n"

args <- commandArgs(trailingOnly = TRUE)

for(argument in args) {
  if(argument == "-h" || argument == "--help") {
    cat(USAGESTRING)
    quit(save = "no", status = 0)
  }
}

if(length(args) != 4) {
  cat("Error: Wrong number of arguments provided.\n")
  cat(USAGESTRING)
  quit(save = "no", status = 1)
}


infile <- args[1]
nind <- as.numeric(args[2])
nHCs <- as.numeric(args[3])
outfile <- args[4]


library(data.table)
library(Matrix)
library(sparsesvd)


cat("Reading data\n")
cl <- fread(infile)

cat("Making sparse matrix\n")
A <- sparseMatrix(cl$V1, cl$V2, x = log10(cl$V3 + 1), dims = c(nind, nind))

cat("Calculating SVD \n")
res <- sparsesvd(A, rank = nHCs)

cat("Calculating HCs\n")
HCs <- res$u %*% diag(sqrt(res$d))

cat("Writing HCs\n")
colnames(HCs) <- paste0("HC", 1:ncol(HCs))
fwrite(HCs, outfile, row.names = FALSE, sep = ",")
