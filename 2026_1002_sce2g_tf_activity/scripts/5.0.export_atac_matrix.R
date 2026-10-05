#!/usr/bin/env Rscript
# 5.0 Export scE2G's Kendall peak x cell ATAC matrix (dgCMatrix .rds) as raw CSC vectors that
# numpy reads directly: i (0-based row), p (column pointers), x (values), plus row / column names.
# Usage: Rscript 5.0.export_atac_matrix.R <atac_matrix.rds> <out_prefix>
suppressPackageStartupMessages(library(Matrix))
args <- commandArgs(trailingOnly = TRUE)
m <- readRDS(args[1])
cat("class", class(m), "dim", dim(m), "nnz", length(m@x), "\n")
if (!is(m, "dgCMatrix")) m <- as(m, "CsparseMatrix")
cat("value table (first 1e6 nonzeros):\n")
print(table(head(m@x, 1e6)))
cat("rownames:", head(rownames(m), 3), "\ncolnames:", head(colnames(m), 3), "\n")
out <- args[2]
writeBin(as.integer(m@i), paste0(out, ".i.int32"), size = 4)
writeBin(as.integer(m@p), paste0(out, ".p.int32"), size = 4)
writeBin(as.numeric(m@x), paste0(out, ".x.float64"), size = 8)
writeLines(rownames(m), paste0(out, ".rows.txt"))
writeLines(colnames(m), paste0(out, ".cols.txt"))
writeLines(as.character(dim(m)), paste0(out, ".dim.txt"))
cat("export complete\n")
