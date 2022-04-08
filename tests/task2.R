setwd("/home/bx2ur/Documents/Courses/BIMS8601/interval_basics")

source("src/intervals.R")

filesList = c(
  "data/helas3_ctcf.narrowPeak.gz",
  "data/helas3_jun.narrowPeak.gz",
  "data/hepg2_ctcf.narrowPeak.gz",
  "data/hepg2_jun.narrowPeak.gz"
)

labels = c("helas3_ctcf", "helas3_jun", "hepg2_ctcf", "hepg2_jun")
files = lapply(filesList, read.table)
granges = lapply(files, function(x) {
  GRanges(seqnames=x$V1, ranges=IRanges(x$V2, x$V3))
})

pwj = pairwiseJaccard(granges)

colnames(pwj) = labels
rownames(pwj) = labels

pwj