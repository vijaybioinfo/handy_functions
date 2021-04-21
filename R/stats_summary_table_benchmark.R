#!/usr/bin/R

setwd("~/ad_hoc/tests")
library(microbenchmark)
library(ggplot2)
library(Matrix)
library(matrixStats)
library(Biobase)
library(Rfast)
n = 1000; mat = matrix(rnorm(n*n), ncol = n)
groups = rep(LETTERS[1:3], each = n / 3, length.out = n)
# no alternative for matrixStats::rowSds for now
res <- microbenchmark(
  base::rowSums(mat > 0, na.rm = TRUE) / ncol(mat) * 100,
  Matrix::rowSums(mat > 0, na.rm = TRUE) / ncol(mat) * 100,
  Matrix::rowMeans(mat > 0, na.rm = TRUE) * 100,
  matrixStats::rowMeans2(mat > 0, na.rm = TRUE) * 100, # faster
  Rfast::rowmeans(mat > 0) * 100
)
res <- microbenchmark(
  Matrix::rowMeans(mat, na.rm = TRUE), # faster in big numbers
  matrixStats::rowMeans2(mat, na.rm = TRUE), # faster in small numbers
  Rfast::rowmeans(mat) # undoublty faster!!
)
res <- microbenchmark(
  Biobase::rowMedians(mat, na.rm = TRUE),
  matrixStats::rowMedians(mat, na.rm = TRUE),
  Rfast::rowMedians(mat, na.rm = TRUE) # undoublty faster
)
res <- microbenchmark(
  Matrix::mean(c(mat), na.rm = TRUE),
  matrixStats::mean2(c(mat), na.rm = TRUE) # faster
)
# idxs_i = as.numeric(factor(groups)) # idxs only works for a SINGLE group/subset
# matrixStats_b = apply(mat, 1, function(vec) matrixStats::mean2(x = vec, idxs = idxs_i, na.rm = TRUE) )
res <- microbenchmark(
  Matrix_a = apply(mat, 1, function(vec) tapply(vec, groups, Matrix::mean, na.rm = TRUE) ),
  matrixStats_a = apply(mat, 1, function(vec) tapply(vec, groups, matrixStats::mean2, na.rm = TRUE) )
) # matrixStats wins
# head(t(Matrix_a))
# head(t(matrixStats_a))
# all.equal(Matrix_a, matrixStats_a) # TRUE
res <- microbenchmark(
  stats_a = apply(mat, 1, function(vec) tapply(vec, groups, stats::sd, na.rm = TRUE) ),
  BiocGenerics_a = apply(mat, 1, function(vec) tapply(vec, groups, BiocGenerics::sd, na.rm = TRUE) )
) # stats_a looks better in big numbers... and it wins
# Maybe check: https://stackoverflow.com/questions/52459711/how-to-find-cumulative-variance-or-standard-deviation-in-r
res
autoplot(res); dev.off()
