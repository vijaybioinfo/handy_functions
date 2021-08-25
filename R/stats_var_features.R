#!/usr/bin/R

#' Most/Highly variable features
#'
#' This function calculates the most variables features from a matrix and plots
#' their distribution in a coef. of variation by mean scatter.
#'
#' @param counts Matrix of counts.
#' @param normalise Normalise with DESeq2.
#' @param lib.size Library size vector matching the columns.
#' @param padjthr Adjusted P-values threshold of significance.
#' @param fitv Minimum variance to derive a fit line from. May need to increase
#' if fit line drops dramatically.
#' @param top_n Top variable features.
#' @param take_sig Take only the significant ones (still tries the 'top_n').
#' @param plot Create the scatter of mean vs sq. coef. of variation.
#' @param verbose Show steps.
#' @references \url{http://pklab.med.harvard.edu/scw2014/subpop_tutorial.html}
#' @keywords Variablility
#' @return Returns Object with highly variable features
#'
#' @importFrom DESeq2 estimateSizeFactorsForMatrix
#' @importFrom stats var
#' @importFrom statmod glmgam.fit
#'
#' @export
#'
#' @examples
# ' addinfo <- getMostVariableGenes(counts = cts, plot = TRUE)
#'

getMostVariableGenes <- function(
  counts,
  normalise = TRUE,
  lib.size = NULL,
  padjthr = 0.05,
  fitv = 0.3,
  top_n = 500,
  take_sig = TRUE,
  plot = FALSE,
  verbose = FALSE
){
  # suppressPackageStartupMessages(library(DESeq2))
  counts <- as.matrix(counts)
  counts <- counts[base::rowMeans(counts) > 0, ]

  if(normalise){
    if(verbose) cat("Normalising\n")
    if(is.null(lib.size)) lib.size <- try(DESeq2::estimateSizeFactorsForMatrix(counts), silent = TRUE)
    if(class(lib.size)[1] == 'try-error'){
      lib.size <- DESeq2::estimateSizeFactorsForMatrix(
        counts,
        geoMeans = apply(counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
      )
    }
    ed <- t(t(counts)/lib.size)
  }else{
    ed <- counts
  }; rm(counts)

  if(verbose) cat("Calculating mean and coefficient of variation\n")
  means <- base::rowMeans(ed)
  vars <- apply(ed, 1, stats::var)
  cv2 <- vars / means^2 # sq. coef. of variation

  if(verbose) cat("Filter mean for fitting: >= .95 Q of variance >", fitv, "\n")
  minMeanForFit <- unname( quantile( means[ which( cv2 > fitv ) ], .95 ) )
  if(verbose) cat("Minimum mean detected", minMeanForFit, "\n")
  useForFit <- means >= minMeanForFit # & spikeins
  fit <- statmod::glmgam.fit(cbind(a0 = 1, a1tilde = 1/means[useForFit]), cv2[useForFit])
  a0 <- unname( fit$coefficients["a0"] )
  a1 <- unname( fit$coefficients["a1tilde"])

  afit <- a1/means+a0
  varFitRatio <- vars/(afit*means^2)
  varorder <- rownames(ed[order(varFitRatio, decreasing = TRUE), ])

  df <- ncol(ed) - 1
  if(verbose) cat("Chi test and its adjusted P-values\n")
  pval <- pchisq(varFitRatio * df, df = df, lower.tail = FALSE)
  adj.pval <- p.adjust(pval, "fdr")
  sigVariedGenes <- adj.pval < padjthr;
  if(verbose) cat("Passed ( <", padjthr, "):", sum(sigVariedGenes), "\n")

  if(isTRUE(take_sig)){
    if(verbose) cat("Taking only significative features\n")
    varorder <- varorder[varorder %in% names(adj.pval[sigVariedGenes])]
  }
  varorder <- head(varorder, top_n)
  if(verbose) cat("Found top", length(varorder), "\n")

  xg <- exp(seq(min(log(means[means>0])), max(log(means)) , length.out = 1000))
  vfit <- a1/xg + a0

  to_return <- list(
    means = means, cv2 = cv2, varFitRatio = varFitRatio, useForFit = useForFit,
    vfit = vfit, df = df, xg = xg, top_n = names(means) %in% varorder,
    pval.adj = adj.pval, is.signif = sigVariedGenes
  )
  if(plot) plot_hvg(to_return)
  return(to_return)
}

plot_hvg <- function(
  ddf, # list returned by getMostVariableGenes
  show_point = c('top_n', 'is.signif') # which points to show: column, is.signif or top_n
){
  show_point <- match.arg(show_point)
  show_point <- ddf[[show_point]]
  df <- unique(ddf$df)
  par(mar=c(3.5,3.5,1,1),mgp=c(2,0.65,0),cex=0.9);
  smoothScatter(
    log(ddf$means), log(ddf$cv2),
    xlab = paste("Logged Mean, fit-min:", round(min(ddf$cv2[ddf$useForFit]), 3)),
    ylab = paste("Logged Variance-to-mean ratio, fit-min:", round(min(ddf$means[ddf$useForFit]), 3))
  );
  lines(log(ddf$xg), log(ddf$vfit), col="black", lwd=3 );
  lines(log(ddf$xg), log(ddf$vfit * qchisq(0.975, df) / df), lty = 2, col = "black");
  lines(log(ddf$xg), log(ddf$vfit * qchisq(0.025, df) / df), lty = 2, col = "black");
  points(log(ddf$means[show_point]), log(ddf$cv2[show_point]), col = 2)
}
