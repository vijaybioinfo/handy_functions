#!/usr/bin/R

###################
# Seurat plotting #
###################

# This functions are designed to be help with visualisation of Seurat object

qc_violin <- function(
  mydat,
  xax = "Data",
  yax = colnames(mydat)[1],
  lb_filt = NULL,
  hb_filt = NULL,
  filtby = NULL
){
  if(!"Data" %in% colnames(mydat)) mydat$Data <- "SET"
  if(!is.factor(mydat[, xax])){
    mydat[, xax] <- factor(mydat[, xax])
  }
  mydat <- mydat[-which.max(mydat[, yax]), ] # when the violin is squashed

  if(!any(is.null(c(lb_filt, hb_filt)))){
    mydat$filter <- mydat[, yax] >= lb_filt[yax] & mydat[, yax] <= hb_filt[yax]
    pcts <- paste0(table(mydat[mydat$filter, xax]), "/", table(mydat[, xax]))
    names(pcts) <- levels(mydat[, xax])
    mydat$pct <- pcts[as.character(mydat[, xax])]
  }
  p <- simple_violin(dat = mydat, x = xax, y = yax) +
    theme(
      axis.text.x = element_text(hjust = 1, angle = 45),
      legend.position = "none",
      axis.title.x = element_blank()
    )
  if('pct' %in% colnames(mydat)){
    ymax <- max(mydat[, yax], na.rm = TRUE)
    ll <- c(lb_filt[yax], min(mydat[, yax], na.rm = TRUE))
    hl <- c(hb_filt[yax], ymax); coly <- ifelse(isTRUE(filtby[yax]), 'red', '#808080')
    p <- p + geom_hline(yintercept = max(ll), linetype = "dashed", color = coly) +
      geom_hline(yintercept = min(hl), linetype = "dashed", color = coly) +
      geom_text(
        aes_string(label = 'pct', y = ymax), # position = position_dodge(0.9),
        nudge_x = ifelse(nlevels(mydat[, xax]) < 3, 0.3, 0),
        vjust = 1.6, angle = 45
      )
  }
}

# plot two variables with a rectangle of cutoffs
GenePlotc <- function(
  object,
  x,
  do.return = FALSE
){
  tvar <- colnames(x)
  void <- FeatureScatter(object = object, feature1 = tvar[1], feature2 = tvar[2])
  if(diff(range(object@meta.data[, tvar[1]])) > maxdiff) void <- void + scale_x_log10()
  if(diff(range(object@meta.data[, tvar[2]])) > maxdiff) void <- void + scale_y_log10()
  for(i in seq(1, nrow(x), 2)){
    # if(vrs < 3){ rect(x[i, 1], x[i, 2], x[i + 1, 1], x[i + 1, 2], border = ifelse(i + 1 == nrow(x), 'red', 'black')); next }
    void <- void + annotate("rect", xmin = x[i, 1], ymin = x[i, 2], xmax = x[i + 1, 1], ymax = x[i + 1, 2],
      alpha = 0.1, color = ifelse(i + 1 == nrow(x), 'red', 'black'), fill = 'white')
  }
  if(isTRUE(do.return)) return(void)
  print(void)
}


#' Scatter plot of features and their thresholds
#'
#' This function creates a scatter plot of two (third optional) features or
#' variables with their respective thresholds
#'
#' @param mat matrix or data.frame of numeric values, additionaly you can give
#' a Seurat object.
#' @param object Annotation (data.frame) or object containing the variables
#' @param variables Vector of variables to plot.
#' @param thresh Vector/table/data.frame with 'variables' columns and two values of
#' @param simple Use a simple thresholding (usually when you have only one
#' threshold to show).
#' @param verbose Show progress.
#' @keywords ROC
#'
#' @return Returns a data.frame with the signature (generates a lot of reports)
#'
#' @importFrom pROC pheatmap
#'
#' @export
#'
#' @examples
# ' addinfo <- get_signature(mat = edata, annot = annot_group, genes = c('IL5', 'IL4', 'IL9'))
#'

feature_scatter <- function(
  object,
  variables,
  thresh = NULL,
  verbose = FALSE
){
  if(casefold(class(object)) == "seurat"){ object <- FetchData(object, vars = variables); gc() }
  p <- plot_corr(
    df = object, var1 = variables[1],  var2 = variables[2], addvar = variables[3],
    add_line = FALSE, log2t = 'auto', return_plot = TRUE, v = verbose
  )
  p <- p + theme_minimal();
  if(!is.null(thresh)){ # when no threshold is given
    if(verbose) cat("Adding boundaries\n")
    for(i in seq(1, nrow(thresh), 2)){
      p <- p + annotate(geom = "rect",
        xmin = thresh[i, 1], ymin = thresh[i, 2], xmax = thresh[i + 1, 1], ymax = thresh[i + 1, 2],
        alpha = 0.1, color = ifelse(i + 1 == nrow(thresh), 'red', 'black'), fill = 'white'
      )
    }
  }
  p; #dev.off()
}

ncounts_outsies <- function(
  object,
  redu = c("PC_1", "PC_2", "nFeature_RNA", "nCount_RNA" , "percent.mt")
){
  tvar <- if(casefold(class(object)) == "seurat") FetchData(object = object, vars = redu) else object
  tvar$cell_name <- rownames(tvar); x <- order(abs(tvar[, 1]))
  y <- order(abs(tvar[, 2])); outsies <- tvar[unique(c(tail(x), tail(y))), ]; rm(x, y)
  void <- ggplot(tvar, aes_string(redu[1], redu[2])) + labs(x = redu[1], y = redu[2]) + theme_classic()
  nplots <- lapply(redu[-c(1:2)], function(x){
    void + geom_point(aes_string(colour = x), size = 0.7) + scale_colour_gradientn(colours = rainbow(7)) +
      ggrepel::geom_text_repel(data = outsies, aes(label = cell_name), size = 0.75)
  })
  print(cowplot::plot_grid(plotlist = nplots))
  return(outsies)
}

# mean-variance from HVFInfo
VariableFeaturePlotp <- function(
  object,
  cols = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000'),
  pt.size = 1.5,
  log = NULL,
  axes = NULL,
  selection.method = NULL,
  assay = NULL
) {
  library(cowplot)
  hvf.info <- HVFInfo(object = object, assay = assay, selection.method = selection.method, status = TRUE)
  var.status <- c("no", "yes")[unlist(x = hvf.info[, ncol(x = hvf.info)]) + 1]
  hvf.info <- hvf.info[, c(1, 3)]
  hvf.info$pct <- round(Matrix::rowMeans(GetAssayData(object, assay = assay)[rownames(hvf.info), ] > 0) * 100, 1)
  # 1-5, 5-10, 10-20, 2-50, 50-100 # try to bin the data
  # hvf.info$bpct <- bin(hvf.info$pct, n = length(cols))
  # hvf.info$bpct <- switch(hvf.info$pct, n = length(cols))
  if(is.null(axes)){
    log <- !is.null(log) || (any(c("variance.standardized", "residual_variance") %in% colnames(x = hvf.info)))
    axes <- c('mean', 'variance.standardized', 'pct')
  }else{ log <- FALSE }
  axis.labels <- sapply(axes, switch,
    variance.standardized = "Standardized Variance",
    mean = "Average Expression",
    dispersion.scaled = "Dispersion",
    gmean = "Geometric Mean of Expression",
    pct = "Percentage expressing",
    "log10(mean)" = "Log10(Mean)",
    residual_variance = "Residual Variance")
  if(is.null(axis.labels[1])) axis.labels[1] <- axes[1]
  axis.labels[3] <- gsub(" ", "\n", axis.labels[3])
  hvf.info$tmp <- hvf.info[, axes[3]]
  hvf.info[var.status == "no", "tmp"] <- NA
  if(grepl("pct", axes[3])) lms <- c(0, 100) else lms <- range(hvf.info[, axes[3]])
  plot <- ggplot(data = hvf.info, aes_string(x = axes[1], y = axes[2], color = "tmp")) +
    geom_point(size = pt.size) +
    scale_color_gradientn(name = axis.labels[3], colours = cols, na.value = 'gray70', limits = lms, guide = "colorbar") +
    theme_cowplot() + theme(plot.title = element_text(hjust = 0.5)) +
    labs(subtitle = paste(c("Non-variable", "Variable"), "count:", table(var.status), collapse = "\n"),
      x = axis.labels[1], y = axis.labels[2])
  if (log) {
      plot <- plot + scale_x_log10()
  }
  return(plot)
}

## Visualise fit
mean_var <- function(object){
  hvf.info <- HVFInfo(object, selection.method = "vst")
  not.const <- hvf.info$variance > 0
  fit <- loess(
     formula = log10(x = variance) ~ log10(x = mean),
     data = hvf.info[not.const, ],
     span = 0.3
  )
  hvf.info$variance.expected[not.const] <- 10 ^ fit$fitted
  not.const <- sample(which(not.const), min(c(100, sum(not.const) * 0.01)))
  p <- ggplot(hvf.info, aes(log10(mean), log10(variance))) + geom_point(alpha = 0.6) +
    geom_smooth(se = FALSE, color = 'red', method = "loess") +
    geom_point(data = hvf.info[not.const, ], aes(log10(mean), log10(variance.expected)), alpha = 0.6, size = 0.1, color = "green") +
    labs(title = 'Global mean-variance LOESS')
  return(p )
}
