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
        vjust = 1.6, angle = 45,
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
  if(log_or_not(object@meta.data[, tvar[1]])) void <- void + scale_x_log10()
  if(log_or_not(object@meta.data[, tvar[2]])) void <- void + scale_y_log10()
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
  thresh = NULL, #
  simple = TRUE,
  verbose = FALSE
){
  if(casefold(class(object)) == "seurat"){ object <- FetchData(object, vars = variables); gc() }
  if(!isTRUE(simple) && (length(thresh) == length(variables))){
    tvar <- thresh[3, head(variables, 3)]
    tmp <- apply(tvar, 2, function(x) all(is.infinite(x)) )
    tvar <- tvar[tmp, apply(tvar, 2, function(x) all(is.infinite(x)) )]
    p <- get_densities(
      mat = t(object[, variables]),
      genes = head(variables, 3),
      log2t = TRUE,
      cuof = tvar,
      pdist = 0, # distance for percentages
      even_axis = FALSE,
      ptsize = 1,
      return_plot = TRUE,
      metas = object,
      v = verbose
    )$scatter
  }else{
    p <- plot_corr(
      df = object, var1 = variables[1],  var2 = variables[2], addvar = variables[3],
      log2t = 'auto', return_plot = TRUE, v = verbose
    )
  }
  p <- p + theme_minimal();
  if(!is.null(thresh)){ # when no threshold is given
    if(verbose) cat("Adding boundaries\n")
    # tvar <- sapply(object[, variables[1:2]], function(x) log_or_not(x) )
    # for(coly in names(tvar[tvar])){ # log it if necessary
    #   thresh[, coly] <- ifelse(is.finite(thresh[, coly]) | thresh[, coly] > 0, log10(thresh[, coly]))
    # }
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

## barplot
plot_pct <- function(
  x,
  groups = c(1, 2), # proportion of 1 per 2
  orderby = FALSE, # group in groups[1] to order by
  normalise = FALSE, # think carefully which variables need it
  print_ptables = FALSE,
  type = c("bar", "pie", "donut"),
  return_table = FALSE,
  v = FALSE
){
  require(dplyr)
  type <- match.arg(type)
  cats <- table(x[, groups[1]])
  ddf <- ddf2 <- table(x[, groups])
  if(v) cat("Type of plot:", type, "\n")
  if(v) cat("Proportions of:", groups[1], "\n")
  if(v) cat("Per:", groups[2], "\n")
  if(isTRUE(print_ptables)) print(ddf2)
  if(isTRUE(normalise)) normalise <- 1
  if(is.numeric(normalise)){
    if(v) cat("Normalising by", groups[normalise], "\n")
    ddf2 <- (ddf2 * (min(table(x[, groups[normalise]])) / rowSums(ddf2))) # factorising/normalising
  }
  if(isTRUE(print_ptables)) print(ddf2)
  prop_table <- prop.table(ddf2, 2); if(isTRUE(print_ptables)) print(prop_table)
  mylevels <- if(is.numeric(orderby) || all(orderby %in% names(cats))){
    direct = FALSE
    orderby <- if(is.numeric(orderby)){
      direct = any(orderby > 0)
      names(cats[abs(orderby)])
    }else{
      orderby
    }
    if(v) cat("Sorting by", orderby, "\n")
    tvar <- if(length(orderby) > 1) colSums(prop_table[orderby, ]) else prop_table[orderby, ]
    names(sort(tvar, decreasing = direct))
  }else{ NULL }
  prop_table <- round(as.data.frame.matrix(rbind(ddf[, colnames(prop_table)], prop_table)), 2)
  prop_table <- cbind(
    Group = rep(rownames(ddf2), times = 2),
    Type = rep(c("Count", "Percentage"), each = nrow(ddf2)),
    prop_table
  )
  rownames(prop_table) <- NULL

  ddf2 <- reshape2::melt(ddf2)
  ddf2$fill_group <- if(!is.factor(x[, groups[1]])) factormix(ddf2[, 1]) else factor(ddf2[, 1], levels(x[, groups[1]]))
  ddf2$in_group <- if(!is.null(mylevels)){
    factor(ddf2[, 2], mylevels)
  }else if(!is.factor(x[, groups[2]])){
    factormix(ddf2[, 2])
  }else{
    factor(ddf2[, 2], levels(x[, groups[2]]))
  }
  ddf2 <- ddf2[, -c(1:2)]
  if(type == "bar"){
    p <- ggplot(data = ddf2, aes(x = in_group, y = value, fill = fill_group)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(fill = make_title(groups[1]), y = "Percentage")
  }else if(type == "pie"){
    p <- ggplot(data = ddf2, aes(x = "", y = value)) +
    geom_bar(aes(fill = fill_group), stat = "identity", position = "fill") +
    coord_polar("y", start=0) + labs(x = NULL, y = NULL) + scale_y_reverse() +
    facet_wrap(~in_group, ncol = fitgrid(unique(ddf2$in_group))[2])
  }else{
    prop_table <- thisdf <- ddf2 %>% group_by(in_group) %>%
      mutate(
        fraction = value / sum(value),
        ymax = cumsum(fraction),# top of each rectangle
        ymin = c(0, head(ymax, n = -1)), # bottom of each rectangle
        position = ymax - fraction / 2,
        Percentage = paste0(round(fraction * 100, 1), "%")
      )
    p <- ggplot(thisdf, aes(ymax = ymax, ymin = ymin, xmax = 4, xmin = 3, fill = fill_group)) +
      geom_rect() +
      coord_polar(theta = "y") + # Try to remove that to understand how the chart is built initially
      facet_wrap(~in_group, ncol = fitgrid(unique(ddf2$in_group))[2]) +
      xlim(c(2, 4))
  }
  p <- p + theme_minimal() + labs(x=NULL)
  if(isTRUE(return_table)) return(list(plot = p, table = prop_table))
  return(p)
}

seu_heatmap <- function(
  object,
  annoc = NULL,
  rnames = NULL,
  cnames = NULL,
  orderby = NULL,
  seu_scale = TRUE,
  feature_order = FALSE,
  regress = c('nCount_RNA', 'percent.mt'),
  v = FALSE
){
  if(is.null(annoc)) annoc <- object@meta.data
  if(is.null(rnames)) rnames <- rownames(object@assays$RNA@scale.data)
  if(is.null(cnames)) cnames <- colnames(annoc)
  cnames <- getfound(cnames, rownames(annoc), v = v)
  cnames <- getfound(cnames, colnames(object), v = v)
  orderby <- orderby[1]
  if(v){
    cat("Features:", commas(rnames), "\n")
    cat("Cells:", commas(cnames), "\n")
    cat("Variables:", commas(colnames(annoc)), "\n")
  }
  annoc <- annoc[cnames, ]
  object <- object[, cnames]
  annoc <- annoc[order(annoc[, orderby]), ]
  tvar <- which(annoc[, orderby] == min(abs(annoc[, orderby])))
  if(length(tvar) != 0) tvar <- 1:nrow(annoc) %in% (tvar - 2):(tvar + 2)
  if(is.numeric(annoc[, orderby]) && length(tvar) != 0) annoc$Threshold = ifelse(tvar, "0", " ")
  anncolist <- lapply(annoc[, sapply(annoc, is.character), drop = FALSE], function(x){
    v2cols(select = x, fw = "gg", v = v)
  })
  if(is.numeric(annoc[, orderby]) && length(tvar) != 0) anncolist$Threshold <- c(" " = "#FFFFFF", "0" = "red")
  anncolist <- rev(anncolist)
  regress <- regress[regress %in% colnames(object@meta.data)]
  if(length(regress) == 0) regress <- NULL
  object <- object[rnames, rownames(annoc)]
  if(isTRUE(seu_scale) && casefold(class(object)) == 'seurat'){
    object <- ScaleData(
      object = object, features = rnames,
      vars.to.regress = regress,
      verbose = v
    )
    mat_to_plot <- object@assays$RNA@scale.data
  }else if(isTRUE(seu_scale)){
    if(v) cat("Scaling\n")
    if(casefold(class(object)) == 'seurat'){
      mat_to_plot <- as.matrix(expm1(object@assays$RNA@data) * 100)
    }
    mat_to_plot <- t(scale(t(mat_to_plot)))
  }
  if(isTRUE(feature_order) && casefold(class(object)) == 'seurat'){
    object <- RunPCA(object = object, npcs = 5, verbose = v)
    thesegenes <- rev(names(sort(Loadings(object, "pca")[, 1], decreasing = FALSE)))
  }else{ thesegenes <- rnames }
  mat_to_plot <- as.matrix(mat_to_plot[thesegenes, rownames(annoc)])
  topz <- max(c(min(abs(c(range(mat_to_plot), 2))), 1))
  mat_to_plot[mat_to_plot > topz] <- topz; mat_to_plot[mat_to_plot < (-topz)] <- -topz;
  range(mat_to_plot)
  NMF::aheatmap(mat_to_plot, annCol = annoc, annColors = anncolist,
    scale = 'none', Rowv = NA, Colv = NA, distfun = 'spearman', subsetCol = NULL,
    col = rev(colorRampPalette(colors = c('yellow', 'black', 'blue'))(7)))
}
