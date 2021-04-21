#' @title Dot plots
#'
#' @description `dot_plot` creates dot/curtain plots from feature measurements.
#' @details This is a function that takes measurements data and it's annotation
#' to calculate the mean and percentage, and displays them in a dot plot.
#'
#' @param edata Measurements matrix (rows: features/var, columns: observations).
#' @param mdata Annotation/metadata.
#' @param features Features from the measurements matrix (rows).
#' @param columns Columns from metadata.
#' @param scale_mean Should the means per feature be z-scored?, Default: TRUE.
#' @param cols Colours from the mean scale, Default: c("#fff4ba", "#ff0000").
#' @param clust_rows Cluster rows, Default: FALSE. Options: "rev"/"-1" to
#' revers the clustering, "right" to put the tree on the right.
#' @param clust_cols Same as 'clust_rows' but use "bottom", Default: FALSE.
#' @param annotate_cols Whether to annotate the columns, Default: FALSE. It can
#' also be columns in the metadata.
#' @param annotate_colours Colours to use for labels in the annotation,
#' Default: NULL. It cal also be a named (columns in metadata) list or a set
#' from RColorBrewer.
#' @param annotate_legend Whether to plot the annotation legend, Default: TRUE.
#' @param dot_legend Whether to plot the dot plot legend, Default: TRUE.
#' @param cols_limits Mean limits, Default: c(-2, 2).
#' @param size_limits Percentage/size of the dots limits, Default: c(10, max).
#' @param size_scale 'range' in scale_size/scale_radius), Default: 4.
#' @param scale_fun 'radius' or 'size', Default: 'radius'.
#' @param features_position Argument in scale_y_discrete, Default: 'left'.
#' @param values_trans Function to apply to the values in 'edata', Default: c.
#' @param verbose Show progress, Default: FALSE.
#' @return A patchwork object.
#'
#' @keywords z-score gene-expression
#'
#' @examples
#' \dontrun{
#'  if(interactive()){
#'   dot_plot(
#'    edata = obj@assays$RNA@data, mdata = obj@meta.data,
#'    columns = "cluster", features = mygenes, values_trans = expm1
#'   )
#'  }
#' }
#' @seealso
#'  \code{\link[reshape2]{c("cast", "melt")}}
#'  \code{\link[data.table]{rbindlist}}
#'  \code{\link[cowplot]{c("theme_cowplot", "plot_grid")}}
#'  \code{\link[base]{scale}}
#'  \code{\link[aplot]{c("xlim2", "xlim2")}}
#'  \code{\link[RColorBrewer]{RColorBrewer}}
#'  \code{\link[ggtree]{ggtree}}
#'  \code{\link[patchwork]{patchwork}}
#' @rdname dot_plot
#' @export
#' @importFrom reshape2 dcast melt
#' @importFrom data.table rbindlist
#' @importFrom cowplot theme_cowplot plot_grid
#' @importFrom scales squish
#' @importFrom aplot xlim2 ylim2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ggtree ggtree
#'
#' @author Ciro Ramírez-Suástegui, \email{ksuasteguic@gmail.com}
#' @references \url{http://davemcg.github.io/post/lets-plot-scrna-dotplots}
#' @references \url{https://divingintogeneticsandgenomics.rbind.io/post/clustered-dotplot-for-single-cell-rnaseq}
#' @references \url{https://uigradients.com}
#' @references \url{https://colorbrewer2.org}
dot_plot = function(
  edata,
  mdata,
  features,
  columns,
  scale_mean = TRUE,
  cols = c('#fff4ba', '#ff0000'),
  clust_rows = FALSE,
  clust_cols = FALSE,
  annotate_cols = FALSE,
  annotate_colours = NULL,
  annotate_legend = TRUE,
  dot_legend = TRUE, # dotplot legend
  cols_limits = c(-2, 2),
  size_limits = 10,
  size_scale = 6,
  scale_fun = "radius", # 'radius' is used by Seurat::DotPlot
  features_position = "left",
  values_trans = c, # expm1
  values_exprthr = 0, # expm1
  verbose = FALSE
){
  suppressPackageStartupMessages({
    library(tidyverse)
    library(patchwork)
  })
  # rev(viridis::magma(20)) # viridis::viridis(20)
  # c("#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#6baed6", "#4292c6", "#2171b5", "#08519c", "#08306b")
  # c("#ffffe5", "#fff7bc", "#fee391", "#fec44f", "#fe9929", "#ec7014", "#cc4c02", "#993404", "#662506")
  cols = colorRampPalette(colors = cols, space = 'Lab')(max(c(length(cols), 20)))
  any_xy = function(x, y) any(as.character(x) %in% as.character(y))
  cluster_axis = function(data, axes = c("Gene", "cluster", "count"), reorder = FALSE){
    mat <- reshape2::dcast(data = data, formula = paste0(axes[1], "~", axes[2]), value.var = axes[3])
    rownames(mat) <- mat[, axes[1]]; mat <- mat[,-1] %>% as.matrix() %>% dist()
    hc <- hclust(mat); if(any_xy(reorder, c("rev", "-1"))) hc$order <- rev(hc$order)
    return(hc)
  }
  set_limits = function(x, limits){
    if(is.null(limits)) limits <- range(x, na.rm = TRUE)
    # tvar <- x >= limits[1]; if(length(limits) == 1) tvar <- tvar & x <= limits[2]
    if(length(limits) == 1) limits <- c(limits, max(x, na.rm = TRUE))
    return(limits) # range(x[which(tvar)], na.rm = TRUE)
  }
  scale_fun <- switch(
    EXPR = scale_fun, 'size' = scale_size, 'radius' = scale_radius,
    stop("'size' or 'radius' are the options for 'scale_fun'")
  )

  # mdata$Identity = ident_combine(mdata, columns, "-") # from utilities
  mdata$Identity = if(length(columns) > 1){
    do.call(paste, c(mdata[, columns, drop = FALSE], sep = "-"))
  }else{ mdata[, columns] }
  efeatures = data.frame(t(as.matrix(edata[features, rownames(mdata)])), id = mdata$Identity)
  # Calculating the summary statistics
  gene_cluster <- lapply(X = levels(x = efeatures$id), FUN = function(ident) {
    slice <- efeatures[efeatures$id == ident, -ncol(efeatures), drop = FALSE]
    as.data.frame(list(
      count = apply(X = slice, MARGIN = 2, FUN = function(x) mean(x = values_trans(x)) ),
      pct = apply(X = slice, MARGIN = 2, FUN = function(x) mean(x = values_trans(x) > values_exprthr) ),
      Gene = colnames(slice), cluster = ident
    ))
  })
  gene_cluster = data.frame(data.table::rbindlist(gene_cluster))
  gene_cluster$pct <- gene_cluster$pct * 100
  gene_cluster$Gene <- factor(gene_cluster$Gene, rev(features))
  scale2 <- function(x) base::scale(x[!is.na(x)])[,1] # scaling mean
  if(scale_mean) gene_cluster <- gene_cluster %>% group_by(Gene) %>% mutate_at("count", scale2)
  if(verbose) print(summary(gene_cluster[, c("count", "pct")]))

  if(!isFALSE(clust_rows)){ # TRUE, rev|-1, right to control the position
    clustr = cluster_axis(data = gene_cluster, axes = c("Gene", "cluster", "count"), reorder = clust_rows)
    gene_cluster$Gene = factor(gene_cluster$Gene, levels = clustr$labels[clustr$order])
  }
  if(!isFALSE(clust_cols)){
    clustc = cluster_axis(data = gene_cluster, axes = c("cluster", "Gene", "count"), reorder = clust_cols)
    gene_cluster$cluster = factor(gene_cluster$cluster, levels = clustc$labels[clustc$order])
  }
  cols_limits = set_limits(gene_cluster$count, cols_limits)
  size_limits = set_limits(gene_cluster$pct, size_limits)

  if(verbose) cat("Creating dot-plot\n")
  dotplot <- gene_cluster %>%
    ggplot(aes(x = cluster, y = Gene, color = count, size = pct)) +
    geom_point() + cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
      axis.text.y = element_text(size = 13, face = "bold.italic"),
      axis.ticks = element_blank(), axis.line  = element_blank()
    ) +
    scale_color_gradientn(colours = cols, limits = cols_limits, oob = scales::squish) +
    scale_fun(range = c(0, size_scale), limits = size_limits) +
    labs(
      x = NULL, y = NULL,
      size = "Percent\nExpressed",
      color = paste0(ifelse(scale_mean, "Z-scored\n", ""), "Mean\nExpression")
    ) +
    scale_y_discrete(position = features_position) +
    scale_x_discrete(labels = function(l) parse(text=l))

  if(verbose) cat("Getting annotation\n");
  tvar = if(length(columns) == 1) c(annotate_cols[!annotate_cols %in% columns]) else annotate_cols
  tvar <- unique(tvar[!tvar %in% c("TRUE", "FALSE")])
  if(isTRUE(annotate_cols) || length(tvar) > 0) annotate_cols = TRUE
  if(length(tvar) == 0) tvar = "Identity"
  labels_l <- list(); labels_legend_l <- list()
  if(verbose) cat("Annotation:", paste0(tvar, collapse = ", "), "\n")
  for(i in tvar){
    labels_tab <- lapply(X = levels(x = mdata$Identity), FUN = function(ident) {
      slice <- mdata[mdata$Identity == ident, , drop = FALSE] # summarise_table from utilities
      if(!is.numeric(slice[, i])){
        paste0(levels(factormix(as.character(slice[, i]))), collapse = "|")
      }else{ mean(slice[, i]) }
    }); names(labels_tab) <- levels(x = mdata$Identity); labels_tab <- reshape2::melt(labels_tab)
    labels_tab$L1 <- factor(labels_tab$L1, levels = levels(gene_cluster$cluster))
    if(all(as.character(labels_tab$value) %in% as.character(mdata[, i]))) # keep order
      labels_tab$value <- factor(as.character(labels_tab$value), levels(factormix(mdata[, i])))
    labels_p <- labels_tab %>%
      ggplot(aes(x = L1, y = 1, fill = value)) + geom_tile() + theme_nothing() + aplot::xlim2(dotplot)
    groups_i = levels(labels_tab$value)
    guide_vars = list("right", "left", element_text(angle = 0))
    if(length(groups_i) > 7 && !is.numeric(mdata[, i])) guide_vars = list("bottom", "left", element_text(angle = 90))
    if(!is.numeric(mdata[, i]) && length(annotate_colours) > 0){ # getting colours
      col_var = if(i %n% names(annotate_colours)){
        annotate_colours[[i]]
      }else if(!is.list(annotate_colours)) annotate_colours
      if(!is.null(col_var)){
        col_var = if(length(col_var) == 1){
          brewer_n = c(Accent=8,Dark2=8,Paired=12,Pastel1=9,Pastel2=8,Set1=9,Set2=8,Set3=12)
          col_var = RColorBrewer::brewer.pal(n = brewer_n[col_var], name = col_var)
          colorRampPalette(colors = col_var, space = 'Lab')(length(groups_i))
        }else{
          col_var = if(any(groups_i %in% names(col_var))) col_var[groups_i] else setNames(col_var, groups_i)
        }; labels_p + scale_fill_manual(values = col_var)
      }
    }
    labels_legend_l[[i]] <- get_legend(
      labels_p + theme(legend.position = "bottom") + labs(fill = i) +
      guides(fill = guide_legend(
        label.position = guide_vars[[1]], label.theme = guide_vars[[3]],
        title.position = guide_vars[[2]], nrow = 1
      ))
    ); labels_l[[i]] <- labels_p + theme(legend.position = "none")
  }
  labels_legend = cowplot::plot_grid(plotlist = labels_legend_l, nrow = length(labels_legend_l))
  annotate_plot = cowplot::plot_grid(plotlist = labels_l, nrow = length(labels_l))
  dotplot_legend <- get_legend(dotplot)
  dotplot <- dotplot + theme(legend.position = "none")

  widths_i = c(row = 0.7, sep = -0.1, dot = 4)
  heights_i = c(col = 0.5, ann = length(labels_legend_l) / 10, sep = -0.1, dot = 4, leg = 1)
  if(!isFALSE(clust_rows)){
    if(verbose) cat("Showing row clusters\n")
    ggtree_plot_r <- suppressMessages(ggtree::ggtree(as.dendrogram(clustr)) + aplot::ylim2(dotplot))
    if(any_xy(clust_rows, c("rev", "-1"))) ggtree_plot_r <- suppressMessages(ggtree_plot_r + scale_y_reverse())
    if(any_xy(clust_rows, c("right"))) ggtree_plot_r <- ggtree_plot_r + scale_x_reverse()
  }else{ widths_i = NULL }
  if(!isFALSE(clust_cols)){
    if(verbose) cat("Showing column clusters\n") # ggtree::layout_dendrogram() # acts as clockwise flip
    ggtree_plot_c <- ggtree::ggtree(as.dendrogram(clustc)) + coord_flip() # face up
    ggtree_plot_c <- suppressMessages(ggtree_plot_c + aplot::xlim2(dotplot))
    if(!any_xy(clust_cols, c("bottom"))) ggtree_plot_c <- ggtree_plot_c + scale_x_reverse()
    if(any_xy(clust_cols, c("rev", "-1"))) ggtree_plot_c <- suppressMessages(ggtree_plot_c + scale_y_reverse())
  }else{ heights_i = heights_i[-1] }

  if(!(!isFALSE(clust_cols) || annotate_cols)) heights_i = heights_i[heights_i != -0.1]
  annotate_legend = annotate_cols && annotate_legend
  if(!annotate_cols) heights_i = heights_i[names(heights_i) != "ann"]
  if(annotate_cols) heights_i[length(heights_i)] = tail(heights_i, 1) + (length(labels_legend_l) / 2)
  if(!annotate_legend) heights_i = head(heights_i, -1)
  if(length(heights_i) == 1) heights_i = NULL
  rows = ifelse(!isFALSE(clust_rows), yes = "plot_spacer()+plot_spacer()+", no = "")
  command = c(
    ifelse(!isFALSE(clust_rows) && !isFALSE(clust_cols), rows, ""), # 1
    ifelse(!isFALSE(clust_cols), "ggtree_plot_c+", ""), # 1
    ifelse(!isFALSE(clust_rows) && annotate_cols, rows, ""), ifelse(annotate_cols, "annotate_plot+", ""), # 2
    ifelse(!isFALSE(clust_rows) && (!isFALSE(clust_cols) || annotate_cols), rows, ""), # 3
    ifelse(!isFALSE(clust_cols) || annotate_cols, "plot_spacer()+", ""), # 3
    ifelse(!isFALSE(clust_rows), "ggtree_plot_r+plot_spacer()+", ""), "dotplot+" # 4
  ); command = panels = paste0(unlist(strsplit(command, "\\+")), "+")
  tvar <- if(!isFALSE(clust_rows)) 3 else 1
  command_grid = init = matrix(command, ncol = tvar, byrow = TRUE)
  # tvar <- if(!isFALSE(clust_rows)) ceiling(sqrt(length(command))) else length(command)
  # command_grid = init = matrix(command, nrow = tvar, byrow = TRUE)
  if(any_xy(clust_cols, c("bottom"))){
    heights_i <- rev(heights_i); command_grid <- command_grid[nrow(command_grid):1, , drop = FALSE]
  }
  if(any_xy(clust_rows, c("right"))){
    widths_i <- rev(widths_i); command_grid <- command_grid[, ncol(command_grid):1, drop = FALSE]
  }
  if(dot_legend){
    tvar <- which(command_grid == "dotplot+", arr.ind = TRUE)
    command_grid = cbind(command_grid, rep("plot_spacer()+", nrow(command_grid)))
    command_grid[tvar[1,1], ncol(command_grid)] <- "dotplot_legend+"
    widths_i = c(if(is.null(widths_i)) 4 else widths_i, 0.7)
  }; command_grid[, ncol(command_grid)] <- paste0(command_grid[, ncol(command_grid)], "\n")

  command = paste0(c("pp <- ", if(any(dim(command_grid)>1)) t(command_grid) else command_grid), collapse = "")
  command = paste0(command, ifelse(!isFALSE(clust_rows) && annotate_legend, rows, ""))
  command = paste0(command, ifelse(annotate_legend, "labels_legend+\n", ""))
  command = paste0(command, "plot_layout(ncol=max(c(length(widths_i),1)),widths=widths_i,heights=heights_i)")
  if(verbose) cat(command, "\n")
  eval(expr = parse(text = command))
  # return(list(plot = pp, panels = panels, grid_init = init, grid = command_grid,
  #   command = command, w = widths_i, h = heights_i))
  return(pp)
}
