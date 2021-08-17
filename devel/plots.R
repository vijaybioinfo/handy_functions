#!/usr/bin/R

suppressPackageStartupMessages({
  library(ggplot2)
  library(cowplot)
})
theme_set(theme_cowplot())

# Colour brewer: colorbrewer2.org
# Colour brewer expanded: https://bl.ocks.org/emeeks/8cdec64ed6daf955830fa723252a4ab3
# Compare range between two colours: https://learnui.design/tools/data-color-picker.html
# Viz your custom palettes: https://projects.susielu.com/viz-palette
# plotrix::color.id("#800020") # hex2col
# grDevices::colorRampPalette("firebrick3")(1) # col2hex
# scales::show_col(sample(x = colors(), size = 7)) # show_colours

couls_opt = list(
  red_gradient = list(
    strong = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000'),
    strong_white = c('#fffffa', '#ffdf32', '#ff9a00', '#ff5719','#EE0000','#b30000', '#670000'),
    white = c('#fffffa', '#fff7cf', '#ffdf32', '#ff9a00', '#EE0000','#b30000', '#670000'), # white
    divisive = c('#fff7cf', '#fcf193', '#ffdf32', '#ff9a00', '#EE0000','#b30000', '#670000'),
    brewer = c("#ffffcc","#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026")),
  green_gradient = list(
    brewer = c("#f7fcf5", "#e5f5e0", "#c7e9c0", "#a1d99b", "#74c476", "#41ab5d", "#238b45", "#006d2c", "#00441b"),
    brewer_yellow = c("#ffffe5", "#f7fcb9", "#d9f0a3", "#addd8e", "#78c679", "#41ab5d", "#238443", "#006837", "#004529")),
  purple_gradient = list(
    brewer = c("#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6", "#8c6bb1", "#88419d", "#810f7c", "#4d004b"),
    brewer_yellow = c("#ffffd9", "#edf8b1", "#c7e9b4", "#7fcdbb", "#41b6c4", "#1d91c0", "#225ea8", "#253494", "#081d58")),
  blue_gradient = list(
    brewer = c("#f7fcf0", "#e0f3db", "#ccebc5", "#a8ddb5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081"))
)

plot_color_list = function(x, col_names = FALSE){
  x = unlist(x, recursive = FALSE)
  n <- max(sapply(x, length))
  colorlist = lapply(x, function(y){
    z <- y[1:n]; z[is.na(z)] <- "#FFFFFF"; z
  })
  names(colorlist) <- gsub("_gradient|gradient", "", names(colorlist))
  n <- rep(n, length(colorlist))
  nr <- length(colorlist)
  nc <- max(n)
  ylim <- c(0, nr)
  oldpar <- par(mgp = c(2, 0.25, 0))
  on.exit(par(oldpar))
  par(mar = c(0.5, 8, 0.5, 0.5))
  plot(1, 1, xlim = c(0, nc), ylim = ylim, type = "n", axes = FALSE,
      bty = "n", xlab = "", ylab = "")
  for (i in 1:nr) {
      nj <- n[i]
      shadi <- colorlist[[i]]
      border_i = ifelse(shadi != "#FFFFFF", "light grey", shadi)
      rect(xleft = 0:(nj - 1), ybottom = i - 1, xright = 1:nj,
          ytop = i - 0.2, col = shadi, border = border_i)
      if(isTRUE(col_names))
        text((1:nc) - 0.5, rep(i - 0.6, nj), labels = shadi, cex = 0.5, xpd = TRUE)
  }
  text(rep(-0.1, nr), (1:nr) - 0.6, labels = names(colorlist), xpd = TRUE, adj = 1)
}
# pdf("lab_colors.pdf", width = 10, height = 7)
# plot_color_list(couls_opt)
# dev.off()

# Packages: ggplot2, cowplot

#' Violin Plot
#'
#' This function creates nice violin plots
#'
#' @param dat data.frame with the data to plot.
#' @param xax Groups.
#' @param yax Y axis.
#' @param dots Whether to show all the data points (can be a number).
#' @param colour_by Column to colour by. Options are: 'pct' (0-100%), 'pct50'
#' (0-50%) or 'means' (takes the range of 'yax') or another column
#' (takes the mean). It also accepts a structure:
#' list(column_in_dat = list(legend_name = '%+', breaks = c(0-9)).
#' @param couls Colours you want to use. It can be a vector length >= 2.
#' If one colour is given, then the gradient is from white to the colour
#' @param vsetting Violin plot setting: width, alpha, trim, adjust, and scale
#' @keywords GO
#' @references \url{https://no/url/given.pdf}
#'
#' @return Returns ggplot object
#'
#' @importFrom ggplot2
#' @importFrom cowplot theme_cowplot
#'
#' @export
#'
#' @examples
# ' p <- violin(dat = iris, xax = "Species", yax = "Sepal.Length")
#'

violin <- function(
  dat,
  xax = 1,
  yax = 2,
  dots = FALSE,
  colour_by = NULL, # pct, means (based on 'xax')
  couls = NULL,
  vsetting = NULL,
  chilli = FALSE
){
  if(is.null(colour_by)) colour_by <- xax
  if(is.numeric(xax)) xax <- colnames(dat)[xax]
  if(is.numeric(yax)) xax <- colnames(dat)[yax]
  colnames(dat) <- make.names(colnames(dat)) # dashes make it fail...
  xax <-  make.names(xax)
  yax <-  make.names(yax)
  colour_by <-  make.names(colour_by)

  dat <- dat[order(dat[, xax]), ]
  dat[, yax] <- as.numeric(dat[, yax])
  ylims <- range(dat[, yax])

  vsetting_master <- list(width = 0.8, alpha = 0.7, trim = TRUE, adjust = 1, scale = 'width')
  #Deafault: trim = TRUE, scale = "area", adjust = 1
  vsetting <- lapply(names(vsetting_master), function(x){
    ifelse(x %in% names(vsetting), vsetting[[x]], vsetting_master[[x]])
  })
  names(vsetting) <- names(vsetting_master)

  if(isTRUE(colour_by %in% colnames(dat))){
    tvar <- length(unique(dat[, colour_by])) > length(unique(dat[, xax]))
    if(is.numeric(dat[, colour_by]) && tvar){
      means <- tapply(dat[, colour_by], dat[, xax], mean, na.rm = TRUE)
      dat$dummy123 <- sapply(as.character(dat[, xax]), function(X) return(means[X]) )
      colnames(dat) <- gsub("dummy123", paste0(colour_by, "_avg"), colnames(dat));
      colour_by <- paste0(colour_by, "_avg")
    }
  }; filname <- colour_by
  pct <- tapply(dat[, yax], dat[, xax], function(pop) round(sum(pop > 0, rm.na = TRUE) / length(pop) * 100, digits = 2) )
  pct[pct>100] <- 100
  dat$pct <- dat$pct50 <- sapply(as.character(dat[, xax]), function(X) return(pct[X]) )
  dat$pct50[dat$pct50>50] <- 50
  means <- tapply(dat[, yax], dat[, xax], mean, na.rm = TRUE)
  dat$mean <- sapply(as.character(dat[, xax]), function(X) return(means[X]) )
  dat$mean_full <- dat$mean_12 <- dat$mean

  if(is.null(couls)) couls = 1
  if(is.numeric(couls)) couls = couls_opt$red_gradient[[couls]]
  brks <- NULL
  master_brks = list( # create a master for breaks in mean and percentage
    pct = list(filname = '%+', brks = c(0, 10, 20, 40, 60, 80, 100)),
    mean = list(filname = 'Mean', make_breaks(dat$means)),
    mean_full = list(filname = 'Mean', make_breaks(dat[, yax])),
    mean_12 = list(filname = 'Mean', brks = c(0, 2, 4, 6, 8, 10, 12)),
    pct50 = list(filname = '%+cells', brks = c(0, 5, 10, 20, 30, 40, 50))
  ) # substitute the master with the given breaks
  if(is.list(colour_by)){
    if(is.null(names(colour_by))) warning("Breaks must be a named list of list.")
    master_brks[[names(colour_by)]] <- colour_by; colour_by <- names(colour_by)
  }
  if(isTRUE(colour_by %in% names(master_brks))){
    filname <- master_brks[[colour_by]][[1]]
    brks <- master_brks[[colour_by]][[2]]
  }; #if(colour_by == "pct50") colour_by <- "pct"

  if(isTRUE(chilli) || is.numeric(chilli)){
    if(isTRUE(chilli)) chilli <- 0
    dat <- dat[dat[, yax] > chilli, ]
  }

  dp <- ggplot(dat, aes_string(x = xax, y = yax, fill = colour_by))
  if(isTRUE(dots) || is.numeric(dots)){
    #round(min(c(ifelse(nrow(d2show) < 20000, nrow(d2show) / 3, nrow(d2show) / 5), 15000)))
    if(isTRUE(dots)) dots <- 1000
    set.seed(27);
    d2show <- dat[unname(sample_even(dat, xax, -dots)), ]
    # jittered_x <- unlist(lapply(make_list(dat, xax), function(g){
    #   vipor::offsetX(y = dat[g, yax], x = dat[g, xax])
    # }))
    dp <- dp + geom_jitter(
      data = d2show, mapping =  aes_string(x = xax, y = yax), inherit.aes = FALSE,
      size = 0.3, shape = 16, position = position_jitter(0.2), color = 'black'
    )
  }
  dp <- dp + geom_violin(
      aes_string(color = colour_by),
      alpha = vsetting[['alpha']],
      trim = vsetting[['trim']],
      adjust = ifelse(length(table(dat[, xax])) > 1, vsetting[['adjust']], 1),
      scale = vsetting[['scale']]
    ) +
    geom_boxplot(width=0.1, fill = "white", alpha = 0.25, outlier.shape = NA)+
    cowplot::theme_cowplot() + ylim(ylims)
  if(is.numeric(dat[, colour_by])){
    if(is.null(brks)){
      brks <- if(colour_by == "pct") c(0, dat[, colour_by], 100) else dat[, colour_by]
      brks <- make_breaks(brks, n = length(couls))
    };# dat[dat[, colour_by] > max(brks), colour_by] <- max(brks)
    couls <- colorRampPalette(colors = couls, space = 'Lab')(length(brks))
    lbrks <- as.character(round(brks, 2))
    if(any(grepl("\\-", lbrks))) lbrks <- ifelse(!grepl("\\-", lbrks), paste0("+", lbrks), lbrks)
    lbrks[grepl("\\+0$", lbrks)] <- "  0"
    arg_list <- list(name = filname, colours = couls, na.value = 'transparent',
      breaks = brks, labels = lbrks, limits = c(min(brks), max(brks)), guide = "colorbar")
    dp <- dp + guides(color = guide_colorbar(barwidth = 0.8, barheight = 6)) +
        do.call(scale_color_gradientn, arg_list) +
        do.call(scale_fill_gradientn, arg_list)
  }else if(length(unique(dat[, colour_by])) < 10){
    dp <- dp + scale_fill_brewer(palette = "Set1")
    dp <- dp + scale_color_brewer(palette = "Set1")
  }
  dp
}

violins <- function(
  dat,
  xax = 1,
  yax = 2,
  combine = TRUE,
  ncols = NULL,
  add_legend = TRUE,
  extra_theme = theme(plot.title = element_text(face = "bold.italic")), # usefull when combining internally
  ...
){
  if(length(yax) > 1){
    pp <- lapply(yax, function(x) violin(dat = dat, xax = xax, yax = x, ...) )
    pp <- lapply(1:length(yax), function(x) pp[[x]] + labs(title = yax[x], x = NULL, y = NULL) )
    if(!is.null(extra_theme)) pp <- lapply(1:length(yax), function(x) pp[[x]] + extra_theme )
    if(isTRUE(combine) || !is.null(ncols)){
      arrangge(plotlist = pp, ncols = ncols)
    }
    if(isTRUE(add_legend)) pp <- addlegend(pp)
    if(isTRUE(combine)) pp <- cowplot::plot_grid(plotlist = pp, ncol = ncols)
    return(pp)
  }else{
    return(violin(dat = dat, xax = xax, yax = yax, ...))
  }
}
arrangge <- function(
  plotlist,
  ncols = 1
) {
  if(is.null(ncols)) ncols <- make_grid(length(plotlist))[2]
  lapply(1:length(plotlist), function(x){
    if(isTRUE(x <= round(length(plotlist) - ncols))){
      plotlist[[x]] + theme(
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()
      )
    }else{ plotlist[[x]] }
  })
}
getlegend <- function(
  x,
  myfile = NULL
){
  # https://wilkelab.org/cowplot/articles/shared_legends.html
  legend <- cowplot::get_legend(x)#graphics.off();
  if(!is.null(myfile)){
    pdf(myfile, 5, 5)
    grid.newpage()
    grid.draw(legend)
    graphics.off()
    return(NULL)
  }
  return(legend)
}
addlegend <- function(x){
  p <- getlegend(x[[1]]); x <- lapply(x, function(y) y + Seurat::NoLegend() ); x$legend <- p
  return(x)
}

custom_heatmap <- function(
  object, # should be normalised
  annoc = NULL,
  rnames = NULL,
  cnames = NULL,
  use_mean = FALSE, # can be a vector of the columns to group
  orderby = NULL, # first calculates the mean if needed
  sample_it = FALSE, # can be c("column name", "-num_per_group"); or will use the first column
  scale_row = TRUE,
  categorical_col = "orig",
  feature_order = FALSE, # logical, 'pca', 'hclust'
  couls = NULL, # colors for columns and rows
  hcouls = c('yellow', 'black', 'blue'),
  regress = c('nCount_RNA', 'percent.mt'),
  topz = NULL,
  verbose = FALSE,
  type = NULL,
  cluster_cols_override = FALSE,
  do_log = FALSE,
  return_objects = FALSE,
  ... # arguments for NMF::aheatmap or pheatmap
){
  if(casefold(class(object)) == 'seurat'){
    if(verbose) cat("Seurat object\n")
    if(is.null(annoc)) annoc <- object@meta.data
  }else if(is.null(annoc)){ stop("Annotation not provided") }
  if(is.null(rnames)) rnames <- rownames(object)
  if(is.null(cnames)) cnames <- colnames(object)
  cnames <- show_found(cnames, rownames(annoc), verbose = FALSE)
  cnames <- show_found(cnames, colnames(object), verbose = FALSE)
  if(!any(categorical_col %in% colnames(annoc))){
    categorical_col <- filters_columns(annoc, include = categorical_col, v = verbose)
  }; categorical_col <- categorical_col[categorical_col %in% colnames(annoc)]
  if(is.null(orderby)) orderby <- categorical_col[1]

  if(isTRUE(sample_it[1]) || is.character(sample_it)){
    nsample <- if(length(sample_it) > 1) as.numeric(sample_it[2]) else NULL
    sample_it <- if(isTRUE(sample_it[1])) categorical_col[1] else sample_it[1]
    cnames <- sample_even(annot = annoc, cname = sample_it, maxln = nsample, v = verbose)
  }

  if(verbose){
    cat("\nCategories:", show_commas(categorical_col), "\n")
    cat("Features:", show_commas(rnames), "\n")
    cat("Cells:", show_commas(cnames), "\n")
    cat("Variables in metadata:", show_commas(colnames(annoc)), "\n\n")
  }

  tvar <- unique(c(orderby, categorical_col, use_mean)); tvar <- tvar[is.character(tvar)]
  annoc <- annoc[cnames, tvar[tvar %in% colnames(annoc)], drop = FALSE]
  for(i in colnames(annoc)) if(is.factor(annoc[, i])) annoc[, i] <- droplevels(annoc[, i])
  object <- object[, cnames]
  if(verbose){ str(annoc); str(object, max.level = 1) }

  if(use_mean %in% colnames(annoc) || isTRUE(use_mean)){
    if(casefold(class(object)) == 'seurat'){
      if(verbose) cat("Taking only the matrix from Seurat object\n")
      object <- expm1(object@assays$RNA@data)[rnames, rownames(annoc)]
    }
    use_mean <- if(isTRUE(use_mean[1])) head(colnames(annoc)[sapply(annoc, is.character)], 1) else use_mean
    if(verbose) cat("Means:", paste0(use_mean, collapse = ", "), "\n")
    # annoc$index123 <- paste0("X", do.call(paste, c(annoc[, use_mean, drop = FALSE], sep = "_")))
    annoc$index123 <- do.call(paste, c(annoc[, use_mean, drop = FALSE], sep = "_"))
    means <- stats_summary_table(
      mat = object,
      rnames = rnames,
      groups = make_list(annoc, colname = 'index123', grouping = TRUE),
      moments = "mn",
      v = verbose
    )
    colnames(means) <- gsub("_mean", "", colnames(means))
    # annoc <- annoc[!duplicated(annoc$index123), ]
    tvar <- reshape2::melt(table(annoc[, rev(colnames(annoc)), drop = FALSE]))
    tvar <- tvar[tvar$value > 0, -ncol(tvar), drop = FALSE]
    for(i in colnames(annoc)) if(is.factor(annoc[, i])) tvar[, i] <- factor(as.character(tvar[, i]), levels(annoc[, i]))
    annoc <- summarise_table(tvar)
    annoc <- annoc[, !colnames(annoc) %in% "index123", drop = FALSE]
    object <- as.matrix(means[, rownames(annoc)])
    if(verbose){ str(annoc); str(object) }
  }

  if(verbose) cat("Ordering:", orderby)
  if(length(orderby) > 1 && any(c("pca", "hc") %in% orderby)){
    if(verbose) cat("\nbased on expression\n")
    tvar <- order_df(
      annot = annoc,
      order_by = ifelse(length(orderby[-1]) == 0, "pca", orderby[-1]),
      cname = orderby[1],
      grupos = names(table(annoc[, orderby[1]])),
      mat =
        if(casefold(class(object)) == 'seurat'){
          as.matrix(object@assays$RNA@data[rnames, rownames(annoc)])
        }else{ object[rnames, rownames(annoc)] },
      verbose = verbose
    )
    annoc <- annoc[unname(unlist(tvar)), ]
  }else{
    if(verbose) cat("\norder(annoc[, ", orderby, "])\n")
    # annoc <- annoc[order(annoc[, orderby]), , drop = FALSE] # will follow the factors order
    data.table::setorderv(x = annoc, cols = orderby)
  }

  tvar <- sapply(annoc, is.character) | sapply(annoc, is.factor)
  if(verbose) cat("Getting colours:", show_commas(colnames(annoc)[tvar]), "\n")
  anncolist <- lapply(annoc[, tvar, drop = FALSE], function(x){
    v2cols(select = x, sour = couls, v = FALSE)
  })
  row_params <- c('annRow', 'annotation_row')
  if(any(row_params %in% names(list(...)))){
    if(verbose) cat("Adding row annotation\n")
    tvar <- data.frame(
      list(...)[[which(names(list(...)) %in% row_params)]],
      check.names = FALSE, stringsAsFactors = FALSE
    ); str(tvar)
    # tmp <- colnames(tvar[, !sapply(tvar, is.numeric), drop = FALSE])
    # for(i in tmp){ # Factor levels on variable X do not match with annotation_colors
    #   anncolist[[i]] <- v2cols(
    #     select = unique(c(tvar[, i], anncolist[[i]])), sour = couls, v = verbose
    #   )
    # }
    if(!any(rownames(tvar) %in% rnames)) stop("row annotation needs row names!")
    anncolist <- c(anncolist,
      lapply(tvar[, !sapply(tvar, is.numeric), drop = FALSE], function(x){
        v2cols(select = x, sour = couls, v = verbose)
      })
    )
  }
  if(verbose){ cat("Colors:\n"); str(anncolist) }

  object <- object[rnames, rownames(annoc)]
  if(isTRUE(do_log) && class(object) == "matrix"){
    if(verbose) cat("log2(x+1)\n")
    object <- log2(object + 1)
  }

  objects2return = list()
  objects2return$annotation_col = annoc

  if(verbose) cat("Scaling\n")
  if(isTRUE(scale_row) && (casefold(class(object)) == 'seurat')){
    regress <- regress[regress %in% colnames(object@meta.data)]
    if(length(regress) == 0) regress <- NULL
    object <- ScaleData(
      object = object, features = rnames,
      vars.to.regress = regress,
      verbose = verbose
    )
    mat_to_plot <- object@assays$RNA@scale.data
  }else if(isTRUE(scale_row)){
    if(casefold(class(object)) == 'seurat'){
      if(verbose) cat("RNA@data from Seurat\n")
      mat_to_plot <- as.matrix(expm1(object@assays$RNA@scale.data)) # safe matrix on subsetted object?
    }else{
      if(verbose) cat("Expression matrix\n")
      mat_to_plot = object; object = "tiny"; gc()
    }
    objects2return$matrix = mat_to_plot
    mat_to_plot <- t(scale(t(mat_to_plot)))
  }

  nnpcs <- min(c(3, dim(mat_to_plot) - 1))
  thesegenes <- if((as.character(feature_order) == "pca") && (casefold(class(object)) == 'seurat')){
    if(verbose) cat("Ordering features Seurat PCA\n")
    object <- RunPCA(object = object, npcs = nnpcs, features = rnames, verbose = verbose)
    rev(names(sort(Loadings(object, "pca")[, 1], decreasing = FALSE)))
  }else if(feature_order[1] == "pca"){
    if(verbose) cat("Ordering features PCA\n")
    pc1 <- irlba::prcomp_irlba(mat_to_plot, n =  )$x[, 1]
    names(pc1) <- rownames(mat_to_plot)
    rev(names(sort(pc1, decreasing = T)))
  }else if(any(feature_order %in% rnames)){
    feature_order
  }else{ rnames }
  feature_order <- if(as.character(feature_order[1]) == "hclust"){
    if(verbose) cat("Ordering features based on hierarchical ordering\n"); TRUE
  }else{ NA }
  objects2return$features = feature_order

  if(verbose) cat("Creating heatmap plot\n")
  thesegenes <- c(thesegenes[thesegenes %in% rnames], rnames[!rnames %in% thesegenes])
  mat_to_plot <- as.matrix(mat_to_plot[thesegenes, rownames(annoc)])

  if(is.null(topz)) topz <- max(c(min(abs(c(range(mat_to_plot), 2))), 1))
  if(verbose) cat("Range:", range(mat_to_plot), "\n")
  mat_to_plot[mat_to_plot > topz] <- topz; mat_to_plot[mat_to_plot < (-topz)] <- -topz;
  if(verbose) str(mat_to_plot)

  mypalette <- colorRampPalette(colors = rev(hcouls), space = 'Lab')
  palettebreaks <- seq(from = -topz, to = topz, by = 0.1)

  if(isTRUE(grepl("NMF", type))){
    cat("Using NMF\n")
    y <- NMF::aheatmap(
      x = mat_to_plot,
      scale = 'none', Rowv = feature_order, Colv = NA,
      annCol = annoc, annColors = anncolist,
      col = mypalette(length(palettebreaks) - 1),#7)),
      ...
    )
  }else if(isTRUE(grepl("gg", type))){
    ddf = reshape2::melt(mat_to_plot)
    ddf$Var1 <- factor(ddf$Var1, rev(rownames(mat_to_plot)))
    ddf$Var2 <- factor(ddf$Var2, colnames(mat_to_plot))
    y <- ggplot(ddf, aes(x = Var2, y = Var1)) +
      geom_tile(aes(fill = value)) +
      scale_fill_gradientn(colours = mypalette(length(palettebreaks) - 1)) +
      labs(fill = NULL, x = NULL, y = NULL)
    # if(){ # for annotation maybe look at dot_plot
    #   ggplot(data = annoc, aes(x = L1, y = 1, fill = value)) + geom_tile() + theme_nothing()
    # }
  }else{
    cat("Using pheatmap\n")
    source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/pheatmapCorrection.R')
    if(is.na(feature_order) || is.character(feature_order)) feature_order <- FALSE
    y <- pheatmap(
      mat = mat_to_plot,
      scale = "none", cluster_rows = feature_order, cluster_cols = cluster_cols_override,
      annotation_col = annoc, annotation_colors = anncolist,
      color = mypalette(length(palettebreaks) - 1),
      breaks = palettebreaks,
      border_color = NA,
      drop_levels = TRUE,
      ...
    )
  }

  if(isTRUE(return_objects)){
    return(objects2return)
  }else{ return(y) }
}

plot_rm_layer <- function(gp, lays = "ext", verbose = FALSE){
  if(verbose) cat("Layers:", sapply(gp$layers, function(x) class(x$geom)[1] ), "\n")
  tvar <- sapply(gp$layers, function(x) !grepl(lays, class(x$geom)[1], ignore.case = TRUE) )
  gp$layers <- gp$layers[tvar]
  gp
}
plot_blank <- function(gp, lays = "ext", ...){
  if(!is.null(gp$facet)) gp <- gp + theme(strip.text.x = element_blank())
  if("patchwork" %in% class(gp) && isTRUE(length(gp$patches$plots) > 0)){
    for(i in 1:length(gp$patches$plots)){
      y <- gp$patches$plots[[i]]; if(length(y$layers) == 0) next
      gp$patches$plots[[i]] = plot_rm_layer(y, lays = lays, ...) + theme_axes()
    }
    gp
  }else{
    plot_rm_layer(gp, lays = lays, ...) + theme_axes()
  }
}

plot_flush <- function(x){ graphics.off(); invisible(file.remove('Rplots.pdf')) }

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
  dcount <- dnorm <- table(x[, groups])
  if(v) cat("Type of plot:", type, "\n")
  if(v) cat("Proportions of:", groups[1], "\n")
  if(v) cat("Per:", groups[2], "\n")
  if(isTRUE(print_ptables)) print(dnorm)
  if(isTRUE(normalise)) normalise <- 1
  if(is.character(normalise)) normalise <- which(groups == normalise)
  is_normalised <- NULL
  if(is.numeric(normalise)){
    is_normalised = "Normalised"
    if(v) cat("Normalising by", groups[normalise], "\n") # changed 2021-02-11
    # dnorm <- (dnorm * (min(table(x[, groups[normalise]])) / rowSums(dnorm))) # factorising/normalising
    denominator <- min(table(x[, groups[normalise]]))
    dnorm <- apply(
      X = dnorm,
      MARGIN = normalise,
      function(x){
        y <- x * (denominator / sum(x))
    }); if(normalise == 1) dnorm <- t(dnorm)
    if(isTRUE(print_ptables)) print(apply(dnorm, normalise, sum))
  }
  if(isTRUE(print_ptables)) print(dnorm)
  dprop <- prop.table(dnorm, 2); if(isTRUE(print_ptables)) print(dprop)
  mylevels <- if(is.numeric(orderby) || all(orderby %in% names(cats))){
    direct = FALSE
    orderby <- if(is.numeric(orderby)){
      direct = any(orderby > 0)
      names(cats[abs(orderby)])
    }else{
      orderby
    }
    if(v) cat("Sorting by", orderby, "\n")
    tvar <- if(length(orderby) > 1) colSums(dprop[orderby, ]) else dprop[orderby, ]
    names(sort(tvar, decreasing = direct))
  }else{ NULL }
  if(!is.null(is_normalised)) dprop <- rbind(dnorm, dprop)
  dprop <- rbind(dcount[, colnames(dprop)], dprop)
  dprop <- round(as.data.frame.matrix(dprop), 2)
  dprop <- cbind(
    Group = rep(rownames(dnorm), times = ifelse(is.null(is_normalised), 2, 3)),
    Type = rep(c("Count", is_normalised, "Percentage"), each = nrow(dnorm)),
    dprop
  )
  rownames(dprop) <- NULL

  dnorm <- reshape2::melt(dnorm)
  dnorm$fill_group <- if(!is.factor(x[, groups[1]])) factormix(dnorm[, 1]) else factor(dnorm[, 1], levels(x[, groups[1]]))
  dnorm$in_group <- if(!is.null(mylevels)){
    factor(dnorm[, 2], mylevels)
  }else if(!is.factor(x[, groups[2]])){
    factormix(dnorm[, 2])
  }else{
    factor(dnorm[, 2], levels(x[, groups[2]]))
  }
  dnorm <- dnorm[, -c(1:2)]
  if(type == "bar"){
    p <- ggplot(data = dnorm, aes(x = in_group, y = value, fill = fill_group)) +
      geom_bar(stat = "identity", position = "fill") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(fill = make_title(groups[1]), y = "Percentage")
  }else if(type == "pie"){
    p <- ggplot(data = dnorm, aes(x = "", y = value)) +
    geom_bar(aes(fill = fill_group), stat = "identity", position = "fill") +
    coord_polar("y", start=0) + labs(x = NULL, y = NULL) + scale_y_reverse() +
    facet_wrap(~in_group, ncol = make_grid(unique(dnorm$in_group))[2])
  }else{
    thisdf <- dnorm %>% group_by(in_group) %>%
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
      facet_wrap(~in_group, ncol = make_grid(unique(dnorm$in_group))[2]) +
      xlim(c(2, 4))
  }
  p <- p + theme_minimal() + labs(x=NULL)
  if(isTRUE(return_table)) return(list(plot = p, table = dprop))
  return(p)
}

plot_axis_cross = function(
    gp, column_names = NULL, limit = NULL, verbose = FALSE
){
  column_map <- gsub("~", "", as.character(gp$mapping))
  if(length(column_map) == 0) column_map = colnames(gp$data)
  if(verbose) cat("Mapping:", paste(column_map, collapse = ", "), "\n")
  if(is.null(column_names)) column_names <- column_map
  column_names <- column_names[sapply(gp$data[, column_names], is.numeric)]
  if(verbose) cat("Using:", paste(column_names, collapse = ", "), "\n")
  if(length(column_names) < 2) return(gp)
  dff <- gp$data[, column_names, drop = FALSE]
  if(all(sapply(dff, is.numeric))){
    if(!is.numeric(limit)) lx <- min(dff[, 1], na.rm = TRUE)
    if(!is.numeric(limit)) ly <- min(dff[, 2], na.rm = TRUE)
    len = min(mean(range(dff[, 1]))/3, mean(range(dff[, 2]))/3)
    if(verbose) cat("Cross:", lx, ly, ", length=", len, "\n")
    scatter_sizes = data.frame(
      x = c(lx, lx), xend = c(lx + len, lx),
      y = c(ly, ly), yend = c(ly, ly + len))
    gp <- gp + geom_segment(
        data = scatter_sizes, mapping = aes(x=x, y=y, xend=xend, yend=yend),
        inherit.aes = FALSE
      ) + coord_cartesian(clip = "off")
  }
  gp + theme(
    line = element_blank(),
    axis.title.x = element_text(hjust = 0.1), axis.title.y = element_text(hjust = 0.1),
    axis.text.x = element_blank(), axis.text.y = element_blank()
  )
}

plot_png = function(x){
  # conda install -c conda-forge poppler
  # export R_LD_LIBRARY_PATH=$R_LD_LIBRARY_PATH:/home/ciro/bin/miniconda3/lib
  # install.packages("pdftools")
  # sudo apt-get -y install libmagick++-dev
  # install.packages("magick")
  # devtools::install_github("ropensci/magick")
  img = switch(
    EXPR = casefold(tools::file_ext(x)),
    png = png::readPNG(x),
    jpeg = jpeg::readJPEG(x),
    pdf = magick::image_read_pdf(x)
    # pdf = pdftools::pdf_render_page(x)
  )
  qplot(1:10, 1:10, geom="blank") +
    annotation_custom(grid::rasterGrob(img, interpolate = TRUE),
      xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
}

### Aesthetic ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theme_axes <- function (
  base_size = 11,
  base_family = "",
  base_line_size = base_size / 22,
  base_rect_size = base_size / 22
){
  half_line <- base_size / 2
  t <- theme(
    line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "butt"),
    rect = element_rect(fill = "white",colour = "black", size = base_rect_size, linetype = 1),
    text = element_text(
      family = base_family, face = "plain",
      colour = "black", size = base_size, lineheight = 0.9,
      hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
      debug = FALSE
    ),
    axis.title.x = element_blank(), axis.title.y = element_blank(),
    axis.text.x = element_blank(), axis.text.y = element_blank(),
    legend.position = "none",
    plot.title = element_blank(), plot.subtitle = element_blank(), plot.caption = element_blank(),
    strip.text = element_blank(), strip.background = element_blank()
  )
  # ggplot_global$theme_all_null %+replace% t
}

rotatedAxisElementText = function(angle, position = 'x', cnter = 0.5){
  angle     = angle[1];
  position  = position[1]
  positions = list(x=0,y=90,top=180,right=270)
  if(!position %in% names(positions))
    stop(sprintf("'position' must be one of [%s]",paste(names(positions),collapse=", ")),call.=FALSE)
  if(!is.numeric(angle))
    stop("'angle' must be numeric",call.=FALSE)
  rads  = (angle - positions[[ position ]])*pi/180
  hjust = cnter*(1 - sin(rads))
  vjust = cnter*(1 + cos(rads))
  element_text(angle = angle,vjust = vjust, hjust = hjust)
}

theme_grid <- function(){
  theme_classic() + theme(
    panel.background = element_blank(), panel.grid.major = element_blank(),
    strip.text.x = element_text(face = 'bold'),
    strip.background = element_rect(fill = "#FFFFFF", linetype = 0)
  )
}

# set legend text
SetLegendTextGG <- function(x = 12, y = "bold") {
  return(theme(legend.text = element_text(size = x, face = y)))
}

# set legend point size
SetLegendPointsGG <- function(x = 6) {
  return(guides(colour = guide_legend(override.aes = list(size = x))))
}

# nice theme
theme_highlight <- function(){
  ggplot2::theme_classic() +
   ggplot2::theme(
      axis.text.x = element_text(face = "bold"),
      axis.title.x = element_text(face = "bold"),
      axis.text.y = element_text(angle = 0, face = "bold"),
      axis.title.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = 'bold'),
      strip.background = element_rect(fill = "#FFFFFF", linetype = 1, size = 0.001)
    )
}

# https://www.youtube.com/watch?v=UvRaNJC05Ec
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page
# https://www.r-bloggers.com/2020/05/vignette-generate-your-own-ggplot-theme-gallery
# https://github.com/tidyverse/ggplot2/blob/master/R/theme-defaults.r
# https://stackoverflow.com/questions/27689222/changing-fonts-for-graphs-in-r
# https://rpubs.com/Koundy/71792
theme_pub0 <- function(
  base_size = 12,
  base_family = "",
  base_line_size = base_size / 22,
  base_rect_size = base_size / 22
){
  half_line <- base_size / 2
  t <- theme(
    line = element_line(colour = "black", size = base_line_size, linetype = 1, lineend = "butt"),
    rect = element_rect(fill = "white", colour = "black", size = base_rect_size, linetype = 1),
    text = element_text(
      family = base_family, face = "plain",
      colour = "black", size = base_size, lineheight = 0.9,
      hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
      debug = FALSE
    ),
    axis.title.x = element_text(face = "bold", size = rel(0.8)),
    axis.title.y = element_text(face = "bold", size = rel(0.8)),
    axis.text.x = element_text(face = "plain", size = rel(0.8)),
    axis.text.y = element_text(face = "plain", size = rel(0.8)),
    legend.key.width = unit(0.2, 'cm'),
    legend.text = element_text(face = "plain", size = rel(0.6), hjust = 1),
    legend.title = element_text(face = "plain", size = rel(0.8)),
    plot.title = element_text(face = "bold", size = rel(1)),
    plot.subtitle = element_text(face = "plain", size = rel(0.8)),
    plot.caption = element_text(face = "plain", size = rel(0.6)),
    strip.text = element_blank(),
    strip.background = element_rect(fill = "#FFFFFF", linetype = 1, size = 0.001)
  )
  # ggplot_global$theme_all_null %+replace% t
}

# Quadrant counts
# inspired from
# https://rdrr.io/cran/ggpmisc/src/R/stat-quadrant-counts.R
# https://cran.r-project.org/web/packages/ggplot2/vignettes/extending-ggplot2.html
stat_quadrant_loc = function(qcut, center = FALSE) {
  qposition = list()
  for(i in 1:ncol(qcut)){
    for(j in 1:nrow(qcut)){
      y <- sapply(X = unlist(dimnames(qcut[j, i, drop = FALSE])),
        FUN = function(s){
          patty_i <- if(isTRUE(center)) "\\[|\\]|\\(|\\)" else if(grepl("\\)", s)) ",.*|\\[" else ".*,|\\]"
          z <- gsub(pattern = patty_i, replacement = "", x = s)
          if(isTRUE(center)) sapply(strsplit(z, ","), function(w) mean(as.numeric(w)) ) else as.numeric(z)
      })
      qposition[[paste0(j,i)]] <- y
    }
  }
  qposition <- data.frame(t(sapply(qposition, c)))
  if(!isTRUE(center)){
    qposition <- sapply(qposition, function(x) x + (max(scale(x)) * .2) )
    qposition <- data.frame(qposition)
  }; colnames(qposition) <- names(dimnames(qcut))
  qposition
}
stat_quadrant_fun = function(
  data, scales, params,
  xintercept = 0, yintercept = 0, type = "count", center = FALSE,
  loc_type = c("position", "component")
) {
  loc_type = match.arg(loc_type)
  qlimits = list(x = xintercept, y = yintercept)
  qcounts <- table(lapply(
    X = setNames(names(qlimits), names(qlimits)),
    FUN = function(x){
      y <- data[, x]
      if(min(y, na.rm = TRUE) == qlimits[[x]]){
        lo_i = paste0("[", min(y, na.rm = TRUE) + .001, ",", qlimits[[x]] + .002, ")")
        hi_i = paste0("[", qlimits[[x]] + .003, ",", round(max(y, na.rm = TRUE), 3), "]")
        ifelse(y > qlimits[[x]], hi_i, lo_i)
      }else{ Hmisc::cut2(x = y, cuts = qlimits[[x]]) }
  }))
  qposition <- if(loc_type == "position"){
    stat_quadrant_loc(qcut = qcounts, center = center)
  }else{
    colnames(qcounts) <- c("n", "p")[1:ncol(qcounts)]
    rownames(qcounts) <- c("n", "p")[1:nrow(qcounts)]
    data.frame(qcounts)[, -3]
  }
  qposition$count <- if(type == "percent"){
    scales::percent(c(qcounts / sum(qcounts)))
  }else{ c(qcounts) }
  qposition
}

StatQuadrant <- ggproto(
  "StatQuadrant", Stat,
  required_aes = c("x", "y"),
  compute_panel = stat_quadrant_fun,
  default_aes = ggplot2::aes(
    x = stat(x), y = stat(x),
    label = sprintf("%s", stat(count))
  )
)

stat_quadrant <- function (
  mapping = NULL, data = NULL, geom = "text", position = "identity",
  xintercept = 0, yintercept = 0, type = "count", center = FALSE,
  na.rm = FALSE, show.legend = NA,
  inherit.aes = TRUE, ...
) {
  ggplot2::layer(
    stat = StatQuadrant, data = data, mapping = mapping,
    geom = geom, position = position, show.legend = show.legend,
    params = list(
      na.rm = na.rm, xintercept = xintercept, yintercept = yintercept, type = type, center = center, ...
    )
  )
}

# `%>%` = dplyr::`%>%`
# centroid = data %>% dplyr::group_by(center_by) %>%
#   dplyr::summarize(x = type(x = x), y = type(x = y))
StatCentroid <- ggplot2::ggproto(
  "StatCentroid", ggplot2::Stat,
  required_aes = c("x", "y"),
  compute_panel = function(
    data, scales, params, center_by = "colour", type = median
  ) {
    data$center_by = data[, center_by]
    aggregate(cbind(x,y) ~ center_by, data = data, FUN = type)
  },
  default_aes = ggplot2::aes(x = stat(x), y = stat(x), label = stat(center_by))
)

stat_centroid <- function (
  mapping = NULL, data = NULL, geom = "text", position = "identity",
  center_by = "colour", type = median, na.rm = FALSE,
  show.legend = FALSE, inherit.aes = TRUE, ...
) {
  ggplot2::layer(
    stat = StatCentroid, data = data, mapping = mapping,
    geom = geom, position = position, show.legend = show.legend,
    params = list(
      na.rm = na.rm, center_by = center_by, type = type, ...
    )
  )
}

plot_add_quadrants = plots_add_quadrants = function(
  plot, limits = list(0, 0), ...
) {
  plot +
    geom_vline(xintercept = limits[[1]], linetype = "dashed", color = "gray50") +
    geom_hline(yintercept = limits[[2]], linetype = "dashed", color = "gray50") +
    stat_quadrant(xintercept = limits[[1]], yintercept = limits[[2]], ...)
}

plot_squared = plots_squared = function(
  gp, column_names = NULL, limit = NULL, verbose = FALSE
){
  column_map <- gsub("~", "", as.character(gp$mapping))
  if(verbose) cat("Mapping:", paste(column_map, collapse = ", "), "\n")
  if(is.null(column_names)) column_names <- column_map[1:2]
  if(verbose) cat("Using:", paste(column_names, collapse = ", "), "\n")
  dff <- gp$data[, column_names]
  if(all(sapply(dff, is.numeric))){
    if(!is.numeric(limit)) limit <- max(abs(dff), na.rm = TRUE)
    if(verbose) cat("Squared:", limit, "\n")
    gp <- gp + xlim(c(-limit, limit)) + ylim(c(-limit, limit))
  }
  gp
}

plot_add_densities = function(object, group){
  cowplot::plot_grid(
    margin_density(data=object$data, x=names(object$data)[1], group=group) , NULL,
    p, margin_density(data=object$data, x=names(object$data)[2], group=group) + coord_flip(),
    rel_widths = c(.8, .2), rel_heights = c(.2, .8)
  )
}
margin_density = function(data, x, group){
  ggplot(data, aes_string(x, fill=group)) +
    geom_density(alpha=.5) + theme_void() +
    theme(legend.position = "none")
}

plot_facet_wrap = function(p, facet, ...) {
  p + facet_wrap(facets = paste("~", facet), ...) +
    theme(
      strip.background = element_rect(fill = NA),
      strip.text = element_text(face = "bold")
    )
}

scatter_contour = function(
  data, axis_x, axis_y, axis_z = NULL, dp = TRUE
) {
  if(is.null(data$Density) && is.null(axis_z)){
    data$Density <- 0
    dps = if(isTRUE(dp)){
      rowSums(data[, c(axis_x, axis_y)] > 0) > 1
    }else{ rownames(data) }
    tmp = try(MASS_kde2d(
      x = data[dps, axis_x],
      y = data[dps, axis_y]
    ), silent = TRUE)
    if(class(tmp) != "try-error"){
      data[dps, ]$Density <- tmp; axis_z = "Density"
    }else{ data$Density = NULL }
  }
  aesy = if(!is.null(axis_z)){
    aes_string(x = axis_x, y = axis_y, color = axis_z)
  }else{ aes_string(x = axis_x, y = axis_y) }
  p <- ggplot(data = data, mapping = aesy) + geom_point(size = 0.5)
  if(!is.null(axis_z))
    p <- p + geom_density2d(data = data[dps, ], colour = "#c0c5ce")
  return(p)
}

plot_coord_radar = function (theta = "x", start = 0, direction = 1, ...) {
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") "y" else "x"
  ggproto(
    "CordRadar", CoordPolar, theta = theta, r = r, start = start,
    direction = sign(direction),
    is_linear = function(coord) TRUE, ...)
}

plot_radar = function(
  data, # feature, group, value
  enriched = TRUE,
  range_expand = 0.5,
  cex = 1,
  threshold = NULL
) {
  if(is.matrix(data)){
    data = as.data.frame(data)
    data = suppressWarnings(suppressMessages(
      reshape2::melt(cbind(data, rownames(data)))
    ))
  }
  colnames(data) <- c("feature", "group", "value")
  if(!isTRUE(enriched)){
    data[which(data[, 3] > 0), 3] <- 0; data[, 3] <- abs(data[, 3])
  }else{ data[which(data[, 3] < 0), 3] <- 0 }
  if(is.numeric(data[, 2]))
    data[, 2] <- as.character(data[, 2])
  if(!is.factor(data[, 2]))
    data[, 2] <- factor(data[, 2], gtools::mixedsort(unique(data[, 2])))

  # Make the circle not have a gap between the extreme points on X
  bridges <- data[data$group == levels(data[, 2])[1], ]
  bridges$group <- NA
  data <- rbind(data, bridges)

  d0 = data; d0[is.na(d0[, 3]), 3] <- 0
  range_0 = range_i = range(d0[, 3])
  range_i[1] <- range_i[1]-range_expand[1]
  if(!is.na(range_expand[2])) range_i[2] <- range_i[2]+range_expand[2]

  aesy = aes_string(x = "group", y = "value", color = "feature", group = "feature")
  groups <- setNames(1:nlevels(data[, 2]), levels(data[, 2]))
  p <- ggplot(data = data, mapping = aesy) +
    geom_path(color = "#BEBEBE", size = cex*1.0) +
    geom_line(color = "#ededed", size = cex*0.5) +
    geom_polygon(data = d0, color = "#919191", size = cex*0.1, fill = NA) +
    geom_point(size = cex*4.5) + ylim(range_i) +
    scale_x_discrete(expand = c(0,0), breaks = levels(data[, 2])) +
    labs(x = NULL, y = NULL) + theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.text.y = element_blank(), axis.ticks = element_blank(),
      axis.text.x = element_text(vjust = 3, size = rel(1.5)),
      panel.border = element_blank(), # panel.border = element_rect(fill = NA)
      panel.grid.major.x = element_line(size = rel(0.1), colour = "#dbdbdb"),
      panel.grid.major.y = element_line(size = rel(0.1), colour = "#dbdbdb")
    )
  dtext <- rbind(
    head(d0[which(d0[, 3] == range_0[1]), , drop = FALSE], 1),
    head(d0[which(d0[, 2] == levels(d0[, 2])[1]), ], 1),
    head(d0[which(d0[, 3] == range_0[2]), ], 1)
  ); dtext$value_str <- ""
  dtext[2:3, ]$value_str = paste(c("Min =", "Max ="), dtext[2:3, ]$value)
  if(is.numeric(threshold)){
    hline = data.frame(gr = seq(1,length(groups)+.9, by=0.1))
    hline$value = rep(threshold, nrow(hline))
    dtext[1, ]$value = threshold
    dtext[1, ]$value_str = paste("Threshold =", threshold)
    p <- p + geom_polygon( # can't be hline bc of plot_coord_radar
      data = hline, mapping = aes(x = gr, y = value), inherit.aes = FALSE,
      colour = "#919191", linetype = "dashed", fill = NA)
  }else{ dtext <- dtext[-1, ] }
  p <- p + geom_text(
    data = dtext, mapping = aes(x=group,y=value,label=value_str),
    nudge_y = -min(c(0.3, range_expand[1])), size = cex*2, inherit.aes = FALSE
  )
  p <- p + plot_coord_radar()
  return(p)
}

# devtools::install_git('https://bitbucket.org/nicholasehamilton/ggtern')
scatter_compare <- function(
  data,
  groups = NULL,
  label = NULL # Number of features to show or a vector of features
){
  if(is.null(groups)) groups = head(colnames(data), 3)
  if(!is.null(label)){
    maxs = if(is.numeric(label)){
      unique(unlist(lapply(groups, function(x){
        rownames(data[head(order(data[, x]), label), ])
      })))
    }else{ label }
    datatr <- data[maxs, ]
    datatr$gene_name <- rownames(datatr)
    # datatr$gene_name[!rownames(data) %in% maxs] <- ""
  }

  p <- if(!is.na(groups[3])){
    fun = ggtern::ggtern # ggtern::coord_tern(expand = FALSE)
    aesy = aes_string(x = groups[1], y = groups[2], z = groups[3])
    fun(data = data, mapping = aesy) +
      ggtern::stat_density_tern(geom = 'polygon',
        aes(fill  = ..level.., alpha = ..level..)) +
      guides(color = "none", alpha = "none") +
      geom_point() + scale_fill_gradientn(colours = viridis::viridis(25)) +
      theme(
        tern.axis.arrow.show = TRUE, legend.position = "bottom",
        plot.margin = margin(-0.5, -0.5, -0.5, -0.5, "cm"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey95"),
        panel.grid = element_line(colour = "grey97", size = rel(0.7)),
        strip.background = element_rect(fill = "#FFFFFF")
      )
  }else{
    scatter_contour(data, groups[1], groups[2])
  }

  if(!is.null(label)){
    # position = position_dodge(0.9), # nudge not approved by ggtern
    # fun_text = if(!is.na(groups[3])) ggtern::geom_text_viewport else ggplot2::geom_text
    # if(requireNamespace("ggrepel", quietly = TRUE)) fun_text = ggrepel::geom_text_repel
    p <- p + ggplot2::geom_text(data = datatr, mapping = aes(label = gene_name), color = 'black')
    # p <- p + ggtern::annotate(geom  = 'text', # same as ggplot2::geom_text
    #   x = datatr[, groups[1]], y = datatr[, groups[2]],
    #   z = datatr[, groups[3]], label = datatr$gene_name)
  }
  return(p)
}
# set.seed(1)
# ddf = data.frame(Ab = runif(100), An = runif(100), Or = runif(100))
# scatter_compare(ddf)

### Utilities ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
scatter_summary_df = function(
  data, features, threshold = 0, group = NULL
) {
  if(is.null(group)){ data$Data = "Data"; group = "Data" }
  groups = levels(factor(data[, group]))
  features_stats = lapply(
    X = 1:ncol(features),
    FUN = function(i){
      z <- c(as.character(features[1, i]), as.character(features[2, i]))
      pdata_i = data[, z]
      colnames(pdata_i) <- c("x", "y")
      lapply(
        X = setNames(nm = groups),
        FUN = function(x){
          y <- stat_quadrant_fun(
            pdata_i[data[, group] == x, ],
            xintercept = threshold, yintercept = threshold,
            type = "percent", loc_type = "component"
          ); y[, 1] <- paste0(features[1, i], y[, 1])
          y[, 2] <- paste0(features[2, i], y[, 2]); y
      })
    }); suppressMessages(reshape2::melt(features_stats))
}

# create numeric breaks for ggplot legends
make_breaks <- function( # pretty brakes!!!
  x,
  n = 10,
  push = c(min = -0, max = 0),
  roundy = 2
){
  if(length(push) == 1 && is.null(names(push))) push <- c(min = -push, max = push)
  punames <- c('min', 'max')
  push <- push[punames]; names(push) <- punames; push[is.na(push)] <- 0
  lb <- suppressWarnings(min(x, na.rm = TRUE)) + push['min']
  hb <- suppressWarnings(max(x, na.rm = TRUE)) + push['max']
  if(is.infinite(lb)) lb <- 0
  if(is.infinite(hb)) hb <- 1
  breaky <- round(seq(lb, hb, length.out = n), roundy)
  breaky[1] <- lb
  breaky[length(breaky)] <- hb
  return(breaky)
}

# how to fit a number of objects in a 2D dimention
make_grid <- function(x, n_col = NULL, verbose = FALSE){
  xlen <- if(length(x) > 1) length(x) else if(is.numeric(x)) x
  if(isTRUE(xlen == 0) || is.null(xlen)) return(c(1, 1))
  xlen <- length(seq(xlen))
  if(verbose) cat('Length:', xlen, '\n')
  if(xlen %in% 1) return(c(1, 1))
  if (!is.null(x = n_col)) {
    if(verbose) cat('Adjust for n_col:', n_col, '\n')
    n_row <- floor(x = xlen / n_col - 1e-5) + 1
  }else if(!sqrt(xlen) %% 1){
    n_row <- n_col <- sqrt(xlen)
  }else{
    n_col <- round(sqrt(xlen) + .5)
    n_row <- max(c(round(xlen / n_col), 1))
  }
  if(n_row * n_col < xlen) n_row <- n_row + 1
  if(verbose) cat('n_row:', n_row, '\nn_col', n_col, '\n')
  return(c(n_row, n_col)) # columns increase faster than rows
}

plot_size = function(x, n_base = 7, n_max = 25)
  sapply(x, function(y) min(c(n_max, y * n_base)) )

make_title <- function(x){
  y <- gsub("orig|\\.", "_", casefold(x, upper = TRUE), ignore.case = TRUE)
  y <- gsub("_", " ", gsub("_{2,}", "_", y))
  y <- gsub(" {1,}", " ", y)
  y <- gsub(" $|^ ", "", y)
  gsub("_$|^_", "", y)
}
