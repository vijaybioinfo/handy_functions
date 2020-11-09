#!/usr/bin/R

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
# ' p <- violin(dat = iris, xax = "Species", yax = "Sepal.Length", colour_by = 'Sepal.Length')
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
  vsetting <- lapply(names(vsetting_master), function(x) ifelse(x %in% names(vsetting), vsetting[[x]], vsetting_master[[x]]) )
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
  dat$pct <- dat$pct50 <- sapply(as.character(dat[, xax]), function(X) return(pct[X]) )
  means <- tapply(dat[, yax], dat[, xax], mean, na.rm = TRUE)
  dat$means <- dat$mean <- sapply(as.character(dat[, xax]), function(X) return(means[X]) )

  if(is.null(couls)) couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
  brks <- NULL
  master_brks = list( # create a master for breaks in mean and percentage
    pct = list(filname = '%+', brks = c(0, 10, 20, 40, 60, 80, 100)),
    mean = list(filname = 'Mean', make_breaks(dat[, yax])),
    means = list(filname = 'Mean', brks = c(0, 2, 4, 6, 8, 10, 12)),
    pct50 = list(filname = '%+cells', brks = c(0, 5, 10, 20, 30, 40, 50))
  ) # substitute the master with the given breaks
  if(is.list(colour_by)){
    if(is.null(names(colour_by))) warning("Breaks must be a named list of list.")
    master_brks[[names(colour_by)]] <- colour_by; colour_by <- names(colour_by)
  }
  if(isTRUE(colour_by %in% names(master_brks))){
    filname <- master_brks[[colour_by]][[1]]
    brks <- master_brks[[colour_by]][[2]]
  }; if(colour_by == "pct50") colour_by <- "pct"

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
    }; couls <- colorRampPalette(couls, space = 'Lab')(length(brks))
    lbrks <- as.character(round(brks, 2))
    if(any(grepl("\\-", lbrks))) lbrks <- ifelse(!grepl("\\-", lbrks), paste0("+", lbrks), lbrks)
    lbrks[grepl("\\+0$", lbrks)] <- "  0"
    arg_list <- list(name = filname, colours = couls, na.value = 'transparent',
      breaks = brks, labels = lbrks, limits = c(min(brks), max(brks)), guide = "colorbar")
    dp <- dp + guides(color = guide_colorbar(barwidth = 0.8, barheight = 6)) +
        do.call(scale_color_gradientn, arg_list) +
        do.call(scale_fill_gradientn, arg_list)
  }else if(length(unique(dat[, xax])) < 10) dp <- dp + scale_fill_brewer(palette = "Set1")
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
  legend <- cowplot::get_legend(x)
  graphics.off(); #tmp <- file.remove('Rplots.pdf')
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
  p <- getlegend(x[[1]]); x <- lapply(x, function(y) y + NoLegend() ); x$legend <- p
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
  type = c("NMF", "pheatmap"),
  ... # arguments for NMF::aheatmap or pheatmap
){
  type <- match.arg(type)
  if(casefold(class(object)) == 'seurat'){
    if(verbose) cat("Seurat object\n")
    if(is.null(annoc)) annoc <- object@meta.data
  }else if(is.null(annoc)){ stop("Annotation not provided") }
  if(is.null(rnames)) rnames <- rownames(object)
  if(is.null(cnames)) cnames <- colnames(object)
  cnames <- show_found(cnames, rownames(annoc), verbose = FALSE)
  cnames <- show_found(cnames, colnames(object), verbose = FALSE)
  if(!any(categorical_col %in% colnames(annoc))){
    categorical_col <- extract_grouping_cols(annoc, ckeep = categorical_col, v = verbose)
  }; categorical_col <- categorical_col[categorical_col %in% colnames(annoc)]
  if(is.null(orderby)) orderby <- categorical_col[1]

  if(isTRUE(sample_it[1]) || is.character(sample_it)){
    nsample <- if(length(sample_it) > 1) as.numeric(sample_it[2]) else NULL
    sample_it <- if(isTRUE(sample_it[1])) categorical_col[1] else sample_it[1]
    cnames <- sample_grp(annot = annoc, cname = sample_it, maxln = nsample, v = verbose)
  }

  if(verbose){
    cat("\nCategories:", show_commas(categorical_col), "\n")
    cat("Features:", show_commas(rnames), "\n")
    cat("Cells:", show_commas(cnames), "\n")
    cat("Variables:", show_commas(colnames(annoc)), "\n\n")
  }

  tvar <- unique(c(orderby, categorical_col, use_mean))
  annoc <- annoc[cnames, tvar[tvar %in% colnames(annoc)], drop = FALSE]
  object <- object[, cnames]

  if(is.character(use_mean) || isTRUE(use_mean)){
    if(casefold(class(object)) == 'seurat'){
      if(verbose) cat("Taking only the matrix from Seurat object\n")
      object <- expm1(object@assays$RNA@data)
    }
    use_mean <- if(isTRUE(use_mean[1])) head(colnames(annoc)[sapply(annoc, is.character)], 1) else use_mean
    annoc$index123 <- do.call(paste, c(annoc[, use_mean, drop = FALSE], sep = "_"))
    if(verbose) cat("Means:", paste0(use_mean, collapse = ", "), "\n")
    means <- stats_summary_table(
      mat = object,
      groups = make_list(annoc, colname = 'index123', grouping = TRUE),
      moments = "mn",
      v = verbose
    )
    colnames(means) <- gsub("_mean", "", colnames(means))
    annoc <- annoc[!duplicated(annoc$index123), ]
    annoc <- annoc[, !colnames(annoc) %in% "index123", drop = FALSE]
    tvar <- sapply(annoc, function(x) length(unique(x)) )
    if(verbose){ print(tvar); str(annoc) }
    annoc <- annoc[, tvar > 1, drop = FALSE]
    rownames(annoc) <- do.call(paste, c(annoc[, use_mean, drop = FALSE], sep = "_"))
    means <- means[, rownames(annoc)]
    head(annoc); head(means)
    object <- log2(as.matrix(means) + 1)
  }
  if(!is.null(orderby)) annoc[, orderby[1]] <- droplevels(annoc[, orderby[1]])

  if(!is.numeric(annoc[, orderby[1]]) && (length(orderby) > 1)){
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
    str(tvar)
    annoc <- annoc[unname(unlist(tvar)), ]
  }else{
    annoc <- annoc[order(annoc[, orderby]), , drop = FALSE] # will follow the factors order
  }

  tvar <- sapply(annoc, is.character) | sapply(annoc, is.factor)
  anncolist <- lapply(annoc[, tvar, drop = FALSE], function(x){
    v2cols(select = x, sour = couls, v = FALSE)
  })
  row_params <- c('annRow', 'annotation_row')
  if(any(row_params %in% names(list(...)))){
    if(verbose) cat("Adding row annotation\n")
    tvar <- data.frame(
      list(...)[[which(names(list(...)) %in% row_params)]],
      check.names = FALSE, stringsAsFactors = FALSE
    )
    anncolist <- c(anncolist,
      lapply(tvar[, sapply(tvar, is.character), drop = FALSE], function(x){
        v2cols(select = x, sour = couls, v = verbose)
      })
    )
  }
  if(verbose) str(anncolist)

  object <- object[rnames, rownames(annoc)]

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
    mat_to_plot <- t(scale(t(mat_to_plot)))
  }

  nnpcs <- min(c(3, dim(mat_to_plot) - 1))
  thesegenes <- if((isTRUE(feature_order) || as.character(feature_order) == "pca") && (casefold(class(object)) == 'seurat')){
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
    if(verbose) cat("Ordering based on hierarchical ordering\n"); TRUE
  }else{ NA }

  if(verbose) cat("Creating heatmap plot\n")
  thesegenes <- c(thesegenes[thesegenes %in% rnames], rnames[!rnames %in% thesegenes])
  mat_to_plot <- as.matrix(mat_to_plot[thesegenes, rownames(annoc)])

  if(is.null(topz)) topz <- max(c(min(abs(c(range(mat_to_plot), 2))), 1))
  if(verbose) cat("Range:", range(mat_to_plot), "\n")
  mat_to_plot[mat_to_plot > topz] <- topz; mat_to_plot[mat_to_plot < (-topz)] <- -topz;
  if(verbose) str(mat_to_plot)

  mypalette <- colorRampPalette(colors = rev(hcouls), space = 'Lab')
  palettebreaks <- seq(from = -topz, to = topz, by = 0.1)

  if(type == "NMF"){
    cat("Using NMF\n")
    NMF::aheatmap(
      mat_to_plot,
      scale = 'none', Rowv = feature_order, Colv = NA,
      annCol = annoc, annColors = anncolist,
      col = mypalette(length(palettebreaks) - 1),#7)),
      ...
    )
  }else{
    cat("Using pheatmap\n")
    source('/home/ciro/scripts/handy_functions/devel/pheatmapCorrection.R')
    if(is.na(feature_order) || is.character(feature_order)) feature_order <- FALSE
    pheatmap(
      mat = mat_to_plot,
      scale = "none", cluster_rows = feature_order, cluster_cols = FALSE,
      annotation_col = annoc, annotation_colors = anncolist,
      color = mypalette(length(palettebreaks) - 1),
      breaks = palettebreaks,
      border_color = NA,
      drop_levels = TRUE,
      ...
    )
  }
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
    facet_wrap(~in_group, ncol = make_grid(unique(ddf2$in_group))[2])
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
      facet_wrap(~in_group, ncol = make_grid(unique(ddf2$in_group))[2]) +
      xlim(c(2, 4))
  }
  p <- p + theme_minimal() + labs(x=NULL)
  if(isTRUE(return_table)) return(list(plot = p, table = prop_table))
  return(p)
}

### Utilities ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create numeric breaks for ggplot legends
make_breaks <- function(
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
make_grid <- function(x, ncol = NULL, verbose = FALSE){
  vecL <- if(is.numeric(x)) x else if(is.character(x)) as.numeric(x) else length(x)
  vecL <- length(seq(vecL))
  if(verbose) cat('Length:', vecL, '\n')
  if (!is.null(x = ncol)) {
    if(verbose) cat('Adjust for ncol:', ncol, '\n')
    if (vecL == 1) ncol <- 1
    # if (vecL > 6) ncol <- 3
    # if (vecL > 9) ncol <- 4
    nRow <- floor(x = vecL / ncol - 1e-5) + 1
  }else{
    if(vecL %in% 1:2) ncol <- 1 else if(!is.numeric(ncol)) ncol <- round(sqrt(vecL) + .5)
    if(vecL == 1) nRow <- 1 else nRow <- round(vecL / ncol)
    if(nRow * ncol < vecL) nRow <- nRow + 1
    if(!sqrt(vecL) %% 1) nRow <- ncol <- sqrt(vecL) # if sqrt is integer
  }
  if(verbose) cat('nRow:', nRow, '\nncol', ncol, '\n')
  return(c(nRow, ncol, ifelse(vecL <= 2,  8, 4)))
}

make_title <- function(x){
  y <- gsub("orig|\\.", "_", casefold(x, upper = TRUE), ignore.case = TRUE)
  y <- gsub("_", " ", gsub("_{2,}", "_", y))
  y <- gsub(" {1,}", " ", y)
  y <- gsub(" $|^ ", "", y)
  gsub("_$|^_", "", y)
}
