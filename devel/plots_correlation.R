#!/usr/bin/R

correlation_heatmap <- function(
  expdata,
  cortype = 'spearman',
  breaks = c(-.7, .7),
  title = 'Correlation',
  caption = NULL,
  corthr = 0, # correlation threshold
  ntop = 0.1, # top correlating features
  unsigned = FALSE, # plot absolute correlations two
  return_mat = FALSE,
  return_plot = FALSE,
  features_show = NULL, # show only these names in the plot
  features_cols = NULL, # NA to turn off in any case
  feature_position = FALSE, # show the position of the feature
  clust_cols = NULL, # cluster colours
  # Cosmetic
  cexs = NULL,
  axes_order = TRUE,
  verbose = TRUE
){
  if(verbose) cat('--- Correlation ---\n')
  if(verbose) cat(class(expdata), "object given\n")
  if(casefold(class(expdata)) %in% c("saver", "list")){
    datasaver <- expdata; expdata <- datasaver$estimate
  }
  matdata <- as.matrix(expdata)
  if(verbose) cat('Features:', nrow(matdata), '=>', commas(rownames(matdata)), '\n')
  if(verbose) cat('Samples:', ncol(matdata), '=>', commas(colnames(matdata)), '\n')

  if(verbose) cat('Correlation:', cortype, '\n')
  if(exists('datasaver')){
    saver1 <- lapply(datasaver[1:2], function(x) x[rownames(matdata), colnames(matdata)] )
    matdata <- SAVER::cor.genes(saver1)
  }else{
    matdata <- cor(t(matdata), method = cortype)
  }
  matdata[is.na(matdata)] <- 0
  diag(matdata) <- 0
  matdata.abs <- abs(matdata)
  maxabs <- max(matdata.abs, na.rm = T)
  min_all <- min(matdata, na.rm = T)
  max_all <- max(matdata, na.rm = T)

  # taking the highly correlating features only
  caption <- if(is.null(caption)) cortype else caption
  caption <- paste0(caption, '\n')
  if(corthr[1] != 0){
    if(is.character(corthr[1])){
      if(!sum(is.na(features_cols)) && is.null(features_cols)){ # when there NAs and is null?
        features_cols <- rep('black', nrow(matdata)) # specify colours
        names(features_cols) <- rownames(matdata)
        features_cols[names(features_cols) %in% corthr] <- 'red'
      }else if(sum(is.na(features_cols))){ features_cols <- NULL }
      if(verbose) cat(corthr, '\n')
      tvar <- matrixStats::rowMaxs(abs(matdata[corthr, , drop = F]))
      corthr <- min(tvar) - 0.0005
      ntop <- 1
    }
    if(ntop <= 0 || (ntop >= ncol(matdata) - 1)) ntop <- 0.1
    if(ntop < 1){
      ntop <- round(ntop*ncol(matdata))
    }
    keep = (colSums(abs(matdata) >= corthr) >= ntop)
    tvar <- colSums(abs(matdata) >= corthr)
    tmp <- paste('Filtering:', round(corthr, 3), 'cor in >=', ntop, 'features')
    if(verbose) cat(tmp, '- ')
    if(sum(keep)){
      if(verbose) cat(nrow(matdata), 'down to', sum(keep), '\n')
      matdata <- matdata[keep, keep]
      caption <- paste(caption, tmp, '\n')
    }else{
      if(verbose) cat('threshold out of bound:', round(maxabs, 2), '/', max(tvar[tvar!=0]), '\n')
    }
  }

  # breaks
  if(breaks[1] == 'seq'){
    if(verbose) cat('Sequence breaks\n')
    palette.breaksabs <- seq(0, maxabs, 0.15)
    palette.breaks <- seq(min_all, max_all, 0.15)
  }else if(breaks[1] == 'quantile'){
    if(verbose) cat('Quantile breaks\n')
    palette.breaksabs  <- quantile_breaks(matdata.abs, n = 11)
    palette.breaks <- quantile_breaks(matdata, n = 11)
  }else if(length(breaks) == 2 && sum(breaks <= 1 & breaks >= -1) == 2){
    if(verbose){
      cat('Rank breaks\n')
      cat(' Upper bound:', max_all,'\n')
      cat(' Lower bound:', min_all,'\n')
    }
    if(breaks[1] == breaks[2] && sum(breaks == 1 & breaks == -1) != 2) breaks <- c(breaks[1], 1)
    breaks <- sort(abs(breaks))
    palette.breaksabs <- unique(c(seq(0, breaks[1], 0.05),
                                  seq(breaks[1], breaks[2], 0.15),
                                  seq(breaks[2], 1, 0.05)))
    palette.breaks <- c(seq(-breaks[2], -breaks[1], 0.15),
                                  seq(-breaks[1], breaks[1], 0.05),
                                  seq(breaks[1], breaks[2], 0.15))
    palette.breaksabs <- sort(unique(round(palette.breaksabs,4)))
    palette.breaks <- sort(unique(round(palette.breaks,4)))
  }else if(is.numeric(breaks) && sum(breaks <= 1 & breaks >= -1)){
    if(verbose){
      cat('Cutoff breaks\n')
      cat(' Upper bound:', max_all,'\n')
      cat(' Lower bound:', min_all,'\n')
    }
    breaks <- abs(breaks)
    palette.breaksabs <- unique(c(seq(0, breaks, 0.05),
                                  seq(breaks, 1, 0.15)))
    palette.breaks <- c(-1, seq(-1, -breaks, 0.15),
                                  seq(-breaks, breaks, 0.05),
                                  seq(breaks, 1, 0.15), 1)
    palette.breaks <- sort(unique(round(palette.breaks,4)))
  }else{
    if(verbose) cat('Sequence breaks\n')
    palette.breaksabs <- seq(0, maxabs, 0.05)
    palette.breaks <- seq(min_all, max_all, 0.05)
  }
  rgb.paletteabs <-  colorRampPalette(c("ghostwhite", "#8B0000"), space = "Lab")
  rgb.palette <-  colorRampPalette(c("#00008B", "ghostwhite", "#8B0000"), space = "Lab")

  # features to show
  if(!is.null(features_show[1])){
    if(verbose) cat('Showing only ')
    tvar <- rownames(matdata)
    features_show <- features_show[features_show != '']
    tvar[!tvar %in% features_show] <- ''
    features_show <- tvar
    if(verbose) cat(sum(features_show != ''), 'features\n')
  }else{ features_show <- rownames(matdata) }

  # features to colour
  if(!is.null(features_cols) && !is.numeric(features_cols)){
    # can't be passed without being affected by row reordering
    if(verbose) cat('Colouring only ')
    features_cols <- features_cols[rownames(matdata)]
    names(features_cols) <- rownames(matdata)
    if(!is.null(features_show[1])) tvar <- sum(features_show != '') else tvar <- length(features_cols)
    if(verbose) cat(sum(!is.na(features_cols)), 'features\n')
  }

  if(is.null(cexs) || (cexs > 0.2 && ncol(matdata) < 100)) cexs <- .1 + 1/(3.5 * log10(nrow(matdata)))
  if(cexs < 0.5 && ncol(matdata) < 70) cexs <- 0.6
  if(verbose) cat('Text size: ', round(cexs, 3), '\n')

  # identifying lcusters
  if(is.numeric(clust_cols)){
    if(verbose) cat("Detecting clusters\n")
    hr <- hclust(dist(matdata))
    if(verbose) print(hr)
    if(clust_cols %% 1 != 0){
      if(verbose) cat("Defining clusters at height", clust_cols, "\n")
      mycl <- cutree(hr, h = clust_cols)#max(hr$height / 2))
    }else{
      if(verbose) cat("Defining", clust_cols, "clusters\n")
      mycl <- cutree(hr, k = clust_cols)
    }
    clust_cols <- get_cols(mycl)[mycl]; names(clust_cols) <- names(mycl)
    if(is.null(features_cols)) features_cols <- clust_cols
    if(verbose)str(clust_cols)
  }

  tmp <- paste(nrow(matdata), 'features')
  if(verbose) cat("Building plot", tmp, '\n')
  caption <- paste0(caption, tmp)
  if(sum(dim(matdata) >= 2) == 2){
    if(return_plot || feature_position) pdf(file = NULL)
    # source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/heatmap2Correction.R')
    # d <- heatmap.2c(
    d <- gplots::heatmap.2(
      x = matdata, scale = "none", col = rgb.palette, main = title,
      keysize = 1, key.title = NA, #lhei = c(1, 7), lwid = c(1, 1),
      Rowv = axes_order, Colv = axes_order, dendrogram = c("row", "none")[1], trace = "none",
      cexRow = cexs, cexCol = cexs,
      breaks = palette.breaks,
      symm = TRUE, revC = axes_order, add = FALSE,
      labRow = features_show, colRow = features_cols,
      RowSideColors = clust_cols
    )
    if(unsigned && !return_plot){
      heatmap.2c(matdata.abs, scale = "none",col = rgb.paletteabs, main = title,
                keysize = 1, key.title = NA, #lhei = c(1, 7), lwid = c(4, 7),
                Rowv = axes_order, Colv = axes_order, dendrogram = "row", trace = "none", cexRow = cexs, cexCol = cexs,
                breaks = palette.breaksabs,
                symm=T, revC=axes_order,
                labRow = features_show, colRow = features_cols)
      if(caption[1] != '') legend('top', legend = caption, 8, cex=.9, bty='n', horiz=T) # Groups labels
    }
    if(feature_position || !is.null(features_cols)){
      matdata <- matdata[rev(colnames(d$carpet)), rev(colnames(d$carpet))]
      tvar <- rownames(matdata); features_show <- features_show[features_show != '']; tvar[!tvar %in% features_show] <- ''
      features_show <- tvar # reorder features_show according to rownames
      if(feature_position){
        colnames(matdata) <- paste0(1:ncol(matdata), " - ", colnames(matdata))
        features_show <- paste0(1:length(features_show), " - ", features_show)
        features_show <- sub(".*- $", "", features_show) # eliminate empty ones
        if(all(features_show != '')){
          tvar <- round(quantile(1:length(features_show), c(.25, .5, .75))) # get Qs
          features_show[tvar] <- paste(features_show[tvar], ' <-- Q', c(25, 50, 75))
          if(is.null(features_cols)){
            features_cols[tvar] <- 'red' # mark the quantiles
            features_cols[is.na(features_cols)] <- features_cols[nrow(matdata)] <- 'black' # complete
            names(features_cols) <- rownames(matdata)
          }
        }
      }
      features_cols <- rev(features_cols[rownames(matdata)])
      if(!is.null(clust_cols)) clust_cols <- clust_cols[names(features_cols)]
      # print(features_cols)
      if(!is.null(features_cols)) features_show <- rev(features_show) # colour doesnt reverse the names
    }
    if(return_plot){
      if(verbose) cat("Returning plot\n")
      dev.off(); return(d)
    }else if(feature_position){
      if(verbose) cat("Plotting\n")
      dev.off(); eval(d$call)
    }
    legend('top', legend = caption, 8, cex = .9, bty = 'n', horiz = TRUE)
    rdf <- data.frame(feature_name = rev(colnames(d$carpet)), stringsAsFactors = FALSE)
    if(exists("mycl")) rdf$cluster <- mycl[rdf$feature_name]
    if(return_mat){
      if(verbose) cat("Returning ordered features, correlation matrix and plot\n")
      return( list(features = rdf, mat_cor = matdata, mat_heat = d$carpet) )
    }else{
      if(verbose) cat("Returning ordered features only\n")
      return( rdf )
    }
  }else{ warning('Nothing to plot'); return( 0 ) }
}

## Utilities
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
commas <- function(x) paste0(head(x), collapse = ", ")
gg_color_hue <- function(n) {
  if(length(n) > 1){
    mynames <- n
    n <- length(n)
  }else{ mynames <- NULL }
  hues = seq(15, 375, length = n + 1)
  tmp <- hcl(h = hues, l = 65, c = 100)[1:n]
  if(!is.null(mynames)) names(tmp) <- mynames
  return(tmp)
}
get_cols <- function(x){
  y <- unique(x)
  z <- gg_color_hue(length(y))
  names(z) <- y
  z
}
