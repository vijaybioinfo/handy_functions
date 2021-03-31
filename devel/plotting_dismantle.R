#!/usr/bin/R

######################
# Plotting functions #
######################

# suppressPackageStartupMessages(library("ggExtra"))
# Highly customed plot functions
# suppressPackageStartupMessages(library("pack"))
# depends <- c("plotly", "ggrepel", "gtable" ,"doParallel", "gridExtra", "plyr",
#   "dplyr", "cowplot", "ggpubr", "GGally", "grid", "Hmisc",
#   "ggExtra", "ks", "ggsignif")
# load_packs(depends, v = FALSE)
suppressPackageStartupMessages(library(ggplot2))
#### General features #### -----------------------------------------------------
# Transform normal plot to grob
gwrapp <- function(x) {
    arbol <- grid.grabExpr( grid.echo(x) )
    y <- unit(1, 'null')
    gtable_col(NULL,list(arbol), y, y)
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

set_big_text <- function(x, sizey = 24, facey = "plain"){
  x <- x + theme(
    plot.title = element_text(size = sizey * 1.2, face = facey),
    text = element_text(size = sizey, face = facey),
    legend.key.size = unit(1, "cm")
  ) #+
  # guides(colour = guide_legend(override.aes = list(size = sizey)))
  return(x)
}

# limit lines
limit_lines <- function(x, y){
  up <- rep(-log10(y), 2)
  bottom <- c(x, x)
  top <- c(-log10(y), Inf)
  sides <- c(x, Inf)
  geom_line(data = data.frame(x = -rev(sides), y = up), aes(x = x, y = y), linetype = "dotted") +
  geom_line(data = data.frame(x = sides, y = up), aes(x = x, y = y), linetype = "dotted") +
  geom_line(data = data.frame(x = bottom, y = top), aes(x = x, y = y), linetype = "dotted") +
  geom_line(data = data.frame(x = -bottom, y = top), aes(x = x, y = y), linetype = "dotted")
}

# Add percentage per quadrant
pct_per_quadrant <- function(
  gp,
  cnames = NULL,
  limits = c(0, 0),
  dist_push = 0,
  simetric = TRUE
){
  if(class(gp)[1] == "gg"){
    if(is.null(cnames)){
      cnames <- gsub("~", "", as.character(gp$mapping))[1:2]
    }
    x <- gp$data
  }else{
    x = gp
  }
  if(length(limits) == 1) limits <- c(limits, limits)

  qposition <- t(combn(sapply(x[, cnames], range), 2)[, c(6:4, 2)])
  # limits[3] <- limits[[1]] * ifelse(any(qposition > 0), -1, 1)
  limits[3] <- limits[[1]] * ifelse(all(qposition > 0), 1, -1)
  qposition <- data.frame(qposition)# + scale(qposition))
  rownames(qposition) <- c('lt',  'rt', 'rb', 'lb')
  colnames(qposition) <- c('x', 'y')
  print(qposition)
  qposition$d <- c(
    np = sum(x[, cnames[1]] <= limits[3] & x[, cnames[2]] > limits[2], na.rm = TRUE),
    pp = sum(x[, cnames[1]] > limits[1] & x[, cnames[2]] > limits[2], na.rm = TRUE),
    pn = sum(x[, cnames[1]] > limits[1] & x[, cnames[2]] <= limits[3], na.rm = TRUE),
    nn = sum(x[, cnames[1]] <= limits[3] & x[, cnames[2]] <= limits[3], na.rm = TRUE)
  )
  qposition$d <- round(qposition$d / nrow(x) * 100, 2)
  qposition$percent <- paste0(qposition$d, "%")
  if(dist_push != 0){
    outergrid <- c(
      xmin = min(x[, cnames[1]]) - dist_push, xmax = max(x[, cnames[1]]) + dist_push,
      ymin = min(x[, cnames[2]]) - dist_push, ymax = max(x[, cnames[2]]) + dist_push
    )
  }else if(isTRUE(simetric)){
    mymax <- max(abs(x[, cnames]), na.rm = TRUE)
    outergrid <- c(
      xmin = -mymax - dist_push, xmax = mymax + dist_push,
      ymin = -mymax - dist_push, ymax = mymax + dist_push
    )
  }
  qposition['lt', 1:2] <- outergrid[c('xmin', 'ymax')]
  qposition['lb', 1:2] <- outergrid[c('xmin', 'ymin')]
  qposition['rt', 1:2] <- outergrid[c('xmax', 'ymax')]
  qposition['rb', 1:2] <- outergrid[c('xmax', 'ymin')]
  print(qposition); print(outergrid)
  if(class(gp)[1] == "gg"){
    if(dist_push > 0 && isTRUE(simetric)) gp <- gp + xlim(outergrid[1:2]) + ylim(outergrid[3:4])
    gp <- gp + geom_text(data = qposition, inherit.aes = FALSE, mapping = aes(x = x, y = y, label = percent))
  }else{
    gp <- qposition
  }
  return(gp)
}

blank_theme <- theme_minimal()+
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

shut_up <- theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.x = element_blank(),
  axis.text.y = element_blank(),
  # axis.ticks = element_blank(),
  plot.title = element_blank(),
  plot.subtitle = element_blank(),
  plot.caption = element_blank(),
  legend.position = "none",
  strip.text = element_blank(),
  strip.background = element_blank()#,
  # panel.background = element_blank(),
  # panel.grid.minor = element_blank(),
  # panel.grid.major = element_blank()
)

drawsq <- theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(size = 2)
  )

theme_Publication <- function(
  x,
  base_size = 14,
  base_family = "helvetica"
) {
  require(ggthemes)
  x <- x + #theme_foundation(base_size = base_size, base_family = base_family) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.2), hjust = 0.5),
      text = element_text(),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      axis.title = element_text(face = "bold", size = rel(1)),
      axis.title.y = element_text(angle = 90, vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.text = element_text(),
      axis.line = element_line(colour="black"),
      axis.ticks = element_line(),
      panel.grid.major = element_line(colour = "#f0f0f0"),
      panel.grid.minor = element_blank(),
      legend.key = element_rect(colour = NA),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.key.size = unit(0.2, "cm"),
      legend.margin = unit(0, "cm"),
      legend.title = element_text(face = "italic"),
      plot.margin = unit(c(10, 5 , 5, 5), "mm"),
      strip.background = element_rect(colour = "#f0f0f0", fill = "#f0f0f0"),
      strip.text = element_text(face="bold")
    )
  mapps <- gsub("factor\\(|\\)|~", "", as.character(x$mapping)[-(1:2)])
  if(length(x = mapps) && length(x$scales$scales)) {
    for(mapp in mapps) {
      type_scale <- ifelse(is.numeric(x$data[, mapp]), continuous_scale, discrete_scale)
      mypalette <- ifelse(is.numeric(x$data[, mapp]), "YlOrBr", "Paired")
      x <- x +
      type_scale(
        aesthetics = names(x$mapping)[grep(mapp, as.character(x$mapping))],
        scale_name = "Publication",
        palette = scales::manual_pal(values = RColorBrewer::brewer.pal(n = 9, name = mypalette)),
      )
    }
  }
  return(x)
}

# draw a vertical line outside (left) the plot
draw_vline <- function(gp){
  dff <- gp$data[, gsub("~", "", as.character(gp$mapping))[1:2]]
  ys <- range(dff[, 2])
  p1 <- gp + annotation_custom(grob = linesGrob(), xmin = -0.1, xmax = -0.1, ymin = ys[1]-3, ymax = ys[2]+4)
  # p1 <- gp + geom_segment(aes(x = -0.1, xend = -0.1, y = ys[1]-3, yend = ys[2]+4), arrow = arrow(length = unit(0.6, "cm")))
  pdf(file = NULL); gt <- ggplot_gtable(ggplot_build(p1)); dev.off()
  gt$layout$clip[grep("panel", gt$layout$name)] <- "off"
  # gt$layout[grep("ylab\\-l", gt$layout$name), c('t', 'b')] <- "off"
  grid.draw(gt)
}

## dot-lines plot
dot_lines <- function(
    cpmdata,
    metadata = NULL,
    rnames = NULL,
    gg = 'IL5',
    titl = NULL,
    orderby = NULL,
    maing = NA,
    covi = NA,
    facetting = TRUE,
    datatype = 'CTS',
    mysca = c('fixed', 'free_y'),
    legendt = c('posiv', 'mean', 'median'),
    rotx = FALSE,
    legendpos = 'right',
    return_plot = FALSE,
    v = FALSE
){
  set.seed(27)
  legendt <- legendt[1]
  mysca <- match.arg(mysca)
  # if(class(cpmdata) == 'seurat'){
  #   if(is.null(metadata)) metadata <- cpmdata@meta.data
  #   cpmdata <- cpmdata@data
  #   datatype <- 'Seurat Normalised'
  # }
  if(casefold(class(cpmdata)) == 'seurat'){
    if(is.null(metadata)) metadata <- cpmdata@meta.data
    if(v) cat('Seurat object\n')
    cpmdata <- GetAssayData(cpmdata, slot = "data")
    datatype <- 'Seurat'
  }
  if(is.null(orderby)){
    tvar <- colnames(metadata)
    orderby <- tvar[tvar %in% 'res.0.6']
    if(!length(orderby)) orderby <- head(tvar[grepl('res\\.', tvar)], 1)
    if(!length(orderby)) orderby <- head(tvar[grepl('orig\\.', tvar)], 1)
  }
  if(v) cat('Ordering by:', orderby, '\n')
  if(!class(cpmdata) %in% c('matrix', 'data.frame')) cpmdata <- as.matrix(cpmdata)
  if(!is.null(rnames)) metadata <- metadata[rnames, ]
  metadata <- as.data.frame(metadata)
  data.table::setorderv(metadata, orderby); rnames <- rownames(metadata) # order after selection
  cpmdata <- cpmdata[, rnames]
  if(v && !is.na(covi[1])) cat('Covariates:', commas(covi), '\n')
  if(sum(gg %in% colnames(metadata))){
    tvar <- gg[gg %in% colnames(metadata)]
    if(v) cat('Taking', commas(tvar), 'from annotation\n')
    tvar <- tvar[sapply(tvar, function(x) Hmisc::all.is.numeric(metadata[colnames(cpmdata), x]) )]
    tmp <- sapply(tvar, function(x) as.numeric(metadata[colnames(cpmdata), x]) )
    rownames(tmp) <- colnames(cpmdata); tmp <- t(tmp); rownames(tmp) <- tvar
    cpmdata <- rbind(cpmdata, tmp)
  }
  gg <- unique(getfound(x = gg, sour = rownames(cpmdata), element = 'Genes', v = v))
  for(g in gg){
    if(v) cat('.')
    allexp <- unlist(cpmdata[g, rnames])
    acats <- levels(metadata[, orderby])
    if(is.null(acats)) acats <- unique(metadata[, orderby])
    if(all(covi %in% colnames(metadata))){
      if(sum(is.na(maing))) maing <- acats
      tmeta <- metadata
      tmeta$covg <- do.call('paste', c(tmeta[covi], sep=''))
      setorder(tmeta, covg)
      acats <- unique(tmeta$covg) # saving all cats to preserve same number
      tmeta <- tmeta[tmeta[, orderby] %in% maing, ] # filter only for the main group
      dfexp <- data.frame(x = unique(tmeta$covg),
        mean = tapply(allexp[rownames(tmeta)], tmeta$covg, mean, na.rm = TRUE),
        median = tapply(allexp[rownames(tmeta)], tmeta$covg, median, na.rm = TRUE),
        posiv = tapply(allexp[rownames(tmeta)], tmeta$covg, function(x){
          round(sum(x > 0) / length(x) * 100, digits = 2) }),
        ncells = tapply(allexp[rownames(tmeta)], tmeta$covg, function(x) length(x))
      )
      if(sum(levels(dfexp$x) %in% acats) != length(acats)){ # make sure all levels will be there
        dfexp <- rbind(dfexp, data.frame(x = as.factor( acats[!acats%in%levels(dfexp$x)] )),
          mean = 0,
          posiv = 0,
          ncells = 0
        )
      }
      group <- sapply(acats, function(x) unique(tmeta[tmeta$covg == x, orderby]), simplify = F)
      dfexp <- cbind(group = unlist(group), dfexp[rep(names(group), sapply(group, length)), ])
    }else{
      dfexp <- data.frame(x = acats,
        mean = round(tapply(allexp, metadata[, orderby], mean, na.rm = TRUE), digits = 3),
        median = round(tapply(allexp, metadata[, orderby], median, na.rm = TRUE), digits = 3),
        posiv = tapply(allexp, metadata[, orderby], function(x){
          round(sum(x > 0) / length(x) * 100, digits = 2) }),
        ncells = tapply(allexp, metadata[, orderby], function(x) length(x))
      )
    }
    dfexp$x <- factor(dfexp$x, levels = gtools::mixedsort(levels(dfexp$x)))
    dfexp$gene <- g
    if(g == gg[1]) tdfexp <- dfexp else tdfexp <- rbind(tdfexp, dfexp)
  }; if(v) cat('\n')
  dfexp <- tdfexp; rm(tdfexp)
  dfexp$x <- factor(dfexp$x, levels = gtools::mixedsort(levels(dfexp$x)))
  dfexp$gene <- factor(dfexp$gene, levels = gg)

  if(v) cat('Building plot\n')
  vln <- ggplot(dfexp, aes_string(x = 'x', y = legendt, color = 'gene') ) +
    geom_point(alpha=0.7) + geom_line(aes(group = gene), alpha=0.7) +
    ylim(c(0, 100))
  # if(nlevels(dfexp$gene) == 1) vln <- vln + geom_bar(stat = 'identity', position = 'identity', fill = 'red', width = 0.2, alpha=0.3)
  if(rotx){ if(v) cat('Rotate\n'); rotx <- 90 }else rotx <- 0
  etext <- element_text(face = "bold", colour = "#990000")
  vln <- vln + theme_classic() +
    theme(
      axis.text.x = element_text(angle = rotx, face = "bold"),
      axis.title.x = etext,
      axis.text.y = etext,
      axis.title.y = etext,
      legend.position = ifelse(facetting, "none", legendpos),
      plot.title = element_text(face = "bold", hjust = 0.5),
      strip.text = element_text(face = 'bold'),
      strip.background = element_blank())

  # further modifications
  datatype <- casefold(paste(ifelse(isTRUE(legendt == 'posiv'), "%", 'mean'), datatype), upper = T)
  if(facetting){
    if(is.na(maing[1]) && length(gg) > 1){
      tvar <- fitgrid(1:length(gg))[2]
      if(v) cat('Faceting genes\n')
      vln <- vln + facet_wrap(~ gene, ncol = tvar, scale = mysca)
    }else if(is.na(covi[1]) && length(unique(dfexp$group)) > 1){
      tvar <- fitgrid(1:length(unique(dfexp$group)))[2]
      if(v) cat('Faceting groups\n')
      vln <- vln + facet_wrap(~ group, ncol = tvar, scale = mysca)
    }else if(length(unique(dfexp$group)) > 1){
      if(v) cat('Faceting groups ~ gene\n')
      # tvar <- fitgrid(1:(length(gg)+length(unique(dfexp$group))))[1]
      vln <- vln + facet_wrap(gene ~ group, ncol = length(unique(dfexp$group)), scale = mysca)
    }
  }
  noncero <- FALSE
  tmp <- ifelse(noncero, 'Expressing cells', 'All cells')
  if(!is.null(titl)) titl <- paste(titl, '-')
  titl <- paste(titl, ifelse(length(gg) > 1, 'Dot-line plots', gg))
  vln <- vln + labs(title = titl, subtitle = tmp, y = datatype, x = 'Identity')
  if(v) cat('Title:', titl, '\n')
  if(v) cat('Subtitle:', tmp, '\n')
  if(return_plot){ if(v) cat('Returning\n'); return(vln) }else{if(v) cat('Printing\n'); print(vln)}
}

# plot 'top' genes from a ranked list
getsiggenes <- function(
  plotdata,
  tannot,
  tcatg,
  tres,
  padjthr = 0.05,
  fcthr = 0.25,
  grpn = 'X',
  cellsn = NULL,
  ordr = 'log2FoldChange',
  dec = TRUE, ngenes = 12,
  datatype = "CTS"
){
  tres <- as.data.frame(tres); tres <- tres[order(tres[, ordr], decreasing = dec), ]
  entrez_to_plot <- head(tres[, 'gene'], ngenes)
  tmp <- getDEGenes(tres, pv = padjthr, fc = fcthr, gene_name = "gene")
  tmp <- length(entrez_to_plot[entrez_to_plot %in% tmp])
  # siggenes <- entrez_to_plot[entrez_to_plot %in% tres[tres$padj <= pv, 'gene']]
  tvar <- paste(tmp, 'significative')
  if(tmp == length(entrez_to_plot)) tvar <- 'all significative'
  titl <- paste("Comparison:", compid, '- top', ngenes, 'for', grpn, '-', tvar)

  tannot <- tannot[colnames(plotdata), ]
  ggs <- as.character(entrez_to_plot)
  suppressWarnings(vlnplot(cpmdata = plotdata, metadata = tannot, plotdots = nrow(tannot) < 100,
    gg = ggs, orderby = tcatg, titl = titl, noncero = FALSE, datatype = datatype, rotx = 35, v = F
  ))
  suppressWarnings(vlnplot(cpmdata = plotdata, metadata = tannot, plotdots = nrow(tannot) < 100,
    gg = ggs, orderby = tcatg, titl = titl, noncero = TRUE, datatype = datatype, rotx = 35, v = F
  ))
  void <- try(dot_lines(cpmdata = plotdata, metadata = tannot, v = F,
      gg = ggs, orderby = tcatg, datatype = datatype, mysca = 'fixed'))
  return(as.character(entrez_to_plot))
}

#### t-SNE plots #### ----------------------------------------------------------
grid.tsne <- function(
    mydata,
    redtype = 'tsne',
    dims = c('tSNE_1', 'tSNE_2'),
    group_by = 'orig.group',
    nres = NULL,
    ctag = 'no_tag',
    subst = c('class', 'in'),
    showonly = NA, # show in identities
    hide = NA, # hide in identities
    plot_return = FALSE,
    # ...
    ttl = 'all',
    sampler = FALSE,
    couls = NULL,
    dolabel = FALSE,
    ptitle = TRUE,
    ncells = TRUE,
    labsize = NULL,
    legendp = 'none',
    xlims = NULL,
    ylims = NULL,
    myseed = 27,
    v = FALSE
  ){
  if(casefold(class(mydata)) == 'seurat'){
    if(v) cat('Seurat object\n')
    metadata <- mydata@meta.data
    # coordata <- cell_embedding(object = mydata, reduction = redtype)
    if(get_version(mydata) < 3){
      coordata <- mydata@dr$tsne@cell.embeddings
    }else{
      coordata <- Embeddings(object = mydata, reduction = redtype)
    }
    # print(str(coordata))
  }else{
    metadata <- mydata
    coordata <- mydata[, dims]
  }
  if(v) cat("Dimentions:", commas(dims), "\n")
  metas <- colnames(metadata)
  if(subst[1] %in% metas){ # subsetting data
    if(sum(subst[-1]%in%metadata[,subst[1]])){ # make sure class and column is in data
      celluse <- getsubset(subst, metadata, v = v)
    }else{ celluse <- rownames(metadata) }
  }else{ celluse <- rownames(metadata) }

  # choose  Identities
  if(is.null(nres)) nres <- group_by
  if(is.numeric(nres)){
    if(sum(grepl('res', metas)) == 0){ cat('No clusters in object!\n'); break }
    nres <- paste0('res.', nres)
  }
  # print(headmat(metadata))
  # print(nres)
  # print(head(celluse))
  dimdata <- data.frame(x = coordata[celluse, dims[1]],
                      y = coordata[celluse, dims[2]])
  if(nres %in% metas){
    dimdata$Ident <- metadata[celluse, nres]
  }else{
    dimdata$Ident <- metadata[celluse, metas[grep('res\\.',metas)][1] ]
  }

  if(v) cat('Ident:', nres, '\n')
  if(v) cat('Group:', group_by, '\n')

  # show or hide in identities
  if(!is.na(showonly)) if(sum(showonly %in% dimdata$Ident)) dimdata$Ident[!dimdata$Ident %in% showonly] <- NA
  if(!is.na(hide)) if(sum(hide %in% dimdata$Ident)) dimdata$Ident[dimdata$Ident %in% hide] <- NA
  if(group_by %in% metas) dimdata$Group <- metadata[celluse, group_by] else dimdata$Group <- dimdata$Ident
  if(ctag!=''&ctag %in% metas) dimdata$ctag <- metadata[celluse, ctag] # add column for tags

  if(Hmisc::all.is.numeric(dimdata$Ident) && !dolabel) dimdata$Ident <- paste('Cluster', dimdata$Ident)
  if(Hmisc::all.is.numeric(dimdata$Group) && !dolabel) dimdata$Group <- paste('Cluster', dimdata$Group)
  dimdata$Ident <- factor(dimdata$Ident, gtools::mixedsort(unique(dimdata$Ident)))
  dimdata$Group <- factor(dimdata$Group, gtools::mixedsort(unique(dimdata$Group)))
  cls <- as.character(gtools::mixedsort(unique(dimdata$Ident)))

  if(sum(couls[1] == 'pal')){ # selecting colours
    if(v) cat('Colours from customed palette\n')
    couls <- c('#9e379f','#beba46','#00b159','#e86af0','#00aedb','#854442','#000000',
         '#493267','#49796b','#ae0001','#d11141','#ffc425','#f37735','#3b5998',
         '#be9b7b','#8caba8','#bf9b30','#18392b','#8ae429')
    couls <- couls[1:length(cls)]
    names(couls) <- cls
  }else{
    couls <- v2cols(select = cls, sour = couls, v = v)
  }

  gps <- unique(dimdata$Group)
  if(sampler){
    if(v) cat('Sampling ')
    cells <- vector()
    mincells <- min(table(dimdata$Group))
    for(gr in gps){
      sset <- dimdata$Group == gr
      if(v) cat('|')
      set.seed(myseed)
      cells <- c(cells, sample(rownames(dimdata[sset, ]), mincells))
    }; if(v) cat('\n')
    dimdata <- dimdata[cells, ]
  }
  dimdata <- dimdata[gtools::mixedorder(dimdata$Ident), ]

  # if(length(gps) %% 2) gps<-c('all', gps)
  if(plot_return == 'struct'){
    if(v) cat('Returning structure\n')
    return(dimdata)
  }
  if(v) cat('Building plots ')
  if(grepl('facet', plot_return)){
    if(v) cat('facet\n')
    ploties <- tsne_facet(
        dimdata = dimdata,
        sampler = FALSE,
        couls = couls,
        dolabel = dolabel,
        legendp = legendp,
        myseed = myseed,
        v = v
      )
  }else{
    if(v) cat('tsnep\n')
    registerDoParallel(cores = parallel::detectCores())
    # ploties <- foreach(gp = gps) %dopar% {
    ploties <- lapply(gps, function(gp){
      tsnep(dimdata = dimdata, ttl = gp, sampler = FALSE, couls = couls,
        dolabel = dolabel, ptitle = ptitle, ncells = ncells, labsize = labsize,
        legendp = legendp, xlims = NULL, ylims = NULL, myseed = 27, v = v)
      })
  }

  if(plot_return == 'facet_return' || isTRUE(plot_return)){
    if(v) cat('Returning plots\n')
    return(ploties)
  }
  if(v) cat('Plotting\n')
  if(grepl('facet', plot_return)){
    print(ploties)
  }else{
    ssize <- fitgrid(ploties)
    do.call("grid.arrange", c(ploties, nrow = ssize[1], ncol = ssize[2]))
  }
}
#
tsne_facet <- function(
    dimdata,
    sampler = FALSE,
    couls = NULL,
    dolabel = FALSE,
    ptsize = NULL,
    legendp = 'none',
    ncells = TRUE,
    myseed = 27,
    v = FALSE
  ){
  gps <- unique(dimdata[, 'Group'])

  if(ncells){
    if(v) cat(' .nCells ')
    for(gr in gps){
      sset <- dimdata$Group == gr
      dimdata[sset, 'nCells'] <- paste0(dimdata[sset, 'Group'],'\nT: ', sum(sset), ' cells')
    }; if(v) cat('\n')
  }
  if(sampler){
    if(v) cat(' .sampling\n')
    cells <- vector()
    mincells <- min(table(dimdata$Group))
    for(gr in gps){
      sset <- dimdata$Group == gr
      if(v) cat('|')
      set.seed(myseed)
      cells <- c(cells, sample(rownames(dimdata[sset, ]), mincells))
      if(ncells) dimdata[sset, 'Group'] <- paste0(dimdata[sset, 'Group'],'\n', mincells, ' cells')
    }; if(v) cat('\n')
    dimdata <- dimdata[cells, ]
  }
  if(!is.factor(dimdata[, 'Ident'])){
    if(v) cat(' .factors\n')
    dimdata[, 'Ident'] <- factor(dimdata[, 'Ident'], gtools::mixedsort(unique(dimdata[, 'Ident'])))
  }

  ploty <- ggplot(dimdata, aes(x, y) ) + geom_point(shape = 20, aes(colour = Ident))
  if(!is.null(ptsize)) ploty <- ploty + geom_point(shape = 20, aes_string(colour = 'Ident'), size = ptsize)
  if(v) cat(' .')
  couls <- v2cols(select = levels(dimdata[, 'Ident']), sour = couls, v = v)
  ploty <- ploty + scale_colour_manual(values = couls, na.value = '#CCCCCC') +
    guides(colour = guide_legend(override.aes = list(size = 6)))
  # if(dolabel){
  #   if(v) cat(' .labs\n')
  #   dimdata %>%
  #     group_by(Ident) %>%
  #     summarize(x = median(x = x), y = median(x = y)) -> centers
  #   ploty <- ploty + geom_text(data = centers, mapping = aes(label = Ident))
  # }
  ploty <- ploty +
    facet_wrap(as.formula(paste('~', 'Group')), ncol = fitgrid(gps)[1]) +
    labs( x= "tSNE_1", y= "tSNE_2") +
    theme_classic() + theme(panel.background = element_blank(),
          axis.line = element_line(colour = "grey80"),
          legend.position = legendp,
          strip.background = element_rect(fill = "#FFFFFF", linetype = 0))
  return(ploty)
}
# individual tSNE
tsnep <- function(
    dimdata,
    ttl = 'all',
    Ident = 'Ident',
    Group = 'Group',
    ctag = 'no_tag',
    dimname = 'DIM',
    sampler = FALSE,
    couls = NULL,
    dolabel = FALSE,
    ptitle = FALSE,
    ncells = FALSE,
    ptsize = NULL,
    labsize = NULL,
    legendp = 'none',
    xlims = NULL,
    ylims = NULL,
    myseed = 27,
    return_plot = TRUE,
    v = FALSE
  ){
  if(class(dimdata) != 'data.frame') dimdata <- as.data.frame(dimdata, stringsAsFactors = F)
  if(v) cat(' .tags - ')
  if(Group %in% colnames(dimdata)){
    tvar <- unique(dimdata[, Group])
    if(v) cat(Group, 'in', commas(tvar), '\n')
    if(ttl %in% tvar){
      if(v) cat('  +', commas(ttl), '\n')
      sset <- dimdata[, Group] %in% ttl
      if(sampler){
        ncells <- TRUE
        mincells <- min(table(dimdata[, Group]))
        set.seed(myseed)
        cells <- sample(rownames(dimdata[sset,]), mincells)
      }else{
        cells <- rownames(dimdata[sset,])
      }
    }else{
      if(v) cat('  - ')
      if(v && !sum(ttl %in% 'all')) cat(commas(ttl), '\n') else if(v) cat('all cells\n')
      sset <- dimdata[, Group] != ttl
      cells <- rownames(dimdata[sset,])
    }
  }else{
    if(v) cat('no grouping\n')
    sset <- dimdata != ttl
    cells <- rownames(dimdata)
  }

  # Get tag
  if(ctag %in% colnames(dimdata)){
    if(v) cat('\n .tags\n')
    tvar <- paste(unique(dimdata[sset, ctag]),collapse='-')
    ttl <- paste(ttl, tvar, sep=': ')
  }
  if(ttl[1] == 'all') ttl <- 'All cells'

  if(is.null(xlims)) xlims <- c(min(dimdata$x), max(dimdata$x))
  if(is.null(ylims)) ylims <- c(min(dimdata$y), max(dimdata$y))
  ploty <- ggplot(dimdata[cells, ], aes(x, y) ) +
    geom_point(shape = 20, aes_string(colour = Ident)) +
    coord_cartesian(xlim = xlims, ylim = ylims) +
    theme_classic() +
    theme(legend.position = legendp,
      plot.title = element_text(hjust = 0.5, face = 'bold')
    ) + labs(x = paste0(dimname, "_1"), y = paste0(dimname, "_2"))
  if(!is.null(ptsize)) ploty <- ploty + geom_point(shape = 20, aes_string(colour = Ident), size = ptsize)
  if(ncells) tvar <- paste0(length(cells),' cells') else tvar <- ''
  if(sampler){
    if(v) cat('  .sample\n')
    tvar <- paste0(tvar, '(',sum(sset),')')
  }
  if(ptitle){
    if(v) cat(' .title\n')
    ploty <- ploty + ggtitle(label = ttl, subtitle = tvar)
  }#else{ ploty <- ploty + theme(title = element_blank()) }

  if(!is.factor(dimdata[, 'Ident'])){
    if(v) cat(' .factors\n')
    dimdata[, Ident] <- factor(dimdata[, Ident], gtools::mixedsort(unique(dimdata[, Ident])))
  }
  if(v) cat(' .'); couls <- v2cols(select = levels(dimdata[, Ident]), sour = couls, v = v)
  # if(sum(couls[1] != 'r')) ploty <- ploty + scale_colour_manual(values = couls, breaks = unique(dimdata$Ident), labels = names(couls), na.value = '#CCCCCC')
  ploty <- ploty + scale_colour_manual(values = couls, na.value = '#CCCCCC') +
    guides(colour = guide_legend(override.aes = list(size = 6)))
  if(sum(is.na(dimdata[, Ident]))) ploty <- ploty + scale_color_discrete(na.value = "#CCCCCC")
  # if (dolabel) {
  #   if(v) cat(' .labels\n')
  #   dimdata %>%
  #     group_by(Ident) %>%
  #     summarize(x = median(x = x), y = median(x = y)) -> centers
  #   ploty <- ploty +
  #     geom_point(data = centers, mapping = aes(x = x, y = y), size = 0, alpha = 0)
  #   if(is.null(labsize)){
  #     ploty <- ploty + geom_text(data = centers, mapping = aes(label = Ident))
  #   }else{
  #     ploty <- ploty + geom_text(data = centers, mapping = aes(label = Ident), size = labsize)
  #   }
  # }
  if(return_plot) return(ploty)
  if(v) cat(' .plotting\n')
  print(ploty)
}

# volcano plot
volplot <- function(x,
  pvalth = 0.05,
  lfcth = 0.5,
  pvaltype = 'padj',
  lfctype = 'log2FoldChange',
  gene_name = NULL,
  group = NULL,
  interact = FALSE,
  ngenes = 5,
  legends = FALSE,
  titl = "Volcano Plot",
  subt = NULL,
  check_genes = FALSE,
  gname = "Markers",
  only_degs = FALSE,
  return_plot = FALSE,
  clipp = FALSE,
  mycolours = c("#b2c5ca", "#008080", "#b7dac5", "#b22028", "#18f6e7"),
  v = FALSE
){
  x <- data.frame(x, stringsAsFactors = F)
  # add a grouping column; default value is "not significant"
  if(is.null(gene_name) && !'gene_name' %in% colnames(x)){
    x["gene_name"] <- rownames(x)
  }else{
    colnames(x) <- sub(gene_name, "gene_name", colnames(x))
  }
  rownames(x) <- x[, "gene_name"]
  colnames(x) <- sub(pvaltype, "FDR", colnames(x))
  colnames(x) <- sub(lfctype, "Fold", colnames(x))
  if(v) cat("Setting groups\n")
  if(is.null(group)){
    x["group"] <- "Not_significant"
    # change the grouping for the entries with significance but not a large enough Fold change
    x[which(x['FDR'] < pvalth & abs(x['Fold']) > lfcth ), "group"] <- "Significant"

    # change the grouping for the entries a large enough Fold change but not a low enough p value
    x[which(x['FDR'] > pvalth & abs(x['Fold']) > lfcth ), "group"] <- "FC"

    # change the grouping for the entries with both significance and large enough fold change
    x[which(x['FDR'] < pvalth & abs(x['Fold']) > lfcth ), "group"] <- "DEG"
  }else{ x["group"] <- x[group] }
  x[, 'FDR'] <- (-log10(x[, 'FDR']))
  if(v) cat(commas(unique(x[, "group"])), "\n")
  # x <- x[order(x['FDR']), ]
  # if(max(abs(diff(head(x[, 'FDR'], 40)))) > 100) # checking if there's a huge difference between the top

  if(v) cat("Find and label the top peaks\n")
  if(check_genes[1] != "" && is.character(check_genes[1])){
    if(v) cat("Given genes\n")
    tvar <- check_genes %in% rownames(x)
    if(sum(!tvar)) warning("Some genes are not in data: ", commas(check_genes[!tvar]))
    tvar <- check_genes[tvar]
    if(length(tvar)) x[tvar, "group"] <- gname else tvar <- ''
  }else{ tvar <- '' }
  if(ngenes > 0){
    if(v) cat("Top", ngenes, "genes\n")
    check_genes <- head(rownames(x[with(x, order(Fold, FDR)),]), ngenes)
    check_genes <- unique(c(check_genes, head(rownames(x[with(x, order(-Fold, FDR)),]), ngenes)))
    check_genes <- unique(c(rbind(check_genes, head(rownames(x[with(x, order(-FDR)), ]), ngenes))))
  }else{ check_genes <- rownames(x[with(x, order(-FDR)), ])[1] }
  check_genes <- check_genes[!is.na(check_genes)]
  if(tvar[1] != "") check_genes <- unique(c(head(tvar, ngenes), check_genes))
  check_genes <- x[check_genes, ]
  # check_genes <- x[x$group == 'DEG', ]

  mycats <- unlist(unique(x["group"]))
  if(v) cat("Setting colours and labels\n")
  if(is.null(names(mycolours))){
    names(mycolours) <- c("Not_significant", "Significant", "FC", "DEG", gname)
  }else if(length(mycolours) != length(unique(x["group"])))
    mycolours <- v2cols(mycats, mycolours)
  signf <- paste("Significance:", pvalth,"/ FC:", lfcth)
  signf <- paste(signf, "\n No. of DEGs:", sum(x$group %in% 'DEG'))
  if(!is.null(subt)){
    subt <- paste(subt, " | ", signf)
  }else{ subt <- signf }

  # make the Plot.ly plot
  if(!interact){
    if(v) cat("Static plot\n")
    p <- ggplot(x, aes(x = Fold, y = FDR)) +
      xlim(c(-max(abs(x$Fold)), max(abs(x$Fold)))) +
      geom_point(data = x, aes(color = group)) +
      geom_hline(yintercept = -log10(pvalth), linetype = "dashed", color = "gray") +
      geom_vline(xintercept = -lfcth, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = lfcth, linetype = "dashed", color = "gray")
    # clipit <- quantile(-log10(res$padj), prob = .9)
    # sum(-log10(res$padj) > clipit) / length(res$padj) < 0.1
    # clipit <- quantile(x[, 'FDR'], prob = .9) # if less than 10 % are over this value
    # if(sum(x[, 'FDR'] > clipit) / length(x[, 'FDR']) < 0.1 && clipp){ # clip it to the 3rd Quantile
    if(isTRUE(clipp)) clipp <- 0.999
    if(is.numeric(clipp)){
      clipp <- if(clipp < 1) quantile(x[, 'FDR'], prob = clipp) else x[tail(order(x[, 'FDR']), clipp)[1], 'FDR']
      p <- p + coord_cartesian(ylim = c(0, clipp))
    }
    if(only_degs){
      if(v) cat("DEGs inly colouring\n")
      # subdata <- dplyr::filter(x, abs(Fold) <= lfcth | FDR >= pvalth)
      p <- p + stat_density2d(data = x[!x$group %in% c('DEG', gname), ],
          aes(fill = ..density..^0.1), geom = "raster",  contour = FALSE, show.legend = FALSE) +
        scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256), breaks = NULL) +
        geom_point(data = x[x$group %in% c('DEG', gname), ], aes(color = group))
    }else{
      p <- p + geom_point(data = x[rownames(check_genes), ],
        aes(x = Fold, y = FDR), size = 2, color = names(mycolours)[gname]) + theme_bw()
    }
    if(length(rownames(check_genes))){
      if(v) cat("Repeling names\n")
      p <- p + geom_text_repel(
          data = x[rownames(check_genes), ],
          aes(label = gene_name),#  size = 7,
          box.padding = unit(0.35, "lines"),
          point.padding = unit(0.3, "lines")
        )
    }
    p <- p + scale_color_manual(values = mycolours) +
      labs(title = titl, subtitle = subt, x = 'Fold Change', y = 'FDR') +
      theme(legend.position = ifelse(legends, 'right', 'none'),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.background = element_blank(),
        panel.border = element_rect(fill = F))
      if(!legends) p <- p + theme(plot.margin = unit(c(0.1, 2, 0.1, 0.1), "cm"))
  }else{
    if(v) cat("Interactive plot\n"); a <- list()
    if(v) cat("Bulding layout\n")
    for (i in seq_len(nrow(check_genes))) {
      m <- check_genes[i, ]
      a[[i]] <- list(
        x = m[["Fold"]],
        y = m[["FDR"]],
        text = m[["gene_name"]],
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 0.5,
        ax = 20,
        ay = -40
      )
    }
    titl <- paste0(titl, '\n', subt)
    p <- plot_ly(data = x, x = ~Fold, y = ~FDR, text = ~gene_name, mode = "markers", color = ~group) %>%
      layout(title = titl) %>%
      layout(annotations = a)
  }
  if(return_plot || interact) return(p) else print(p)
}

# For pca plots
pcaplot <- function(x, main.f, coldata, batchCol, cd,
                    ptsize = 3, printplot = TRUE, logp = T,
                    cols = c('condition','batch'), shapes = 'batch',
                    grcols = NULL, v = FALSE, show_var = F, hidleg = FALSE){
  if(isTRUE(logp)){ x <- log2(x + 1); if(v) cat('log2(cts+1)\n') }
  # coldata$condition <- coldata[, cd]
  # coldata$batch <- coldata[, batchCol]
  coldata[, cd] <- coldata[, cd]
  if(batchCol %in% colnames(coldata)) coldata[, batchCol] <- coldata[, batchCol] else coldata[, batchCol] <- "1"
  if(shapes == 'batch') coldata[, 'shapes'] <- coldata[, batchCol]
  cols <- c(cd, batchCol)
  if(v) cat('PCA calculation\n')
  # pc <- prcomp(x)
  pc <- irlba::prcomp_irlba(x)
  pcdf <- as.data.frame(pc$x)
  pcaVars=signif(((pc$sdev)^2)/(sum((pc$sdev)^2)),3)*100
  pcdf <- data.frame(cbind(pcdf, coldata), stringsAsFactors = F)
  if(is.null(shapes)){
    svalues <- NULL
  }else if(shapes[1] %in% colnames(pcdf)){
    svalues <- 1:length(unique(pcdf[,shapes[1]]))
  }else{ svalues <- 1:length(unique(pcdf[,'batch'])) }
  if(length(svalues) <= 4) svalues <- 16:(16+length(svalues))
  p <- list()
  if(v) cat('Building plots\n')
  for(k in 1:length(cols)){
    for(i in 1:3){ for(j in 1:3){ if(i<j){
      p1 <- ggplot(pcdf, aes_string(x=paste0('PC',i),y=paste0('PC',j), color=cols[k])) +
        theme_bw() + theme(plot.title=element_text(hjust=0.5),
          panel.background = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
      if(!is.null(shapes)){
        p1 <- p1 + geom_point(aes_string(shape = shapes), size = ptsize) +
          scale_shape_manual(values=svalues)
      }else{
        p1 <- p1 + geom_point(size = ptsize)
      }
      if(show_var){
        p[[paste(k,i,j)]] <- p1 + xlab(paste0('PC',i,': ',pcaVars[i],' % of Variance Explained')) +
        ylab(paste0('PC',j,': ',pcaVars[j],' % of Variance Explained'))
      }else{
        p[[paste(k,i,j)]] <- p1
      }
      if(i == 2) p1 <- p1 + ggtitle(main.f) else p1 <- p1 + theme(legend.position="none")
      if(isTRUE(hidleg)) p1 <- p1 + theme(legend.position="none")
      p[[paste(k,i,j)]] <- p1
    }}}
  }
  if(isTRUE(printplot)){
    if(v) cat('Plotting\n')
    invisible(lapply(p, print))
  }else if(printplot == 'grid'){
    if(v) cat('Plotting grid\n')
    pca <- data.table::data.table(pcdf[,1:3], as.data.frame(coldata[,cd]))
    tvar <- unique(coldata[, cd])
    if(is.null(grcols)){
      grcols <- gg_color_hue(length(tvar))
      names(grcols) <- tvar
    }
    cols <- grcols[tvar]
    colnames(pca)=c("PC1","PC2","PC3",cd)
    p <- ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', cd), cardinality_threshold = 27,
      mapping = aes_string(color = cd), upper = list(continuous = 'density'),
      title = main.f)
    if(sum(is.na(cols))==0){
      names(cols) <- tvar
      try(for(i in 1:p$nrow) {
        for(j in 1:p$ncol){
          p[i,j] <- p[i,j] +
              scale_fill_manual(values = cols) +
              scale_color_manual(values = cols)
        }
      })
    }
    suppressMessages( print(p,progress = F))
  }else{ return(p); if(v) cat('Returning plot\n') }
}

# copied from deseq2
plotPCa <- function(
  object,
  intgroup = "condition",
  ntop = 500,
  returnData = FALSE,
  legend = FALSE,
  titl = ""
){
  rv <- rowVars(assay(object))

  # select the ntop genes by variance
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]

  # perform a PCA on the data in assay(x) for the selected genes
  pca <- prcomp(t(assay(object)[select,]))

  # the contribution to the total variance for each component
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )

  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }

  if(length(intgroup) > 1){
    batch <- colData(object)[[intgroup[1]]]
    intgroup <- intgroup[-1]
  }else{
    batch <- factor('no')
  }
  intgroup.df <- as.data.frame(colData(object)[, intgroup, drop=FALSE])

  # add the intgroup factors together to create a new grouping factor
  group <- if (length(intgroup) > 1) {
    factor(apply( intgroup.df, 1, paste, collapse=":"))
  } else {
    colData(object)[[intgroup]]
  }

  # assembly the data for the plot
  d <- data.frame(PC1 = pca$x[, 1], PC2=pca$x[, 2],
                  group = group, batch = batch,
                  intgroup.df, name=1:ncol(object))
  if (returnData) {
    attr(d, "percentVar") <- percentVar[1:2]
    return(d)
  }
  tvar <- nlevels(d$batch)
  myvals <- if(length(tvar) < 11){
    14:(tvar + 14)
  }else{ 1:tvar }

  p <- ggplot(data=d, aes_string(x="PC1", y="PC2", color="group")) +
    xlab(paste0("PC1: ",round(percentVar[1] * 100),"% variance")) +
    ylab(paste0("PC2: ",round(percentVar[2] * 100),"% variance")) +
    coord_fixed() + ggtitle(titl)
  if(batch[1] != 'no'){
    p <- p + geom_point(size=3, aes_string(color="group", shape="batch")) +
    scale_shape_manual(values=myvals)
  }else{
    p <- p + geom_point(size=3, aes_string(color="group"))
  }
  if(!legend) p <- p + theme(legend.position = "none")
  return(p)
}

# another pca plot
plotPCA <- function(
  sca_obj,
  coldata = NULL,
  condition = 'Condition',
  colours = NULL,
  compute = TRUE,
  spike = FALSE,
  return_plot = FALSE,
  v = FALSE
){
  require(GGally)
  if(class(sca_obj) == 'SingleCellAssay'){
    mat <- assay(sca_obj)
    coldata <- as.data.frame(colData(sca_obj))
  }else if(casefold(class(sca_obj)) == 'seurat'){
    vrs <- get_version(sca_obj)
    if(!compute){
      compute <- FALSE;
      if(vrs < 3) projection <- sca_obj@dr$pca@cell.embeddings else projection <- Embeddings(sca_obj, "pca")
    }
    if(vrs < 3) mat <- sca_obj@data else mat <- GetAssayData(sca_obj)
    coldata <- sca_obj@meta.data
  }else if(sum(dim(sca_obj) > 1)) mat <- as.matrix(sca_obj)
  if(is.null(coldata)){
    coldata <- data.frame(condition = rep('Condition', nrow(projection)))
    colnames(coldata) <- condition <- 'Condition'
  }
  if(isTRUE(spike)){
    if(v) cat('Spiking root(s)\n')
    grps <- coldata[, condition[1]]; names(grps) <- rownames(coldata)
    roots <- get_stat_report(mat, groups = grps, moments = c('bm', 'mn'), v = v)[, "Bmean", drop = F]#-c(1:4)]
    mat <- cbind(mat, roots)
    coldata <- rbind(coldata, mat_names(colnames(roots), colnames(coldata)))
    coldata[colnames(roots), ] <- rep(colnames(roots), ncol(coldata))
  }
  condition <- getfound(condition, colnames(coldata), element = 'column(s)', v = v)
  for(i in condition) coldata[, i] <- factor(coldata[, i])
  grps <- levels(coldata[, condition[1]])
  mat <- mat[, rownames(coldata)]
  if(v) cat('Groups:', commas(grps), '\n')
  if(v) cat('Genes:', nrow(mat), '\nCells:', ncol(mat), '\n')
  if(compute){
    projection <- irlba::prcomp_irlba(t(mat), retx = TRUE, center = TRUE, scale. = FALSE)$x
  }
  rownames(projection) <- colnames(mat)
  cols <- v2cols(grps, colours); firstpcs <- c('PC1', 'PC2', 'PC3')
  pca <- data.table::data.table(projection, coldata[rownames(projection), condition])
  colnames(pca) <- c(paste0("PC", 1:ncol(projection)), condition)
  myplots <- list()
  for(cnd in condition){
    p <- ggpairs(pca, columns = c(firstpcs[firstpcs %in% colnames(pca)], cnd),
            cardinality_threshold = 27, mapping = aes_string(color = condition[1]),
            upper = list(continuous = 'density'), diag = list(continuous = 'density')
          ) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
    for(i in 1:p$nrow) { for(j in 1:p$ncol){
        p[i, j] <- p[i,j] + scale_fill_manual(values = cols) + scale_color_manual(values = cols)
    }}
    if(return_plot){ myplots[[cnd]] <- p; next }
    suppressMessages(print(p, progress = F))
  }
  if(length(myplots) == 1) myplots <- myplots[[1]]
  if(return_plot) return(myplots)
  invisible(pca)
}

# proportions
get_props <- function(
    metadata = NULL,
    props = NULL,
    group_by = NULL, # normalises with C
    resolution = NULL, # normalises with R
    group = NULL,
    cluster = NULL,
    facet = FALSE,
    repel_text = TRUE,
    writeto = NULL,
    norm_type = c('n', 'c', 'r', 'b'),
    couls = NULL,
    return_plot = TRUE,
    pies = TRUE,
    ssize = NULL,
    # write_raw = TRUE,
    reverseit = FALSE,
    v = FALSE
){
  norm_type <- match.arg(norm_type)
  if(!is.null(colnames(metadata))){
    if(is.null(group_by)) group_by <- colnames(metadata)[1]
    if(is.null(resolution)){
      tvar <- colnames(metadata)
      resolution <- tvar[tvar != group_by][1]
    }
  }else{
    group_by <- 1
    resolution <- 2
  }
  if(v) cat('Getting',resolution, '/', cluster, 'proportions in', group_by, '/', group, '\n')
  if(is.null(props) && !is.null(metadata)){
    props <- t(table(metadata[, group_by], metadata[, resolution]))
  }
  if(v) print(dim(props))
  if(v) print(head(props))
  if(norm_type == "b") norm_type <- c('c', 'r')
  if(!sum(norm_type %in% c("b", "c", "r"))){
    warning('No normalisation performed')
    props_norm <- props
  }
  if("c" %in% norm_type){
    if(v) cat(' - Normalising by columns\n')
    props_norm <- round(sweep(props, 2, colSums(props),'/') * 100.0, 4)
  }else{ props_norm <- props }
  if("r" %in% norm_type){
    if(v) cat(' - Normalising by rows\n')
    props_norm <- round(sweep(props_norm, 1, rowSums(props_norm),'/')*100.0, 4)
  }

  # pie charts
  mydata <- data.table::melt(props_norm)
  mydata <- mydata[order(mydata[, 1]), ]
  if(is.null(group)) group <- sub('orig\\.', '', group_by)
  if(is.null(cluster)) cluster <- sub('orig\\.', '', resolution)
  colnames(mydata) <- c(cluster, group, 'N')
  mydata$N <- as.numeric(mydata$N)
  if(!is.factor(mydata[, cluster])) mydata[, cluster] <- factormix(mydata[, cluster])
  if(!is.factor(mydata[, group])) mydata[, group] <- factormix(mydata[, group])
  mycouls <- v2cols(select = levels(mydata[, group]), sour = couls, v = v)
  titl <- paste("In", cluster, "-", group, "proportions")

  if(v) cat(' - Bar plots\n')
  p <- ggplot(data = mydata, aes_string(x = cluster, y = 'N', fill = group)) +
    geom_bar(stat = "identity", position = "fill", alpha = 0.7)
  if(isTRUE(facet)) p <- p + facet_wrap(as.formula(paste('~', cluster)), ncol = 1)
  p <- p + labs(fill = group, title = titl, x = resolution, y = 'Percentage') +
    theme_classic() + theme(panel.background = element_blank(), panel.grid.major = element_blank(),
      strip.text.x = element_blank(), strip.background = element_blank()) +
    scale_fill_manual(name = "Group", values = mycouls[levels(mydata[, group])])
  if(is.null(ssize)) ssize <- fitgrid(table(mydata[, cluster]))
  if(pies){
    if(v) cat(' - Pie chart\n')
    mydata$cluster <- mydata[, cluster]
    mydata$group <- mydata[, group]
    thisdf <- mydata %>% group_by(cluster) %>%
      mutate(group = factor(group, rev(levels(group))),
        valuep = round(N / sum(N) * 100, 1),
        cumulative = cumsum(valuep),
        position = cumulative - valuep / 2,
        pct = paste0(valuep, "%")
      )
    levels(thisdf$group)
    thisdf[, 1] <- paste0(mydata[, 1], '\nN = ', rowSums(props_norm)[mydata[, 1]])
    thisplot <- ggplot(thisdf, aes(x = "", y = valuep)) +
      geom_bar(aes_string(fill = group), width = 1, alpha = 0.7, size = 1,
        color = "white", stat = "identity",
        position = position_stack(reverse = reverseit)) + # NOTE: Problems with the stacking
      facet_wrap(as.formula(paste("~", cluster)), ncol = ssize[1]) +
      labs(x = NULL, y = NULL, fill = group, title = titl) +
      scale_fill_manual(name = group, values = mycouls[levels(mydata[, group])]) +
      coord_polar(theta = "y", start = 0)
    if(repel_text){
      thisplot <- thisplot + geom_text_repel(data = thisdf[thisdf[, 'valuep'] > 0, ],
        aes(x = 1.3, y = position, label = pct), show.legend = F, nudge_x = .5, segment.alpha = .5)
    }
    thisplot <- thisplot + theme_classic() +
      theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        strip.text.x = element_text(face = 'bold'),
        strip.background = element_rect(fill = "#FFFFFF", linetype = 0))
  }

  if(!is.null(writeto)){
    dir.create(writeto)
    if(v) cat(' - Writing to:', writeto,'\n')
  }
  if(!return_plot){
    fname <- paste0(writeto, group,'_proportions_in_', cluster,'.pdf')
    if(v) cat(' - Plotting: ', fname,'\n')
    pdf(fname, 4*ssize[3], 4*ssize[3])
    print(plot_grid(p, thisplot, ncol = 1))
    dev.off()
  }

  myprops <- calc_tots(props_norm)
  # if(write_raw){
  # if(!writenorm || !sum(norm_type %in% c("b", "c", "r")) && !is.null(metadata)){
  if(sum(norm_type %in% c("b", "c", "r"))){
    if(v) cat(' - Calculating totals\n')
    myprops <- cbind(myprops, calc_tots(props))
  }

  if(!is.null(writeto)){
    tvar <- paste0(group, '_proportions_in_', cluster,'.csv')
    if(v) cat(' - Writing table', tvar,'\n')
    write.csv(myprops, file = paste0(writeto, tvar))
  }
  if(return_plot) return(list(bars = p, pies = thisplot, table = myprops))
  return(myprops)
}

# SEM for genes
pSEM <- function(
  mygene,
  edata,
  metadata,
  group_by = NULL,
  colour_by = NULL,
  sourcols = NULL,
  group_name = NULL,
  mylevels = NULL,
  ctstype = "CTS",
  titl = NULL,
  l2 = TRUE,
  legendpos = "none",
  ptsize = 1,
  return_plot = FALSE,
  v = FALSE
){
  if(is.null(group_by)) group_by <- colnames(cdata)[ncol(cdata)]
  if(is.null(colour_by)) colour_by <- group_by
  if(is.null(group_name)) group_name <- group_by
  cdata <- metadata[metadata[, group_by] != 'no.class', ]
  edata <- as.matrix(edata)

  tvar <- mygene %in% rownames(edata)
  if(sum(tvar) == 0){ cat('No genes to plot\n'); return(0) }
  mygene <- getfound(mygene, rownames(edata), element = 'Genes', v = FALSE)

  if(l2){
    if(v) cat('Log2 transform\n')
    ctstype <- paste0("log2(", ctstype, "+1)")
  }

  ddf <- data.table::rbindlist(lapply(mygene, function(gg){
    if(!gg %in% rownames(edata)) print(head(edata))
    tvar <- edata[gg, rownames(cdata), drop = F]
    if(l2) tvar <- log2(tvar + 1)
    mydf <- data.frame(
      Expression = as.numeric(tvar),
      Gene = gg,
      cGroups = cdata[[colour_by]],
      xGroups = cdata[[group_by]]
    )
  }))
  if(!is.null(mylevels)){
    if(sum(mylevels %in% ddf$cGroups)){
      ddf$cGroups <- factor(ddf$cGroups, levels = group_order(unique(ddf$cGroups), mylevels))
    }
    if(sum(mylevels %in% ddf$xGroups)){
      ddf$xGroups <- factor(ddf$xGroups, levels = group_order(unique(ddf$xGroups), mylevels))
    }
  }
  gcols <- v2cols(levels(ddf$cGroups), sourcols)
  if(group_by != colour_by) gcols <- c(gcols, v2cols(levels(ddf$xGroups), sourcols))
  if(!is.null(titl)) titl <- paste(titl, '- ')
  titl <- paste0(titl, ifelse(length(mygene) > 1, 'Standard Error of the Mean', mygene))

  g1 <- ggplot(ddf, aes(x = xGroups, y = Expression, color = cGroups)) +
    geom_jitter(position = position_jitter(0.2), size = ptsize) +
    stat_summary(fun.y = mean, geom = "point", color = "black", size = 2) +
    stat_summary(fun.data = mean_se, geom = "errorbar", color = "black") +
    scale_color_manual(values = v2cols(levels(ddf$cGroups), gcols)) +
    theme_bw() + theme(legend.position = legendpos,
      plot.title = element_text(hjust=0.5, face = "bold"),
      panel.background = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
      strip.text = element_text(size=11, face = "bold.italic"),
      axis.text.x = element_text(size=10, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size=13, face = "bold"),
      axis.title.y = element_text(size=12)) +
    labs(x = "Groups", y = ctstype, title = titl, subtitle = group_name)
  if(length(mygene) > 1) g1 <- g1 + facet_wrap(~ Gene, scale = 'fixed')#'free_y')
  if(return_plot) return(g1) else print(g1)
}

# correlation between two variables
plot_corr <- function(
  df,
  var1 = 1,
  var2 = 2,
  addvar = NULL,
  lm_formula = NULL,
  use_residuals = NULL,
  mthod = "spearman",
  log2t = FALSE,
  add_line = TRUE,
  out_path = './',
  suffix = NULL,
  return_plot = FALSE,
  ptsize = 1.5,
  maxdiff = 500,
  v = FALSE
){
  if(v) cat('Correlation:', var1, 'vs', var2, '\n')
  if(!isTRUE(addvar %in% colnames(df))) addvar <- NULL
  myvars_names <- c(var1, var2, addvar)
  df_plot <- df[, myvars_names]
  myvars <- c("X1", "X2")
  if(!is.null(addvar)) myvars <- c(myvars, "Group", "Group2")[1:(length(addvar) + 2)]
  myaes <- aes(x = X1, y = X2)
  if(length(myvars) == 3) myaes <- aes(x = X1, y = X2, color = Group)
  if(length(myvars) == 4) myaes <- aes(x = X1, y = X2, color = Group, shape = Group2)
  colnames(df_plot) <- myvars
  if(isTRUE(log2t)){
    df_plot[, c("X1", "X2")] <- log2(df_plot[, c("X1", "X2")] + 1)
    suffix <- paste0(suffix, '_log2')
  }else if(is.character(log2t) && FALSE){
    if(v) cat("Auto logarithm transformation\n")
    for(myvar in myvars){
      if(diff(range(df_plot[, myvar])) > maxdiff){
        df_plot[, myvar] <- log2(df_plot[, myvar])
        myvars_names[which(myvars_names == myvar)] <- paste0("log10_", myvars_names[which(myvars_names == myvar)])
        if(isTRUE(addvar %in% myvar)) addvar <- paste0("log10_", addvar)
      }
    }
  }
  df_plot$X1 <- as.numeric(df_plot$X1)
  df_plot$X2 <- as.numeric(df_plot$X2)
  df_plot <- joindf(df_plot, df)
  if(is.null(lm_formula)) lm_formula <- "X1 ~ X2"
  lm_formula <- gsub(" ", "", lm_formula)
  if(v) cat("Formula:", lm_formula)
  if(!is.null(use_residuals)){
    tvar <- try(lm(formula = as.formula(paste0("X2~", use_residuals)), data = df_plot), silent = TRUE)
    if(class(tvar) != 'try-error'){
      m1 <- try(lm(formula = as.formula(paste0("X1~X2+", use_residuals)), data = df_plot), silent = TRUE)
      m2 <- try(lm(formula = X1~X2, data = df_plot), silent = TRUE)
      m1rss <- if (is.null(m1$weights)) sum(residuals(m1)^2) else sum(m1$weights * residuals(m1)^2)
      m2rss <- if (is.null(m2$weights)) sum(residuals(m2)^2) else sum(m2$weights * residuals(m2)^2)
      use_residuals <- (m2rss - m1rss)/m2rss
      df_plot$X2 <- summary(tvar)$residuals
    }else{ use_residuals <- NULL }
  }
  fit <- try(lm(formula = as.formula(lm_formula), data = df_plot), silent = TRUE)
  if(class(fit) != 'try-error'){
    coefs <- summary(fit)$coefficients
    corry <- cor(df_plot$X1, df_plot$X2, method = mthod)
    corrtest <- try(cor.test(df_plot$X1, df_plot$X2, method = mthod))
    mypvalue <- if(class(corrtest) != "try-error") corrtest$p.value else coefs[nrow(coefs), ncol(coefs)]
    if(v) print(summary(fit))
    annot_plot <- paste0(
      "r^2 = ", round( summary(fit)$r.squared, 2),
      "\np-value = ", formatC(mypvalue, 3),
      # "\np-value = ", coefs[which(rownames(coefs) %in% c("X2", var2))[1], ncol(coefs)],
      # "\n2-tailed p-value = ", 2*pt(q = -abs(coefs[, 3]), df = fit$df)[2],
      "\n", stringr::str_to_sentence(mthod), " (rho) = ", round( corry, 2)
    )
    if(!is.null(use_residuals)) annot_plot <- paste0(annot_plot, "\nProportional Reductional of Error = ", use_residuals)
    tmp <- ifelse(is.null(use_residuals), var2, paste0("residuals(", var2, ")"))
    tvar <- ifelse(gsub(" ", "", lm_formula) != "X1~X2" && is.null(use_residuals), lm_formula, paste0(var1, "~", tmp))
    annot_plot <- paste0(annot_plot, "\nFormula = ", tvar)
    if(v) cat(annot_plot, '\n')
    df_plot <- joindf(fit$model, df_plot)
  }else{ annot_plot <- NULL }
  g <- ggplot(df_plot, myaes) + geom_point(shape = 20, size = ptsize)
  if(add_line){
    couls <- colorRampPalette(c("white", blues9))(256)
    g <- g + scale_fill_gradientn(colours = couls, breaks = NULL) +
      stat_smooth(method = "lm", col = "black", formula = , se = FALSE) +
      theme_classic() + theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_blank())
  }
  # myrange <- range(df_plot[, c("X1", "X2")])
  # g <- g + xlim(myrange) + ylim(myrange)
  if(diff(range(df_plot$X1)) > maxdiff && log2t == 'auto') g <- g + scale_x_log10()
  if(diff(range(df_plot$X2)) > maxdiff && log2t == 'auto') g <- g + scale_y_log10()
  if(!is.null(addvar)) addvar <- newlines(stringy = addvar, ln = 1, sepchar = "_")
  g <- g + labs(
    title = gsub("Log10_", "", paste0(myvars_names[1], " vs ", myvars_names[2])),
    subtitle = annot_plot, x = myvars_names[1], y = myvars_names[2], color = addvar
  )
  out_path <- ifelse(grepl("_$", out_path), out_path, dircheck(out_path))
  if(sum(!grepl("^_", suffix))) suffix <- paste0('_', suffix)
  if(mthod != "spearman") suffix <- paste0(suffix, '_', mthod)
  if(isTRUE(return_plot)) return(g)
  fname <- paste0(out_path, "correlation_", myvars_names[1], "_", myvars_names[2], suffix, ".pdf")
  if(v) cat("File:", fname, "\n")
  pdf(fname)
  print(g)
  dev.off()
}

# scatter plot from UpSetR
scatter_plot <- function(mydata, x, y){
  tmydata <- mydata[duplicated(mydata$gene_name), ]
  tvar <- unique(c(head(order(tmydata[, x], decreasing = TRUE), 5),
    head(order(tmydata[, y], decreasing = TRUE), 5)))
  subdata <- tmydata[tvar, ]
  att_plot <- (ggplot(data = mydata, aes_string(x = x, y = y, colour = "color"))
               + geom_point(shape=16) + scale_color_identity()
               + theme(panel.background = element_rect(fill = "white"),
                       plot.title = element_text(vjust = 1.3),
                       panel.grid.minor = element_blank(),
                       panel.grid.major = element_blank(),
                       axis.title.y = element_text(vjust = 1.3, size = 8.3),
                       axis.title.x = element_text(size = 8.3),
                       plot.margin=unit(c(0.1,1,0.1,0.5), "cm"))
               + geom_text_repel(
                 data = subdata,
                 aes(label = gene_name)#,
                 # box.padding = unit(0.35, "lines"),
                 # point.padding = unit(0.3, "lines")
               ))
}

#### Specificity score #### ----------------------------------------------------
# Similarly, we used the average expression (after z-score transformation) of 4
# well-defined naïve markers (CCR7, TCF7, LEF1 and SELL) and 12 cytotoxicity
# associated genes (PRF1, IFNG, GNLY, NKG7, GZMB, GZMA, GZMH, KLRK1, KLRB1, KLRD1,
# CTSW, CST7) to define the naïveness score and cytotoxicity score for both CD8+
# and CD4+ T cells, respectively. After delineating the exhaustion, naïveness
# and cytotoxicity scores of each T cell along the trajectory, we used locally
# weighted scatterplot smoothing (LOESS) regression to fit the relationships
# between these scores with Monocle components
specore <- function(
  mat,
  v = FALSE
){
  if(v) cat('-----\nCalculating specificity score\n')
}

marginal_plot = function(
  x,
  y,
  group = NULL,
  data = NULL,
  lm_show = FALSE,
  lm_formula = y ~ x,
  bw = "nrd0",
  adjust = 1,
  alpha = 1,
  plot_legend = T,
  group_colors = NULL,
  ...
){
  # NOTE: I removed "deparse(substitute())" managed variable
  ###############
  # Plots a scatterplot with marginal probability density functions for x and y.
  # Data may be grouped or ungrouped.
  # For each group, a linear fit can be plotted. It is hidden by default, but can be shown by providing lm_show = TRUE.
  # The model can be modified using the 'lm_formula' argument.
  # The 'bw' and 'adjust' argument specify the granularity used for estimating probability density functions. See ?density for more information.
  # For large datasets, opacity may be decreased by setting alpha to a value between 0 and 1.
  # Additional graphical parameters are passed to the main plot, so you can customize axis labels, titles etc.
  ###############
  moreargs = eval(substitute(list(...)))

  # prepare consistent df
  if(missing(group)){
    if(missing(data)){
      if(length(x) != length(y)){stop("Length of arguments not equal")}
      data = data.frame(x = as.numeric(x), y = as.numeric(y))
    } else {
      data = data.frame(x = as.numeric(data[, x]),
                        y = as.numeric(data[, y]))
    }
    if(sum(!complete.cases(data)) > 0){
      warning(sprintf("Removed %i rows with missing data", sum(!complete.cases(data))))
      data = data[complete.cases(data),]
    }
    group_colors = "black"
  } else {
    if(missing(data)){
      if(length(x) != length(y) | length(x) != length(group)){stop("Length of arguments not equal")}
      data = data.frame(x = as.numeric(x), y = as.numeric(y), group = as.factor(group))
    } else {
      data = data.frame(x = as.numeric(data[, x]),
                        y = as.numeric(data[, y]),
                        group = as.factor(data[, group]))
    }
    if(sum(!complete.cases(data)) > 0){
      warning(sprintf("Removed %i rows with missing data", sum(!complete.cases(data))))
      data = data[complete.cases(data),]
    }
    data = subset(data, group %in% names(which(table(data$group) > 5)))
    data$group = droplevels(data$group)
    if(is.null(group_colors)) group_colors = rainbow(length(unique(data$group)))
  }

  # log-transform data (this is need for correct plotting of density functions)
  if(!is.null(moreargs$log)){
    if(!moreargs$log %in% c("y", "x", "yx", "xy")){
      warning("Ignoring invalid 'log' argument. Use 'y', 'x', 'yx' or 'xy.")
    } else {
      data = data[apply(data[unlist(strsplit(moreargs$log, ""))], 1, function(x) !any(x <= 0)), ]
      data[,unlist(strsplit(moreargs$log, ""))] = log10(data[,unlist(strsplit(moreargs$log, ""))])
    }
    moreargs$log = NULL # remove to prevent double logarithm when plotting
  }

  # Catch unwanted user inputs
  if(!is.null(moreargs$col)){moreargs$col = NULL}
  if(!is.null(moreargs$type)){moreargs$type = "p"}

  # get some default plotting arguments
  if(is.null(moreargs$xlim)){moreargs$xlim = range(data$x)}
  if(is.null(moreargs$ylim)){moreargs$ylim = range(data$y)}
  if(is.null(moreargs$xlab)){moreargs$xlab = x}
  if(is.null(moreargs$ylab)){moreargs$ylab = y}
  if(is.null(moreargs$las)){moreargs$las = 1}

  # plotting
  tryCatch(expr = {
    ifelse(!is.null(data$group), data_split <- split(data, data$group), data_split <- list(data))
    orig_par = par(no.readonly = T)
    par(mar = c(0.25,5,1,0))
    layout(matrix(1:4, nrow = 2, byrow = T), widths = c(10,3), heights = c(3,10))

    # upper density plot
    plot(NULL, type = "n", xlim = moreargs$xlim, ylab = "density",
         ylim = c(0, max(sapply(data_split, function(group_set) max(density(group_set$x, bw = bw)$y)))), main = NA, axes = F)
    axis(2, las = 1)
    mapply(function(group_set, group_color){lines(density(group_set$x, bw = bw, adjust = adjust), col = group_color, lwd = 2)}, data_split, group_colors)

    # legend
    par(mar = c(0.25,0.25,0,0))
    plot.new()
    if(!missing(group) & plot_legend){
      legend("center", levels(data$group), fill = group_colors, border = group_colors, bty = "n", title = group, title.adj = 0.1)
    }

    # main plot
    par(mar = c(4,5,0,0))
    if(missing(group)){
      do.call(plot, c(list(x = quote(data$x), y = quote(data$y), col = quote(scales::alpha("black", alpha))), moreargs))
    } else {
      do.call(plot, c(list(x = quote(data$x), y = quote(data$y), col = quote(scales::alpha(group_colors[data$group], alpha))), moreargs))
    }
    axis(3, labels = F, tck = 0.01)
    axis(4, labels = F, tck = 0.01)
    box()

    if(lm_show == TRUE & !is.null(lm_formula)){
      mapply(function(group_set, group_color){
        lm_tmp = lm(lm_formula, data = group_set)
        x_coords = seq(min(group_set$x), max(group_set$x), length.out = 100)
        y_coords = predict(lm_tmp, newdata = data.frame(x = x_coords))
        lines(x = x_coords, y = y_coords, col = group_color, lwd = 2.5)
      }, data_split, rgb(t(ceiling(col2rgb(group_colors)*0.8)), maxColorValue = 255))
    }

    # right density plot
    par(mar = c(4,0.25,0,1))
    plot(NULL, type = "n", ylim = moreargs$ylim, xlim = c(0, max(sapply(data_split, function(group_set) max(density(group_set$y, bw = bw)$y)))), main = NA, axes = F, xlab = "density")
    mapply(function(group_set, group_color){lines(x = density(group_set$y, bw = bw, adjust = adjust)$y, y = density(group_set$y, bw = bw)$x, col = group_color, lwd = 2)}, data_split, group_colors)
    axis(1)
  }, finally = {
    par(orig_par)
  })
}

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
get_densities <- function(
  mat,
  genes,
  log2t = FALSE,
  subname = NULL,
  titl = NULL,
  cuof = 0,
  pdist = 1, # distance for percentages
  usedp = FALSE, # use only double positive
  even_axis = TRUE,
  noncero = FALSE,
  ptsize = 1,
  ptcolor = NULL,
  ptshape = 19,
  metas = NULL,
  groupby = NULL,
  return_plot = FALSE,
  v = FALSE
){
  vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

  if(!isTRUE(grepl("\\.pdf", subname))){ # if not PDF file name defined
    subname <- paste0(subname, ifelse(!grepl("_$", subname) && !is.null(subname), "_", ""))
    fname <- paste0(subname, paste0(genes, collapse = "_"), '_density.pdf')
  }else fname <- subname
  if(v && !isTRUE(return_plot)) cat("File:", fname, "\n")
  if(file.exists(fname) && !isTRUE(return_plot)) return(NULL)

  # mydf <- data.frame(t(mat[genes, ]))
  if(any(genes %in% rownames(mat))){
    mydf <- data.frame(t(mat[genes[genes %in% rownames(mat)], , drop = FALSE]))
  }else{ mydf = data.frame(row.names = rownames(mat)) }
  if(any(genes %in% colnames(metas))){
    mydf2 <- metas[rownames(mydf), genes[genes %in% colnames(metas)], drop = FALSE]
  }
  if(!exists("mydf")) mydf <- mydf2 else if(exists("mydf2")) mydf <- cbind(mydf, mydf2)
  mydf <- mydf[, genes]
  mydf <- mydf[complete.cases(mydf), ]
  # mydf[is.na(mydf)] <- 0
  cuof <- transformations(cuof, genes)
  dpos <- if(isTRUE(usedp)) mydf[, 1] > cuof[1] & mydf[, 2] > cuof[2] else rep(TRUE, nrow(mydf))
  if(sum(log2t) > 0){
    log2t <- transformations(log2t, genes)
    mydf[, log2t > 0] <- log2(mydf[, log2t > 0] + 1); cuof[log2t > 0] <- log2(cuof[log2t > 0] + 1)
  }
  couls <- c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
  # vls <- c(0, 2, 4, 6, 8, 10, 12)
  if(length(genes) < 3){
    ndens <- min(c(nrow(mydf), 100))
    tvar <- try(get_density(mydf[, 1], mydf[, 2], n = ndens), silent = TRUE)
    if(class(tvar) == 'try-error' || isTRUE(usedp)){
      if(v) cat("Double positive density\n")
      tvar <- rep(0, nrow(mydf))
      tvar[dpos] <- try(get_density(mydf[dpos, 1], mydf[dpos, 2], n = ndens), silent = TRUE)
    }
    if(class(tvar) != 'try-error'){
      if(v) cat("Using density\n")
      mydf$Density <- tvar; genes[3] <- 'Density'
    }else{
      if(v) cat("Using means\n")
      mydf$Means <- rowMeans(mydf[, 1:2]); genes[3] <- 'Means'
    }
  }

  # calculate the percentage in each quadrant
  qposition <- t(combn(sapply(mydf[, 1:2], range), 2)[, c(6:4, 2)])
  qposition <- data.frame(qposition)# + scale(qposition))
  rownames(qposition) <- c('lt',  'rt', 'rb', 'lb')
  colnames(qposition) <- c('x', 'y')
  qposition$d <- round(c(
    np = sum(mydf[, 1] <= cuof[1] & mydf[, 2] > cuof[2]),
    pp = sum(mydf[, 1] > cuof[1] & mydf[, 2] > cuof[2]),
    pn = sum(mydf[, 1] > cuof[1] & mydf[, 2] <= cuof[2]),
    nn = sum(mydf[, 1] <= cuof[1] & mydf[, 2] <= cuof[2])
  ) / nrow(mydf) * 100, 2)
  qposition$percent <- paste0(qposition$d, "%")
  even_axis <- if(isTRUE(even_axis)) list(c(1:2), c(1:2)) else c(1, 2)
  outergrid <- c(
    xmin = min(mydf[, 1]) - pdist, xmax = max(mydf[, even_axis[[1]]]) + pdist,
    ymin = min(mydf[, 2]) - pdist, ymax = max(mydf[, even_axis[[2]]]) + pdist
  )
  qposition['lt', 1:2] <- outergrid[c('xmin', 'ymax')]
  qposition['lb', 1:2] <- outergrid[c('xmin', 'ymin')]
  qposition['rt', 1:2] <- outergrid[c('xmax', 'ymax')]
  qposition['rb', 1:2] <- outergrid[c('xmax', 'ymin')]

  # if(v) print(head(mydf))
  genes <- colnames(mydf)
  if(v) cat("Plotting:", commas(genes), "\n")
  str(mydf)

  if(v) cat("Fitting\n")
  fit <- lm(formula = as.formula(paste(genes[1], "~", genes[2])) , data = mydf)#[mydf[, 1] > cuof[1] & mydf[, 2] > cuof[2], ])
  coefs <- summary(fit)$coefficients
  annot_plot <- paste0("r^2 = ", round( summary(fit)$r.squared, 2),
    "\np-value = ", coefs[nrow(coefs), ncol(coefs)],
    "\nSpearman (rho) = ", round( cor(mydf[, genes[1]], mydf[, genes[2]], method = "spearman"), 2))
  if(v) cat(annot_plot, '\n')
  clab <- if(!is.na(genes[3])) newlines(stringy = genes[3], ln = 1, sepchar = "_") else NULL
  p <- ggplot(mydf, aes_string(genes[1], genes[2])) +
    geom_density2d(data = mydf[dpos, ], colour = "#c0c5ce") +
    geom_point(aes_string(color = ifelse(is.null(ptcolor), genes[3], ptcolor)), size = ptsize, shape = ptshape) +
    xlim(outergrid[1:2]) + ylim(outergrid[3:4]) +
    geom_text(data = qposition, aes(x = x, y = y, label = percent)) +
    scale_color_gradientn(colours = couls, name = clab) +#viridis::scale_color_viridis(option = "magma") +
    labs(title = titl, subtitle = annot_plot, color = clab)
  if(sum(cuof > 0)){ p <- p + geom_hline(yintercept = cuof[2], linetype = "dashed", color = "gray") +
    geom_vline(xintercept = cuof[1], linetype = "dashed", color = "gray") }

  vlp <- try(vlnplot(
    cpmdata = t(mydf), metadata = metas[, !colnames(metas) %in% colnames(mydf)], gg = genes[!grepl("Density|Means", genes)],
    legendt = "posiv", orderby = groupby, plotdots = TRUE, v = v, datatype = 'log2(TPM + 1)',
    return_plot = TRUE, noncero = noncero, cuof = cuof, mysca = 'free'
  ) + theme(legend.position = "none"))

  if(isTRUE(return_plot)) return(list(violin = vlp, scatter = p))
  pdf(fname, height = 10, width = 14)
  print(plot_grid(vlp, p, rel_widths = c(1, 2)))
  # void <- try(plot.density2D(mydf[, 1:2], peaks = rownames(mydf)), silent = TRUE)
  # if(class(void) == 'try-error') warning('Density 2D not plotted') else title(titl)
  # print(plot_grid(xplot, NULL, p, yplot, ncol = 2, align = "hv", rel_widths = c(4, 1), rel_heights = c(1, 4)))
  # plot.new(); print(ggMarginal(p))
  # marginal_plot(x = genes[1], y = genes[2], data = mydf)
  graphics.off()
}

plot.density2D <- function(x, peaks = NULL){
  if(is.null(colnames(x))) colnames(x) <- paste0("Dim", 1:2)
  .density <- ks::kde(x)
  par(mar = c(5, 4, 5, 6))
  .zz <- seq(10, 90, 10)
  dnames <- colnames(x)

  plot(.density, display = "filled.contour2", cont = .zz, xlab = dnames[1], ylab = dnames[2])
  points(x[peaks, ,drop = FALSE], pch = 3, cex = 0.5)
  fields::image.plot(zlim = c(0, .zz[length(.zz)]), legend.only = TRUE,
    col = c("transparent", rev(heat.colors(length(.zz)))), axis.args = list(at = .zz,
      labels = sprintf("%s%%", 100 - .zz)), legend.width = 2.0, legend.mar = 4.5)
  return(.density)
}

# count the number of shared features (clone IDs) within a group against others
count_overlap <- function(
  x, # vector of c('groups column', 'of this group', 'groups' 'to compare' 'to', 'features column')
  annot,
  combs = FALSE
){
  clonas <- list(annot[annot[, x[1]] == x[2], x[5]]); names(clonas) <-  x[2]
  tvar <- sapply(x[-c(1:2, length(x))], function(y) annot[annot[, x[1]] == y, x[length(x)]] )
  tmp <- overlap_calc(c(clonas, tvar))
  if(combs == 'r') return(tmp)
  if(combs) return(sapply(tmp, length))
  sapply(tvar, function(y) sum(clonas[[1]] %in% y) )
}

share_counter <- function(
  selection = NULL,
  doby = 1,
  annot,
  name = NULL,
  combs = FALSE
){
  if(is.null(selection)) selection <- c(colnames(annot)[2], unique(annot[, 2]), colnames(annot)[3])
  void <- t(sapply(unique(annot[, doby]), function(x){
    tmeta <- annot[annot[, doby] == x, ]
    count_overlap(selection, tmeta, combs = combs)
  }))
  void <- rbind(void, ALL = count_overlap(selection, annot, combs = combs))
  void <- data.table::melt(void)
  fname <- paste0(tail(selection, 1), "_from_", selection[2], '_in_', paste0(selection[-c(1:2, length(selection))], collapse= "_"))
  fname <- paste0(ifelse(!is.null(name), paste0(name, "_"), ""), fname)
  if(combs) fname <- paste0(fname, "_expanded")
  pdf(paste0(fname, '.pdf'), 10, 10)
  print(ggplot(void, aes(x = Var2, y = value, fill = Var2)) +
    geom_bar(stat = "identity", alpha = 0.7) + rremove('legend') +
    facet_wrap(~ Var1, scale = 'free_y') +
    theme(axis.text.x = element_text(angle = 90)))
  graphics.off()
  return(void)
}

freq_tablep <- function(metadata, cnames, pnames = NULL, dowrite = FALSE){
  com <- metadata[, cnames[1]]
  meta <- metadata[, cnames[2]]
  freq_table <- list(prop.table(x = table(com, meta), margin = 2), prop.table(x = table(meta, com), margin = 2))
  cnames <- sub("integrated_snn_|orig\\.", "", cnames)
  if(is.null(pnames)){
    pnames <- paste(casefold(cnames, T), collapse = ' in ')
    pnames[2] <- paste(casefold(rev(cnames), T), collapse = ' in ')
  }
  names(freq_table) <- pnames
  if(isTRUE(dowrite) || is.character(dowrite)){
    dowrite <- if(is.character(dowrite)) dowrite else "proportion"
    fname <- paste0(c(dowrite, cnames, '.pdf'), collapse = "_")
    pdf(fname, 10, 8)
    for(i in 1:length(freq_table)){
      par(mfrow=c(1, 1), mar=c(5, 5, 4, 4))
      barplot(height = freq_table[[i]], col = v2cols(rownames(freq_table[[i]])), las = 2, legend.text = T,
        xlim=c(0, ncol(freq_table[[i]]) + 3), args.legend = list(x="topright", bty = "n", inset=c(-0.05, 0))
      ); title(names(freq_table)[i])
    }; graphics.off()
  }
  tvar <- table(com, meta)
  freqt <- as.matrix.data.frame(tvar)
  dimnames(freqt) <- dimnames(tvar)
  if(isTRUE(dowrite) || is.character(dowrite)) write.csv(freqt, file = sub("\\.pdf$", ".csv", fname))
}

#' Generate cumulative percentage of highly variable genes
#'
#' This function creates a plot with the cumulative percentage of  the variance/
#' dispersion and the number of genes in order to visualise a plateau
#'
#' @param object Seurat object or data.frame with HVG information (need a
#' logical column 'variable' specifying the if a genes is a highly variable gene)
#' @param smethod Method used, 'mvp', 'vst', etc.
#' @param cname Column with variance values.
#' @param cname_ord Column to order genes with, it usually takes 'cname.'
#' @param cutoff Mean threshold. It can also be a named (mean column).
#' @param frac_cutoff Fraction threshold. It can also be a named (pct column).
#' @param v Verbose
#'
#' @return Returns a plot with the variance explained by the HVGs and generates
#' and the table with the statistics
#'
#' @importFrom
#'
#' @export
#'
#' @examples
#' plot_df <- pct_variance(object = seurat_object, v = TRUE)
#'

cum_pct <- pct_variance <- function(
  object,
  smethod = NULL,
  cname = NULL,
  cname_ord = NULL,
  cutoff = NULL,
  frac_cutoff = NULL,
  v = FALSE
){
  if(class(object) == "Seurat"){
    if(is.null(smethod)){
      smethod <- unique(sub("\\..*", "", names(object@assays$RNA@meta.features)))
    }
    hdf <- HVFInfo(object, selection.method = smethod)
    hdf$variable <- rownames(hdf) %in% VariableFeatures(object)
  }else if(class(object) == "seurat"){
    hdf <- object@hvg.info; hdf$variable <- rownames(hdf) %in% object@var.genes
    smethod <- "mvp"; # cname <- "gene.dispersion"; cname_ord <- "gene.dispersion.scaled"
    object <- UpdateSeuratObject(object)
  }else{
    hdf <- object
  }; smethod <- paste0(smethod, collapse = ", ")
  if(v) cat("Total number genes:", nrow(hdf), "\n")
  if(v) cat("Methods:", ifelse(is.null(smethod), "none", smethod), "\n")
  colnames(hdf) <- sub(".*mean.*", "mean", colnames(hdf)) # re-naming mean column
  if(is.null(cutoff)) cutoff <- min(hdf$mean, na.rm = TRUE)
  # selecting the dispersion/variance column
  if(is.null(cname)) cname <- tail(grep("dispersion|variance", colnames(hdf)), 1)
  if(is.numeric(cname)) cname <- colnames(hdf)[cname]
  if(isTRUE(grepl("scaled", cname))){
    hdf[, paste0("shifted.", cname)] <- hdf[, cname] + abs(min(hdf[, cname]))
    cname <- paste0("shifted.", cname)
  }
  if(is.null(cname_ord)) cname_ord <- cname # column for ordering
  if(v) cat("Dispersion/variance column:", cname, "\n")
  if(v) cat("Ordering by:", cname_ord, "\n") # ' is added to avoid Excel problems
  hdf <- cbind(gene_name = paste0("'", rownames(hdf)), hdf)

  hdf <- hdft <- hdf[order(-hdf[, cname_ord]), ] # ordering by dispersion/variance
  if(v) print(head(hdf[hdf$variable, -1, drop = FALSE], 20)); #print(summary(hdf[, ])) }
  # % of variance out of the TOTAL variance
  # write.table(hdf, file = 'hvGenes.csv', sep = ',', row.names = FALSE)
  hvg_pct_total <- round(sum(hdf[hdf$variable, cname]) / sum(hdf[, cname], na.rm = TRUE) * 100, 2)
  hdf$cumulative_total = cumsum(hdf[, cname]) # cumulative (and %) across ALL genes
  hdf$cumulative_total_pct = round(hdf$cumulative_total / sum(hdf[, cname], na.rm = TRUE) * 100, 2)
  hdf$N_total = 1:nrow(hdf); hdf$lmean <- log2(hdf$mean + 1); hdf$lmean[!hdf$variable] <- NA
  if(!is.null(frac_cutoff) && class(object) == "Seurat"){
    hdf$exprFrac <- get_stat_report(mat = object@assays$RNA@counts, rnames = rownames(hdf), moments = "p", v = TRUE)
    hdf <- hdf[hdf$exprFrac > (frac_cutoff), ]
  }else if(!is.null(frac_cutoff)){
    hdf <- hdf[hdf[, names(frac_cutoff)] > (frac_cutoff), ]
  }
  hdf <- hdf[hdf$mean > cutoff, ] # FILTERING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  hdf$cumulative = cumsum(hdf[, cname]) # cumulative (and %) with genes avobe a threshold of mean expression
  hdf$cumulative_pct = round(hdf$cumulative / sum(hdf[, cname], na.rm = TRUE) * 100, 2)
  hdf$N = 1:nrow(hdf); hdft <- cbind(hdft, hdf[rownames(hdft), !colnames(hdf) %in% colnames(hdft)])
  # write.table(hdf, file = 'hvGenes_pcts.csv', sep = ',', row.names = FALSE)
  hvg_pct <- round(sum(hdf[hdf$variable, cname]) / sum(hdf[, cname], na.rm = TRUE) * 100, 2)
  hvg_label <- paste0(
    "From all ", nrow(hdft), " genes: ",  hvg_pct_total, "% (", sum(hdft$variable),
    " HVGs)\nFrom ", nrow(hdf), " with mean > ", cutoff, ": ", hvg_pct, "% (",
    sum(hdf$variable), " HVGs)\nMethod: ", smethod
  )
  if(any(head(!hdf$variable, sum(hdf$variable)))){
    capt <- paste("Warning: the 'top' genes may not be selected based on'", cname,
    "'\nSeurat v2 dispersions tends to have this effect\nYou could've applied another filter\n")
  }else{ capt <- NULL }
  if(v) cat(hvg_label, "\n")
  coulsn = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000')
  hvg_pct <- max(hdf$cumulative_pct[hdf$variable]) # last HGV genes
  ngenes <- nrow(hdf[hdf$cumulative_pct <= hvg_pct, ]) # max($N)
  p <- ggplot(hdf, aes(x = N, y = cumulative_pct, color = lmean)) + geom_point() +
    labs(
      title = "Cumulative variance of Highly Variable Genes",
      x = paste("Number of genes, order:", cname_ord),
      y = paste("Cumulative percentage:", gsub("\\.", " ", cname)),
      subtitle = hvg_label, color = "Mean", caption = capt
    ) + ylim(c(0, 100)) +
    scale_colour_gradientn(colours = coulsn, limits = c(0, max(hdf$lmean))) + theme_classic() +
    geom_hline(aes(yintercept = hvg_pct), linetype = "dashed", color = "red", size = 1.3) +
    geom_text(
      x = max(hdf$N) / 2, y = hvg_pct - 2.5,
      label = paste(ngenes, "genes = ", hvg_pct, "%"),
      size = 7, color = "black"
    )
  # if(diff(range(hdf$N)) > 1e3) p <- p + scale_x_log10()
  return(list(plot = p, table = hdft))
}

# plot per groups, number of cells
bar_dis <- function(
  annot,
  cnames = 1:3,
  cnorm = TRUE,
  subsamp = NULL,
  test = "poiwald",
  filt = TRUE,
  cols = NULL,
  plotdots = TRUE,
  ncolp = NULL,
  v = TRUE
){ # require ggsignif
  if(is.numeric(cnames)) cnames <- colnames(ddf)[cnames]
  if(v) cat("Variables:", commas(cnames), "\n")
  ddf <- data.table::melt(table(annot[, cnames]))
  ddf <- ddf[ddf$value > 0, ]
  yname <- casefold(sub("orig\\.", "", cnames[1]), upper = TRUE)
  xname <- casefold(sub("orig\\.", "", cnames[2]), upper = TRUE)
  if(!is.null(cnorm)){
    if(!is.character(cnorm)) cnorm <- cnames[1]
    if(v) cat("Percentage to:", cnorm, "\n")
    ddf$nvalue <- as.vector(ddf$value / table(annot[, cnorm])[as.character(ddf[, cnorm])])
    ddf$nvalue <- round(ddf$nvalue * 100, 1); myvalue <- 'nvalue'; yname <- paste0("%", yname)
  }else{ myvalue <- 'value'; yname <- "Number of cells"}
  rmp <- table(annot[, cnames[1]]); thresh <- mean(rmp) * 0.05 # filtering
  if(isTRUE(filt)){
    if(v) cat("Filtering <", round(thresh), "\n")
    tvar <- ddf[, cnames[1]] %in% names(rmp[rmp >= thresh])
    # tvar <- ddf[, 'value'] >= thresh
    if(v) cat("Out:", sum(!tvar), "\n")
    ddf <- ddf[tvar, ]
  }
  ddf$Type <- if(!is.na(cnames[3])) ddf[, cnames[3]] else factor("GROUP")
  tddf <- if(!is.null(subsamp)) ddf[getsubset(c(subsamp), ddf, v = v), ] else ddf
  ddf <- remove.factors(ddf)
  mypairs <- create_pairs(ddf[, cnames[2], drop = FALSE])[, 1:2]
  if(sum(test == "poiwald")){
    if(v) cat("Performing pairwise Wald test\n")
    ddf <- cbind(ddf, total = ave(ddf$value, as.character(ddf[, cnames[1]]), FUN = sum))
    ddf$num_other_cells <- ddf$total - ddf$value; #str(ddf)
    myf <- paste("value ~", cnames[2]); if(v) cat("Formula:", myf, "\n")
    ltype_tests <- lapply(unique(ddf$Type), function(x){
      d <- ddf[ddf$Type == x, ]
      group_tests <- apply(mypairs, 1, function(y){
        if(!all(unlist(y) %in% d[, cnames[2]])) return(1) # if only one group
        dtest <- d[d[, cnames[2]] %in% unlist(y), ]
        psn = glm(formula = as.formula(myf), family = poisson, offset = log(dtest$total), data = dtest)
        round(summary(psn)$coefficients[2, "Pr(>|z|)"], 6)
      })
      names(group_tests) <- apply(mypairs, 1, paste0, collapse = "vs")
      group_tests
    })
    names(ltype_tests) <- unique(ddf$Type); set.seed(27)
    type_tests <- round(unlist(ltype_tests), 7)
  }; set.seed(27)
  ncolp <- ifelse(is.null(ncolp), fitgrid(levels(tddf$Type))[2], ncolp)
  g1 <- ggplot(tddf, aes_string(x = cnames[2], y = myvalue, fill = cnames[2])) +
    stat_summary(fun.data = mean_se, geom = "errorbar", alpha = 0.7,
      width = 0.2, size = 1, aes_string(color = cnames[2]), show.legend = FALSE) +
    stat_summary(fun.y = mean, geom = "bar", size = 1.5, alpha = 0.7) +
    scale_fill_manual(values = v2cols(ddf[, cnames[2]], cols)) +
    scale_color_manual(values = v2cols(ddf[, cnames[2]], cols)) +
    theme_bw() + mytheme +
    theme(legend.position = "none", legend.title = element_blank(),
      panel.background = element_blank(), panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 17),
      strip.background = element_rect(fill = "#FFFFFF", linetype = 0),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12),
    plot.title = element_text(face = "bold", hjust = 0.5)) +
    labs(x = xname, y = yname) + facet_wrap(~ Type, scales = 'free_y', ncol = ncolp)
  if(plotdots) g1 <- g1 + geom_jitter(position = position_jitter(0.2), show.legend = FALSE)
  if(!is.null(test)){ # geom_signif or stat_compare_means
    if(v) cat("Test:", test, "\n")
    conames <- unname(as.list(data.frame(t(mypairs), stringsAsFactors = FALSE)))
    g1 <- g1 + geom_signif(
      test = ifelse(test != "poiwald", test, function(a, b) { list(p.value = 1) }),
      comparisons = conames, map_signif_level = FALSE,
      step_increase = 0.1#, textsize = 2, size = 0.2
    )
    pg <- ggplot_build(g1) #disassemble plot and obtain information
    slevel <- c("****" = 0.0001, "***" = 0.001, "**" = 0.01, "*" = 0.05, "ns" = 1) # less than
    aindex <- which(sapply(pg$data, function(x) 'annotation' %in% colnames(x) ))
    pvdf <- remove.factors(pg$data[[aindex]][seq(1, nrow(pg$data[[aindex]]), 3), c("group", "annotation", "PANEL")])
    pvdf$group <- sub("\\-.*", "", sub("\\-", "vs", pvdf$group))
    pvdf[, ifelse(is.na(cnames[3]), "Type", cnames[3])] <- levels(tddf$Type)[as.numeric(pvdf[, 3])]
    if(test == "poiwald"){
      pvdf$pval <- type_tests[paste0(pvdf[, 4], ".", pvdf[, 1])]
    }else{ pvdf$pval <- as.numeric(pvdf$annotation) }; set.seed(27)
    pvdf$padj = round(p.adjust(pvdf$pval, method = "fdr"), 7)
    pvdf <- pvdf[, c(4, 1, 5, 6)]
    pvdf$Level <- sapply(as.numeric(pvdf$padj), function(x) head(names(slevel[slevel >= x]), 1)  )
    pvdf$Test <- ifelse(test == "poiwald", "Wald test", test)
    pg$data[[aindex]]$annotation <- factor(rep(pvdf$Level, each = 3))
    pdf(file = NULL); g1 <- ggplot_gtable(pg); dev.off() #reassemble the plot
    return(list(bars = g1, res = tableGrob(pvdf)))
  }
  return(g1)
}
