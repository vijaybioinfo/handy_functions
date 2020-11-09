#!/usr/bin/R

################
# Cluster tree #
################

# This script will create plots from clustering analyses using 'clustree'
# package

# .libPaths('~/R/newer_packs_library/3.5')
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')
suppressPackageStartupMessages(library(clustree))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))

generate_tree <- function(
  dpath,
  patfile = 'metadata_',
  cprefix = 'RNA_snn_res.',
  libname = 'origlib',
  sobject = NULL,
  patty = NULL,
  exclude = NULL,
  outdir = NULL,
  myorder = NULL,
  v = FALSE
){
  if(v) cat("-- Building cluster tree --\n")
  ## Finding files ## --
  fnames <- if(!any(file.exists(dpath)) || dir.exists(dpath)){
    dnames <- sort(list.files(path = dpath, pattern = patty, full.names = TRUE))
    dnames <- gsub("\\/{2,}", "/", dnames)
    if(!is.null(exclude)) dnames <- dnames[!grepl(exclude, dnames)]
    sort(list.files(path = dnames, pattern = patfile, full.names = TRUE, recursive = TRUE))
  }else{ dpath }
  fnames <- fnames[file.exists(fnames)]
  if(!is.null(exclude)) fnames <- fnames[!grepl(exclude, fnames)]

  if(is.null(outdir)) outdir <- paste0(head(dnames, 1), '/clustree')
  outdir <- dircheck(outdir)
  # fnames <- fnames[order(file.info(fnames)$ctime)] # ordering file by time of creation

  if(!is.null(myorder)){
    if(v) cat("Ordering\n")
    fnames <- unlist(lapply(myorder, function(x) fnames[grepl(paste0("/", x, "/"), fnames)] ))
  }

  if(v){
    cat("Writing output to:", outdir, '\n')
    cat("Taking files:", fnames, sep = '\n')
  }; if(!dir.exists(outdir)) dir.create(outdir)

  ## Reading files ## --
  metas <- lapply(fnames, function(x){
    if(v) cat(".")
    mypat <- if(grepl('/markers/', x)) '.*kers.([0-9]{1,})PCs.*' else '.*ta_([0-9]{1,})PC.*'
    pcs <- gsub(mypat, '\\1', x)
    y <- remove.factors(readfile(x, row.names = 1, v = FALSE))
    if(all(c('cluster', 'gene') %in% colnames(y))){
      yy <- column_collapse(
        metab = y[, c('cluster', 'gene_name')],
        rname = 'gene_name',
        v = v
      )
      yy <- setrows(yy)
      yy$ngenes <- table(yy[, 'cluster'])[yy[, 'cluster']]
      splat <- strsplit(yy[, 'cluster'], "&")
      yy$nclusters <- sapply(splat, length)
      nlayers <- 1 #max(yy$nclusters)
      newy <- data.frame(sapply(1:nlayers, function(iclust){
        gsub(" {1,}", "", sapply(splat, "[", iclust))
      }), row.names = rownames(yy))
      colnames(newy) <- paste0(cprefix[1], 1:ncol(newy))
      cat(paste0("c", ncol(newy)))
      y <- cbind(yy, newy)
      # length(table(yy[yy$ngenes > 3, 'cluster'])) # too many "clusters"
      # twogroups <-  < 3
      # length(table(yy[twogroups, 'cluster'])) # too many "clusters"
    }else{
      if(libname %in% colnames(y)) rownames(y) <- paste0(gsub('\\-.*', '-', rownames(y)), y[, libname])
      y <- y[, grepl(cprefix[1], colnames(y)), drop = FALSE]
    }
    if(!is.na(cprefix[2])) y <- y[, grepl(cprefix[2], colnames(y)), drop = FALSE]
    colnames(y) <- gsub(cprefix[1], paste0(cprefix[1], pcs), colnames(y))
    tvar <- grep(cprefix[1], colnames(y), value = TRUE)
    if(any(grepl("\\.[0-9]{1,}\\.", tvar))) colnames(y) <- sub("\\.", "", colnames(y))
    return(y)
  }); if(v) cat("\n")
  # str(metas)
  head(metas[[1]])
  aremarkers <- all(c('nclusters', 'ngenes') %in% colnames(metas[[1]]))
  suffixf <- if(aremarkers) 'markers_' else ''
  # metas <- metas[sapply(metas, nrow) == max(sapply(metas, nrow))]

  # Adding number for IDs
  if(v) cat("Creating IDs\n")
  metas <- lapply(1:length(metas), function(x){
    y <- metas[[x]]
    y$barcode <- rownames(y)
    colnames(y) <- gsub(cprefix[1], paste0(cprefix[1], x, "0"), colnames(y))
    return(y)
  })
  # str(metas)

  summary_clust <- metas %>% Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = 'barcode'), .)
  summary_clust <- remove.factors(summary_clust)
  cnames <- colnames(summary_clust)
  summary_clust <- summary_clust[, c('barcode', cnames[cnames != 'barcode'])]
  rownames(summary_clust) <- summary_clust$barcode
  # str(summary_clust)
  if(aremarkers){
    summary_clust <- summary_clust[complete.cases(summary_clust), ]
  }
  summary_clust[is.na(summary_clust)] <- 'N'

  # tvar <- overlap_calc(lapply(metas, rownames), v = TRUE)
  # tvar <- tvar[sapply(tvar, length) > 0]
  redus <- NULL
  if(!is.null(sobject) && !aremarkers){
    if(v) cat("Using object\n")
    if(casefold(class(sobject)) != "seurat") sobject <- theObjectSavedIn(sobject)
    tvar <- paste0(gsub('\\-.*', '-', rownames(sobject@meta.data)), sobject@meta.data[, libname])
    if(any(rownames(summary_clust) %in% tvar)){
      redus <- list(tsne = c('tSNE_1', 'tSNE_2'), umap = c('UMAP_1', 'UMAP_2'))
      ddf <- FetchData(object = sobject, vars = c(unlist(redus), libname))
      redus <- redus[sapply(redus, function(x) all(x %in% colnames(ddf)) )]
      if(length(redus) > 0){
        rownames(ddf) <- paste0(gsub('\\-.*', '-', rownames(ddf)), ddf[, libname])
        missingcells <- rownames(summary_clust)[!rownames(summary_clust) %in% rownames(ddf)]
        ddf <- rbind(ddf, mat_names(missingcells, colnames(ddf)))
        summary_clust <- cbind_repcol(summary_clust, ddf)
      }
    }
  }
  if(v) str(summary_clust)

  # Clustering order
  summary_df <- data.frame(
    N = 1:length(fnames),
    Name = fnames,
    PCs = gsub('.*ta_([0-9]{1,})PC.*', '\\1', fnames)
  )
  write.table(summary_df, file = paste0(outdir, suffixf, 'cluster_order.txt'), quote = FALSE, sep = "\t", row.names = FALSE)
  write.table(summary_clust, file = paste0(outdir, suffixf, 'cluster_annnotation.txt'), sep = "\t", row.names = FALSE)

  pp <- list()
  if(v) cat("Building trees...\nLayers:", sum(grepl(cprefix[1], colnames(summary_clust))), "\n")
  pp[['tree']] <- clustree(summary_clust, prefix = cprefix[1], layout = )
  if(aremarkers){
    pp[['sugi_tree']] <- clustree(summary_clust, prefix = cprefix[1], layout = "sugiyama")
    pp[['sugi_nc_tree']] <- clustree(summary_clust, prefix = cprefix[1], layout = "sugiyama", use_core_edges = FALSE)
  }
  if(!aremarkers && (length(redus) > 0)){
    for(i in 1:length(redus)){
      pname <- paste0(names(redus[i]), '_overlay'); if(v) cat(pname)
      if(!all(redus[[i]] %in% colnames(summary_clust)) || pname %in% names(pp)){
        if(v) cat("done\n"); next
      }; if(v) cat('!\n')
      pp[[pname]] <- clustree_overlay(
        summary_clust,
        prefix = cprefix[1],
        x_value = redus[[i]][1], y_value = redus[[i]][2],
        label_nodes = FALSE
      )
    }
  }

  pp_dims <- lapply(gsub(".*_", "", names(pp)), function(x){
    if(grepl('^tree', x)) return(c(14, 14))
    if(grepl('^overlay', x)) return(c(14, 14))
  })

  if(v) cat("Plotting\n")
  for(i in 1:length(pp)){
    fname <- paste0(outdir, suffixf, names(pp[i]), '.pdf')
    # cat(fname, "\n")
    # if(file.exists(fname)) next
    pdf(fname, width = pp_dims[[i]][1], height = pp_dims[[i]][2]);
    print(pp[[i]])
    dev.off()
  }
}
