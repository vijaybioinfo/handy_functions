#!/usr/bin/R

##############################
# Clustering handy functions #
##############################

# This functions are designed to be help with Seurat object's hanndling and continues to expand

get_source_data <- function(
  xpath,
  pj = '10XData',
  merge_counts = FALSE,
  v = FALSE
) {
  if(v) cat("Path:", xpath, "\n")
  if(isTRUE(merge_counts) && file.exists(xpath)){
    if(dir.exists(xpath)){
      fname <- list.files(dirname(xpath), pattern = "aggregation")
    }else{ fname <- xpath }
    lnames <- read.csv(fname, stringsAsFactors = FALSE)
    lnames <- lnames[, apply(lnames, 2, function(x) all(file.exists(x)) )]
    lnames <- unname(sapply(lnames, function(x){
      if(!dir.exists(x)) return(paste0(dirname(x), "/filtered_feature_bc_matrix"))
      set_dir_pats(x)
    }))
    names(lnames) <- basename(dirnamen(lnames, 2))
    fnames <- lnames #gsub("beegfs/ciro", "BioAdHoc/Groups/vd-vijay/cramirez", lnames)
    if(v) cat(paste0(fnames, collapse = "\n"), "\n")
    fnames <- fnames[dir.exists(fnames)]
    if(v) cat("Merging libraries:", length(fnames), "/", length(lnames), "\n")
    tgottendata <- lapply(names(fnames), function(x){
      if(v) cat("."); y <- Seurat::Read10X(fnames[x])
      colnames(y) <- sub(paste0(x, "_"), "", paste0(colnames(y), "-", which(names(lnames) == x)))
      y
    }); if(v) cat(" ")
    allgenes <- unique(unlist(lapply(tgottendata, rownames)))
    if(v) cat("All genes:", commas(allgenes), "\n")
    allgenes <- allgenes[rowSums(sapply(tgottendata, function(x) allgenes %in% rownames(x) )) == length(tgottendata)]
    if(v) cat("Overlaping genes:", commas(allgenes), "\n")
    tgottendata <- lapply(tgottendata, function(x){
      x[allgenes, ]
    })
    gottendata <- do.call('cbind', tgottendata)
    rm(tgottendata)
  }else if(dir.exists(xpath)){
    xpath <- set_dir_pats(xpath, pj = pj, v = v)
    if(v) cat(' - 10X directory:', xpath, " ")
    gottendata <- Seurat::Read10X(data.dir = xpath)
  }else if(file.exists(xpath)){
    if(v) cat(' - Object... '); gottendata <- readfile(xpath, v = v)
    otype <- casefold(class(gottendata)); if(v) cat(otype, '\n');
    if(otype == 'seurat'){
      this_annot <- gottendata@meta.data
      gottendata <- gottendata[["RNA"]]@counts
    }else if(otype == 'saver'){
      if(v) cat(' - SAVER '); gottendata <- gottendata$estimate
    }; xpath <- dirname(xpath)
  }else{ stop("data not found\n") }
  if(v) cat('DONE!\n'); gc()
  if(exists('this_annot')) return(list(mycellsdata = gottendata, source_file = xpath, annottab = this_annot))
  return(list(mycellsdata = gottendata, source_file = xpath))
}

# Get top markers
get_top_n <- function(
  x,
  cname = 'cluster',
  dpct = 'Dpct',
  orderby = 'p_val_adj',
  patties = "^rps|^rpl|^mt-",
  filter_neg = TRUE,
  n = Inf,
  v = FALSE
) {
  if(v) cat("Getting top", n, "features per", cname, "\n")
  groups <- names(table(x[, cname]))
  if(is.null(dpct)) x[, dpct] <- 100 * (x$pct.1 - x$pct.2)
  y <- as.data.frame(data.table::rbindlist(lapply(groups, function(z){
    if(v) cat('.')
    dat <- x[which(as.character(x[, cname]) == z), ]
    if(isTRUE(filter_neg) && any(dat[, dpct] > 0)) dat <- dat[dat[, dpct] > 0, ]
    tvar <- grepl(patties, casefold(dat$gene))
    if(sum(tvar) && sum(!tvar)) dat <- dat[!tvar, ] # Removing RPS and RPL genes
    dat$ribomito <- sum(tvar)
    data.table::setorderv(dat, orderby)
    head(dat, n)
  }))); if(v) cat('\n')
  return(y)
}

ClassifyScoring <- function(
  object,
  features,
  set.ident = FALSE,
  name = 'Set',
  verbose = FALSE,
  ...
) {
  cc.columns <- grep(pattern = name, x = colnames(x = object[[]]), value = TRUE)
  if(is.null(names(features))) names(features) <- paste0("S", 1:length(features))
  if(name %in% cc.columns){ warning(name, " pre-computed"); return(object) }
  if(verbose) str(features)

  # Renaming columns colliding with previous signatures; first the 'name'
  cc.columns <- make.names(c(cc.columns, name), unique = TRUE); name <- cc.columns[length(cc.columns)]
  cc.columns <- make.names(c(colnames(x = object[[]]), names(features)), unique = TRUE); # now the 'classes'
  names(features) <- cc.columns[tail(1:length(cc.columns), length(names(features)))]

  classes <- names(features)
  if(verbose) cat("Calculating scores\n")
  ctrl_n <- min(vapply(X = features, FUN = length, FUN.VALUE = numeric(length = 1)))
  if(ctrl_n > 1000) ctrl_n <- 100
  if(verbose) cat("Background size:", ctrl_n, "\n")
  object <- AddModuleScore(
    object = object,
    features = features,
    name = name,
    ctrl = ctrl_n,
    ...
  )
  cc.columns <- grep(pattern = name, x = colnames(x = object[[]]), value = TRUE)
  cc.scores <- object[[cc.columns]]
  object@meta.data <- object@meta.data[, !colnames(object@meta.data) %in% cc.columns]
  # rm(object.cc); gc(verbose = FALSE) # 'object' was duplicated from AddModuleScore
  if(verbose) cat("Classification based on score.\nNone: all < 0; Undecided: max score > 1.\n")
  assignments <- apply(
    X = cc.scores,
    MARGIN = 1,
    FUN = function(scores) {
      if (all(scores < 0)) {
        return("None")
      } else {
        if (length(which(x = scores == max(scores))) > 1) {
          return('Undecided')
        } else {
          return(classes[which(x = scores == max(scores))])
        }
      }
    }
  )
  cc.scores <- merge(x = cc.scores, y = data.frame(assignments), by = 0)
  colnames(x = cc.scores) <- c('rownames', classes, name)
  rownames(x = cc.scores) <- cc.scores$rownames
  cc.scores <- cc.scores[, c(classes, name)]
  if(verbose) cat("Adding to meta data:", paste0(c(classes, name), collapse = ", "), "\n")
  object[[colnames(x = cc.scores)]] <- cc.scores
  if (set.ident) {
    if(verbose) cat("Setting classes as identities\n")
    object[['old.ident']] <- Idents(object = object)
    Idents(object = object) <- 'Class'
  }
  return(object)
}

# centered log-ratio (CLR) normalization - for each feature!
clr_function <- function(x) {
  apply(X = x, MARGIN = ifelse(nrow(x) < ncol(x), 1, 2), FUN = function(x){
    return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
  })
}

## Statistics and extra information
markers_summary <- function(
  marktab,
  annot,
  resolut = NULL,
  covar = NULL,
  datavis = NULL,
  datatype = "CTS",
  v = FALSE
){
  marktab <- remove.factors(marktab)
  annot <- remove.factors(annot)
  if(casefold(class(datavis)) == "seurat"){
    if(v) cat("Taking data from seurat object\n")
    datavis <- as.matrix(GetAssayData(datavis))
    datatype <- "SN"
  }
  if(is.null(resolut)) resolut <- colnames(annot)[grep("res", colnames(annot))[1]]
  if(v) cat("Cluster column:", resolut, "\n")

  if(v) cat("Number of clusters per gene\n")
  ncluster <- t(sapply(unique(cmarkers$gene), function(x){
    y <- marktab[marktab$gene == x, "cluster"]
    if(length(y) < 1) return(NA)
    c(nCluster = as.character(length(y)), sCluster = paste0(y, collapse = "&"))
  }))
  marktab <- cbind(marktab, ncluster[marktab$gene, ])

  stattab <- get_stat_report(
    mat = datavis,
    groups = make_list(x = annot, colname = resolut, grouping = TRUE),
    rnames = rownames(datavis),
    do_log = FALSE,
    indiv_values = FALSE,
    datatype = datatype,
    v = v
  )
  if(all(covar %in% colnames(annot)) && !is.null(covar)){
    if(v) cat("Using", commas(covar), "as covariate\n")
    annot$covar <- do.call('paste', c(annot[c(resolut, covar)], sep = "_"))
    stattab <- cbind(stattab, get_stat_report(
      mat = datavis,
      groups = make_list(x = annot, colname = "covar", grouping = TRUE),
      rnames = rownames(datavis),
      datatype = datatype,
      v = v
    ))
    grps <- unique(annot$covar)
  }
  marktab$resolution <- resolut
  marktab$nCluster <- length(unique(marktab$cluster))
  marktab$cluster <- as.character(marktab$cluster)
  marktab$nCells <- as.numeric(table(annot[, resolut])[as.character(marktab$cluster)])
  marktab$cluster_n <- try(paste(marktab$cluster, "- N:", marktab$nCells))
  marktab$cPct <- round(marktab$nCells / nrow(annot) * 100, 4)
  marktab$nMarkers <- as.numeric(table(marktab$cluster)[marktab$cluster])
  if(v) cat("Calculating distances within clusters\n")
  kha <- sapply(unique(marktab$cluster), function(i){ # mean distance per cluster
    thesecells <- getsubset(c(resolut, i), annot, v = FALSE)
    thesecells <- sample(thesecells, min(c(length(thesecells), 500)))
    mydata <- t(datavis[marktab[marktab$cluster == i, 'gene'], thesecells, drop = FALSE])
    tvar <- mean(dist(mydata), na.rm = TRUE)
  });
  marktab$dist <- kha[marktab$cluster]
  if(v) cat("Mean differences per gene\n")
  marktab$gene_mean <- sapply(1:nrow(marktab), function(i){ # mean per gene per cluster
    stattab[as.character(marktab[i, 'gene']), paste0(marktab[i, 'cluster'], '_mean', datatype)]
  }) # mean of the rest of clusters per gene per cluster
  marktab$bmean <- unlist(sapply(unique(marktab$cluster), function(i){
    thesecells <- getsubset(c(resolut, paste0('-', i)), annot, v = F)
    rowMeans(datavis[marktab[marktab$cluster == i, 'gene'], thesecells, drop = FALSE], na.rm = TRUE)
  }))
  marktab$Dmean <- (marktab$gene_mean - marktab$bmean)
  marktab$Dpct <- 100 * (marktab$pct.1 - marktab$pct.2)
  tvar <- table(marktab$cluster, marktab$Dpct < 0)
  marktab$nDpct_neg <- if(ncol(tvar) == 2) tvar[as.character(marktab$cluster), "TRUE"] else 0
  marktab$Significance <- (-log10(marktab$p_val_adj))
  tvar <- max(marktab$Significance[marktab$Significance != Inf])
  marktab$Significance[is.infinite(marktab$Significance)] <- tvar
  marktab <- cbind(marktab, stattab[marktab$gene, ])
  print(headmat(stattab))

  stattab$gene <- rownames(stattab)
  marktab$gene_name <- paste0("'", marktab$gene)
  stattab$gene_name <- paste0("'", rownames(stattab))
  marktab <- data.table::rbindlist(list(marktab, stattab[!rownames(stattab) %in% marktab$gene, ]), fill = TRUE)
  marktab <- data.frame(marktab, stringsAsFactors = FALSE, check.names = FALSE)
  rownames(marktab) <- make.names(marktab$gene, unique = TRUE)
  marktab
}

make_cluster_names <- function(vec){
  vec <- as.character(vec)
  level_names <- names(sort(table(vec), decreasing = TRUE))
  levels <- 0:(length(level_names) - 1)
  names(levels) <- level_names
  factor(x = levels[vec], levels = unname(levels))
}

add_reduction <- function(
  object,
  reduction,
  reduction.key = "Dim_",
  reduction.name = 'reduction'
){
  object@reductions[[reduction.name]] <- object@reductions[[1]]
  scells <- rownames(object@reductions[[1]]@cell.embeddings)
  object@reductions[[reduction.name]]@cell.embeddings <- as.matrix(reduction[scells, ])
  colnames(object@reductions[[reduction.name]]@cell.embeddings) <- paste0(reduction.key, 1:ncol(reduction))
  object@reductions[[reduction.name]]@key <- reduction.key
  return(object)
}
