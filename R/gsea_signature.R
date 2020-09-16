#!/usr/bin/R

###################################
# Clustering analysis: signatures #
###################################

# This script generates signature scores

source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
if(!exists("vlab_signatures")) source('/mnt/BioHome/ciro/scripts/seurat/utilities.R')
# source('/mnt/BioHome/ciro/scripts/seurat/plotting.R')

clean_feature_list <- function(
  mat,
  features, # list of genes
  grps = NULL,
  filterby = NA, # it can also be percentage (p) or mean (mn) [and value] (eg. p2, default is 0)
  return_stats = FALSE,
  v = FALSE
){
  tvar <- gsub("[0-9]{1,}", "", filterby); tvar[2] <- gsub("[a-z]", "", filterby); if(tvar[2] == "") tvar[2] <- 0
  str(tvar)
  if(is.null(grps)){ grps <- rep("SingleCell", ncol(mycells)); names(grps) <- colnames(mycells) }
  mygenes <- unique(unlist(features)); mygenes <- getfound(mygenes, rownames(mat), v = v)
  void <- data.frame(get_stat_report(
    mat = as.matrix(mat[mygenes, ]),
    groups = grps,
    moments = tvar[1],
    v = v
  ), row.names = mygenes); if(v) print(summary(void))
  passed <- rownames(void)[void[, 1] > as.numeric(tvar[2])]
  if(v) str(features)
  features <- sapply(features, function(x) x[x %in% passed], simplify = FALSE)
  if(v) str(features)
  if(isTRUE(return_stats)) list(stats = void, list = features) else return(features)
}


signature_scoring <- function(
  object,
  prefix = NULL,
  lsignatures = list(name = "orig.Cell.Cycle", S.Score = cc.genes$s.genes, G2M.Score = cc.genes$g2m.genes),
  group = FALSE, # group signatures with same name [and description] (NAME_description_source)
  confounders = "RNA.*res",
  reductions = list(pca = c('PC_1', 'PC_2'), tsne = c('tSNE_1', 'tSNE_2'), umap = c('UMAP_1', 'UMAP_2')),
  couls = c("#fffeee", "#ffe080", "#ffc100", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000"),
  v = FALSE
){
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(cowplot))
  theme_set(theme_cowplot())
  str_safe_remove <- function(x, word = "random123") gsub("_{1,}", "_", gsub(word, "", x))

  if(v) cat("-- Signature scoring --\n")
  if(!grepl("signatures", prefix)) prefix <- dircheck(paste0(prefix, 'signatures'))
  if(!dir.exists(prefix) && grepl("\\/$", prefix)) dir.create(prefix, recursive = TRUE)
  if(v) cat("Output:", prefix, "\n")
  # initpath <- getwd()
  # setwdc(prefix, showWarnings = FALSE)
  if(any(!confounders %in% colnames(object@meta.data))){
    confounders <- get_grouping_cols(metadat = object@meta.data, ckeep = confounders, maxn = 56, v = v)
  };
  confounders <- confounders[confounders %in% colnames(object@meta.data)]
  if(v) cat("Confounders", commas(confounders), "\n");
  Idents(object) <- 'orig.ident'

  reductions <- reductions[names(reductions) %in% names(object@reductions)]

  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  signature_list <- list()
  if(length(lsignatures)){
    for(i in 1:length(lsignatures)){
      mynameis <- sub("[[:alnum:]]{1,8}_", "", names(lsignatures[i]))
      if(v && isTRUE(group)) cat("@", mynameis, "\n")
      if(isTRUE(group)){
        class_list <- list(
          lsignatures[grep(mynameis, names(lsignatures))]
        )
        names(class_list) <- mynameis
      }else{
        class_list <- lsignatures[i]
        class_list <- lapply(class_list, function(x) x[x %in% rownames(object)] )
      }
      tvar <- c(gsub("orig.", "", sapply(signature_list, "[[", 1)), names(signature_list))
      if(any(names(class_list) %in% tvar)) next
      print(sapply(class_list, length))
      str(class_list)
      if(any(sapply(class_list, length) == 0)){
        warning("'", mynameis, "' has 0 feature value(s)"); next
      }
      # names(class_list) <- paste0(casefold(gsub(paste0("_", mynameis), "", names(class_list)), upper = TRUE), ".Score")
      names(class_list) <- paste0(casefold(str_safe_remove(names(class_list), "signature"), upper = TRUE), ".Score")
      signature_list <- c(signature_list, list(c(name = paste0("orig.", names(lsignatures[i])), class_list)))
    }
  }
  signame <- paste0(prefix, "signatures.csv")
  if(file.exists(signame)){
    if(v) cat("Pre-computed scores\n")
    signdf <- readfile(signame, stringsAsFactor = FALSE, check.names = FALSE, row.names = 1)
    object@meta.data <- joindf(object@meta.data, signdf, v = v)
  }
  scores_cols <- unname(unlist(sapply(signature_list, function(x) c(names(x)[-1], x[[1]]) )))
  for(scoring in signature_list){
    if(v) cat("\n@", sub("orig\\.", "", scoring[[1]]), "\n")
    if(!scoring[[1]] %in% colnames(object@meta.data)){
      object <- ClassifyScoring(
        object = object, name = scoring[[1]], verbose = v,
        features = scoring[head(2:length(scoring), 2)]
      ); str(scoring)
      tvar <- FetchData(object, vars = scores_cols[scores_cols %in% colnames(object@meta.data)])
      if(exists("signdf")) tvar <- joindf(tvar, signdf)
      write.csv(tvar, file = signame)
    }

    if(v) cat("Plotting\n")
    ddfplot <- FetchData(object, vars = c(scoring[[1]], names(scoring[-1]), unname(unlist(reductions)), confounders))
    ddfplot$Signature <- as.character(ddfplot[, scoring[[1]]])
    ddfplot$Signature <- ifelse(ddfplot$Signature == "None", "No", "Yes")
    # for(i in names(scoring[-1])) ddfplot[, i] <- scales::rescale(ddfplot[, i], to = c(0, 1))
    tvar <- table(ddfplot[, 'Signature'])
    subtitl <- paste0(paste0(names(tvar), ": ", unname(tvar)), collapse = "; ")

    if(v) cat(" # heatmap\n")
    graphics.off()
    fname <- paste0(prefix, gsub("orig\\.", "", casefold(scoring[[1]])), '_heatmap.pdf')
    if(!finished_file(fname)){
      pdf(fname, width = 10, height = 12, onefile = FALSE);
      p <- try(custom_heatmap(
        object = object, # should be normalised
        rnames = unname(unlist(scoring[-1])),
        orderby = scoring[[1]],
        sample_it = c("orig.ident", '3000'), # can be c("column name", "-num_per_group"); or will use the first column
        scale_row = TRUE,
        categorical_col = c(confounders, "Signature"),
        feature_order = "pca",
        couls = NULL, # colors for columns and rows
        hcouls = c('yellow', 'black', 'blue'),
        regress = c('nCount_RNA', 'percent.mt'),
        verbose = FALSE
      ))
      graphics.off()
    }
    for(redu in names(reductions)){
      if(v) cat(" #", redu, "\n")
      fname <- paste0(prefix, gsub("orig\\.", "", casefold(scoring[[1]])), '_', redu, '.pdf')
      if(finished_file(fname)) next
      p <- lapply(names(scoring[-1]), function(x){
        aesy <- aes_string(x = reductions[[redu]][1], y = reductions[[redu]][2], color = x)
        ggplot(data = ddfplot, mapping = aesy) +
        geom_point(size = 0.1) + scale_color_gradientn(colours = couls) +
        labs(colour = NULL, title = make_title(x))
      })
      pdf(fname, width = ifelse(length(p) > 1, 11, 7), height = ifelse(length(p) > 2, 11, 7));
      print(cowplot::plot_grid(plotlist = p)); graphics.off()
    }
    for(confy in confounders){
      if(v) cat(" -", confy, "\n")
      fname <- paste0(prefix, gsub("orig\\.", "", casefold(scoring[[1]])), "_violin_", confy, '.pdf')
      if(finished_file(fname)) next
      tvar <- length(table(ddfplot[, confy]))
      ddfplot[, confy] <- factormix(ddfplot[, confy])
      p <- violins(
        dat = ddfplot,
        xax = confy,
        yax = names(scoring[-1]),
        colour_by = "pct"
      ) + RotatedAxis()
      pdf(fname, width = ifelse(tvar > 4, 11, 7), height = ifelse(length(scoring) > 2, 10, 7));
      print(p)
      graphics.off()
    }
  }
  # setwd(initpath)
  if(v) cat("-- --------- ------- --\n")
  return(0) # return(object)
}
