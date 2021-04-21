#!/usr/bin/R

###################################
# Clustering analysis: signatures #
###################################

# This script generates signature scores

# "devel/utilities.R", "devel/plots.R", "R/stats_summary_table.R", "clustering_utilities.R"

clean_feature_list <- function(
  mat,
  features, # list of genes
  groups = NULL,
  filterby = 'p~0', # Moment~[value] (default is p~0) as indicated in `stats_summary_table`. You can include the group name to filter by, p~0~group1.
  return_stats = FALSE,
  verbose = FALSE
){
  if(verbose) cat("\n")
  filterby <- unlist(strsplit(x = filterby, split = "~"))
  if(is.null(groups)){
    if(verbose) cat("Evaluating globally\n")
    groups <- rep("SingleCell", ncol(mat)); names(groups) <- colnames(mat)
    if(is.na(filterby[3])) filterby[3] <- "SingleCell"
  }
  if(is.na(filterby[3])){
    filterby[3] <- name(sort(table(groups)))[1]
    warning("You may want to specify the group. Taking ", filterby[3], ".")
  }
  if(verbose) cat("Filtering by", filterby[3], "\n")
  mygenes <- unique(unlist(features));
  void <- stats_summary_table(
    mat = as.matrix(mat[rownames(mat) %in% mygenes, ]),
    groups = groups,
    moments = filterby[1],
    v = verbose
  ); if(verbose) print(summary(void))
  passed <- rownames(void)[void[, paste0(filterby[3], moments[filterby[1]])] > as.numeric(filterby[2])]
  if(verbose) str(features)
  features <- sapply(features, function(x) x[x %in% passed], simplify = FALSE)
  if(verbose) str(features)
  if(isTRUE(return_stats)) list(stats = void, list = features) else return(features)
}

# Previous colours
# c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','#EE0000','#b30000', '#670000')
# c("#fffeee", "#ffe080", "#ffc100", "#ff4d00", "#ff0000", "#EE0000", "#a10000", "#670000")
# c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000")
signature_scoring <- function(
  object,
  prefix = NULL,
  lsignatures = list(name = "orig.Cell.Cycle", S.Score = cc.genes$s.genes, G2M.Score = cc.genes$g2m.genes),
  group = FALSE, # group signatures with same name [and description] (NAME_description_source)
  confounders = "RNA.*res",
  reductions = list(pca = c('PC_1', 'PC_2'), tsne = c('tSNE_1', 'tSNE_2'), umap = c('UMAP_1', 'UMAP_2')),
  couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000"),
  violins_color = "mean",
  verbose = FALSE
){
  str_safe_remove <- function(x, word = "random123") gsub("_{1,}", "_", gsub(word, "", x))

  if(verbose) cat("-- Signature scoring --\n")
  prefix <- dircheck(prefix)
  if(verbose) cat("Output:", prefix, "\n")

  if(any(!confounders %in% colnames(object@meta.data))){
    confounders <- filters_columns(object@meta.data, include = confounders, maxn = 56, v = verbose)
  };
  confounders <- confounders[confounders %in% colnames(object@meta.data)]
  if(verbose) cat("Confounders", show_commas(confounders), "\n");
  Idents(object) <- 'orig.ident'

  reductions <- reductions[names(reductions) %in% names(object@reductions)]

  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  signature_list <- list()
  if(length(lsignatures)){
    for(i in 1:length(lsignatures)){
      mynameis <- sub("[[:alnum:]]{1,8}_", "", names(lsignatures[i]))
      if(verbose && isTRUE(group)) cat("@", mynameis, "\n")
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
      # print(sapply(class_list, length))
      # str(class_list)
      if(any(sapply(class_list, length) == 0)){
        warning("'", mynameis, "' has 0 feature value(s)"); next
      }
      names(class_list) <- paste0(casefold(str_safe_remove(names(class_list), "signature"), upper = TRUE), ".Score")
      signature_list <- c(signature_list, list(c(name = paste0("orig.", names(lsignatures[i])), class_list)))
    }
  }
  signame <- paste0(prefix, "signatures.csv")
  if(file.exists(signame)){
    if(verbose) cat("Pre-computed scores\n")
    signdf <- readfile(signame, stringsAsFactor = FALSE, check.names = FALSE, row.names = 1)
    object@meta.data <- joindf(object@meta.data, signdf, v = verbose)
  }
  scores_cols <- unname(unlist(sapply(signature_list, function(x) c(names(x)[-1], x[[1]]) )))
  for(scoring in signature_list){
    if(verbose) cat("\n@", sub("orig\\.", "", scoring$name), "\n")
    check_both <- c(scoring$name, names(scoring[-1])) # check both names
    if(any(!check_both %in% colnames(object@meta.data))){
      object <- ClassifyScoring(
        object = object, name = scoring$name, verbose = verbose,
        features = scoring[head(2:length(scoring), 2)]
      ); str(scoring); str(scores_cols)
      tvar <- FetchData(object, vars = scores_cols[scores_cols %in% colnames(object@meta.data)])
      if(exists("signdf")) tvar <- joindf(tvar, signdf)
      write.csv(tvar, file = signame)
    }

    if(verbose) cat("Plotting\n")
    ddfplot <- FetchData(object, vars = c(scoring$name, names(scoring[-1]), unname(unlist(reductions)), confounders))
    ddfplot$Signature <- as.character(ddfplot[, scoring$name])
    ddfplot$Signature <- ifelse(ddfplot$Signature == "None", "No", "Yes")
    # for(i in names(scoring[-1])) ddfplot[, i] <- scales::rescale(ddfplot[, i], to = c(0, 1))
    tvar <- table(ddfplot[, 'Signature'])
    subtitl <- paste0(paste0(names(tvar), ": ", unname(tvar)), collapse = "; ")

    if(verbose) cat(" # heatmap\n")
    graphics.off()
    fname <- paste0(prefix, gsub("orig\\.", "", casefold(scoring$name)), '_heatmap.pdf')
    if(!is.file.finished(fname)){
      pdf(fname, width = 10, height = 12, onefile = FALSE);
      p <- try(custom_heatmap(
        object = object, # should be normalised
        rnames = unname(unlist(scoring[-1])),
        orderby = scoring$name,
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
      if(verbose) cat(" #", redu, "\n")
      fname <- paste0(prefix, gsub("orig\\.", "", casefold(scoring$name)), '_', redu, '.pdf')
      if(is.file.finished(fname)) next
      p <- lapply(names(scoring[-1]), function(x){
        aesy <- aes_string(x = reductions[[redu]][1], y = reductions[[redu]][2], color = x)
        ggplot(data = ddfplot, mapping = aesy) +
        geom_point(size = 0.1) + scale_color_gradientn(colours = couls) +
        labs(colour = NULL, title = make_title(x))
      })
      tvar <- make_grid(length(p)) # t o determine when it's more
      pdf(fname, width = tvar[1]*tvar[3], height = tvar[2]*tvar[3]);
      print(cowplot::plot_grid(plotlist = p)); graphics.off()
    }
    for(confy in confounders){
      if(verbose) cat(" -", confy, "\n")
      fname <- paste0(
        prefix, gsub("orig\\.", "", casefold(scoring$name)),
        "_violin_", violins_color, "_", confy, '.pdf'
      )
      if(is.file.finished(fname)) next
      tvar <- length(table(ddfplot[, confy]))
      ddfplot[, confy] <- factormix(ddfplot[, confy])
      p <- violins(
        dat = ddfplot,
        xax = confy,
        yax = names(scoring[-1]),
        colour_by = violins_color
      ) + RotatedAxis()
      pdf(fname, width = 12, height = 8);
      print(p)
      graphics.off()
    }
  }

  if(verbose) cat("-- --------- ------- --\n")
  return(object)
}
