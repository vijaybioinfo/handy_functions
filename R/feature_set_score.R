#!/usr/bin/R

###################################
# Clustering analysis: signatures #
###################################

# This script generates signature scores

# "R/stats_summary_table.R", "devel/filters.R", "devel/plots.R", "clustering_utilities.R"


#' @title Clean feature lists
#' @description Cleans a list of features based on the percentage or other
#' summary statistics present in `stats_summary_table`.
#' @param mat Matrix.
#' @param features List of features.
#' @param groups Named (columns in 'mat') vector of groups, Default: NULL.
#' @param filterby Instructions on how to filter 'moment~value', Default: 'p~0'.
#' You have the option of adding ~group at the end to filter based on 'groups'.
#' @param return_stats Whether or not to return summary stats, Default: FALSE.
#' @param verbose Verbosity, Default: FALSE.
#' @param ... Extra parameters for `stats_summary_table`.
#' @return Filtered features list or a list of summary statistics and list of
#' features.
#' @examples
#' \dontrun{
#'  if(interactive()) features_l <- clean_feature_list(edata, features_list)
#' }
#' @rdname clean_feature_list
#' @export

clean_feature_list <- function(
  mat,
  features,
  groups = NULL,
  filterby = 'p~0',
  return_stats = FALSE,
  verbose = FALSE,
  ...
){
  filterby <- unlist(strsplit(x = filterby, split = "~"))
  if(is.null(groups)){
    groups <- setNames(rep("all", ncol(mat)), colnames(mat))
    if(is.na(filterby[3])) filterby[3] <- "all"
  }
  if(is.na(filterby[3])){
    filterby[3] <- names(sort(table(groups)))[1]
    warning("You may want to specify the group. Taking ", filterby[3], ".")
  }
  str(groups)
  if(verbose) cat("Filtering by", filterby[3], "\n")
  if(!exists("stats_summary_table"))
    source(paste0("https://raw.githubusercontent.com/vijaybioinfo/",
      "handy_functions/master/R/stats_summary_table.R"))
  void <- stats_summary_table(
    mat = as.matrix(mat[unique(unlist(features)), ]),
    groups = groups,
    moments = filterby[1],
    v = verbose
  ); if(verbose) print(summary(void))
  tmp <- void[, paste0(filterby[3], moments[filterby[1]])]
  passed <- rownames(void)[tmp > as.numeric(filterby[2])]
  if(verbose) str(features)
  features <- sapply(features, function(x) x[x %in% passed], simplify = FALSE)
  if(verbose) str(features)
  if(isTRUE(return_stats)) list(stats = void, list = features) else return(features)
}

signature_list_process <- function(
  signature_list, grouping = FALSE, verbose = FALSE
){
  list_processed <- list()
  if(verbose) cat("Getting lists ready\n")
  for(i in 1:length(signature_list)){
    mynameis <- names(signature_list[i])
    tmp <- c(names(list_processed), sapply(list_processed, names))
    if(mynameis %in% tmp) next
    if(verbose) cat("@", mynameis, "\n")
    mylists <- signature_list[i]
    if(isTRUE(grouping)){
      # First part of the name can be used to aggregate groups of list
      # e. g., tcell_treg, tcell_tfh, tcell_tfr
      tvar <- sub("^([[:alnum:]]{1,})_.*", "\\1", names(signature_list[i]))
      tmp <- grep(tvar, names(signature_list))
      if(length(tmp) > 1){
        if(verbose) cat(" - grouping in ", tvar, ": ", sep = "")
        mylists <- signature_list[tmp]; mynameis <- tvar
        if(verbose) cat(names(mylists), "\n")
      }
    }#; names(mylists) <- paste0(names(mylists), ".score")
    list_processed[[mynameis]] <- c(list(name = mynameis), mylists)
  }; return(x = list_processed)
}

signature_scoring <- function(
  object,
  metadata = NULL,
  prefix = NULL,
  signature_list = list(name = "Cell.Cycle", S = Seurat::cc.genes$s.genes, G2M = Seurat::cc.genes$g2m.genes),
  confounders = "RNA.*res",
  grouping = FALSE, # grouping lists into one NAME_description_source
  plotting = TRUE,
  reductions = list(pca = c('PC_1', 'PC_2'), tsne = c('tSNE_1', 'tSNE_2'), umap = c('UMAP_1', 'UMAP_2')),
  couls = c("#fffffa", "#fffeee", "#ffe080", "#ffc100", "#ff0000", "#EE0000", "#a10000", "#670000"),
  violins_color = "mean",
  verbose = FALSE,
  ...
){
  if(verbose) cat("-- Signature scoring --\n")
  if(verbose) cat("Output:", prefix, "\n")
  if(!is.null(metadata)){
    object <- object[, rownames(metadata)]
    object@meta.data <- joindf(object@meta.data, metadata)
  }

  if(any(!confounders %in% colnames(object@meta.data))){
    confounders <- filters_columns(object@meta.data, confounders, verbose = verbose)
  };
  confounders <- confounders[confounders %in% colnames(object@meta.data)]
  if(verbose)
    cat("Confounders:", stringr::str_wrap(paste0(confounders, collapse = ", ")), "\n");
  Idents(object) <- 'orig.ident'
  signature_list_processed <- signature_list_process(signature_list, grouping = grouping)
  signame <- paste0(prefix, "signatures.csv")
  if(file.exists(signame)){
    if(verbose) cat("Pre-computed scores\n")
    signdf <- readfile(signame, stringsAsFactor = FALSE, check.names = FALSE, row.names = 1)
    object@meta.data <- joindf(object@meta.data, signdf)
  }
  scores_cols <- unique(unlist(lapply(signature_list_processed,
    function(x) c(x[[1]], names(x)[-1]) )))
  for(scoring in signature_list_processed){
    if(verbose) cat("\n@", scoring$name, "\n")
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

    if(isFALSE(plotting)) next
    if(verbose) cat("Plotting\n")
    ddfplot <- Seurat::FetchData(
      object, vars = c(scoring$name, names(scoring[-1]), unname(unlist(reductions)), confounders)
    )
    ddfplot$Signature <- as.character(ddfplot[, scoring$name])
    ddfplot$Signature <- ifelse(ddfplot$Signature == "None", "No", "Yes")
    # for(i in names(scoring[-1])) ddfplot[, i] <- scales::rescale(ddfplot[, i], to = c(0, 1))
    tvar <- table(ddfplot[, 'Signature'])
    subtitl <- paste0(paste0(names(tvar), ": ", unname(tvar)), collapse = "; ")

    if(verbose) cat(" # heatmap\n")
    fname <- paste0(prefix, gsub("orig\\.", "", casefold(scoring$name)), '_heatmap.pdf')
    if(!(file.exists(fname) && file.size(fname) > 3620)){
      pdf(fname, width = 10, height = 12, onefile = FALSE);
      p <- try(custom_heatmap(
        object = object,
        rnames = unname(unlist(scoring[-1])),
        orderby = scoring$name,
        sample_it = c("orig.ident", '3000'),
        categorical_col = c(confounders, "Signature"),
        feature_order = "pca",
        verbose = FALSE,
        show_colnames = FALSE
      ))
      graphics.off()
    }
    for(i in names(reductions)){
      if(verbose) cat(" #", i, "\n")
      fname <- paste0(prefix, gsub("orig\\.", "", casefold(scoring$name)), '_', i, '.pdf')
      if(file.exists(fname) && file.size(fname) > 3620) next
      p <- lapply(names(scoring[-1]), function(x){
        aesy <- aes_string(x = reductions[[i]][1], y = reductions[[i]][2], color = x)
        ggplot(data = ddfplot, mapping = aesy) +
        geom_point(size = 0.1) + scale_color_gradientn(colours = couls) +
        labs(colour = NULL, title = make_title(x))
      })
      tvar <- plot_size(make_grid(length(p))) # t o determine when it's more
      pdf(fname, width = tvar[2], height = tvar[1]);
      print(cowplot::plot_grid(plotlist = p)); graphics.off()
    }
    for(confy in confounders){
      if(verbose) cat(" -", confy, "\n")
      fname <- paste0(
        prefix, gsub("orig\\.", "", casefold(scoring$name)),
        "_violin_", violins_color, "_", confy, '.pdf'
      )
      if(file.exists(fname) && file.size(fname) > 3620) next
      p <- violins(
        dat = ddfplot,
        xax = confy,
        yax = names(scoring[-1]),
        colour_by = violins_color
      ) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
      pdf(fname, width = 12, height = 8);
      print(p)
      graphics.off()
    }
  }

  if(verbose) cat("-- --------- ------- --\n")
  return(object)
}
