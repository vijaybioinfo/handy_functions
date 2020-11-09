#!/usr/bin/R

libs <- c('liger', 'fgsea', 'parallel', 'mdgsa', 'graphics', 'ggplot2')
libs <- sapply(libs, require, character.only = TRUE, quietly = TRUE)
if(any(!libs)) warning('Missing libraries ', paste0(names(libs[!libs]), collapse = " "))

#' Gene Set Enrichment Analysis.
#'
#' @description `gsea_tests` performs Gene Set Enrichment Analysis.
#'
#' @details This is a function that takes a named vector of fold changes or metrics
#' defined by `gsea_metric` or any other test and performs GSEA given a list
#' of signatures.
#'
#' @param res data.frame or numeric vector of values to use for ranking.
#' @param feature_names Column in which the feature names are.
#' @param metric Column in which the ranking metric is or just to indicate
#' what the metric is.
#' @param gsea_list Vector or list of feature names to test. It supports a path
#' to a file where each column is a gene list or a numeric vector for
#' org.Hs.GO2Symbol.list (liger).
#' @param method Method to use.
#' @param path Output path or prefix.
#' @param plot_it Plot every test.
#' @param verbose Show progress.
#' @keywords enrichment
#'
#' @return A fgsea-like object and creates Enrichment plots (if indicated).
#'
#' @author Ciro Ramírez-Suástegui, \email{ksuasteguic@gmail.com}
#' @references \url{https://en.wikipedia.org/wiki/Gene_set_enrichment_analysis}
#' @seealso \code{\link{gsea_metric}}
#'
#' @importFrom liger org.Hs.GO2Symbol.list bulk.gsea
#' @importFrom ggplot2 labs
#'
#' @examples
#' gsea_res = gsea_tests(res = group1_vs_group2)
#'
#' @export
#'

### Gene Set Enrichment Analysis ### -------------------------------------------
gsea_tests <- function(
  res,
  feature_names = "gene_name",
  metric = NULL,
  gsea_list = 5,
  method = c("fgsea", "liger"),
  path = "./",
  plot_it = TRUE,
  classical_plot = FALSE,
  verbose = FALSE,
  myseed = 27
){
  ncores <- min(c(ceiling(detectCores() / 5), 3))
  method <- match.arg(method)

  if(!grepl("/$", path) && dir.exists(path)) path <- paste0(path, "/")
  if(!grepl("_$", path) && !grepl("/$", path)) path <- paste0(path, "_")
  if(grepl("\\/$", path) && !dir.exists(path)) dir.create(path)
  if(verbose) cat("Output at:", path, "\n")

  if(feature_names %in% colnames(res)){
    res <- res[!is.na(res[, metric]), ]
    vals <- res[, metric]
    guniverse <- unique(sub("'", "", res[, feature_names]))  # get universe
  }else{
    res <- res[!is.na(res)]
    vals <- res; guniverse <- names(res)
  }
  names(vals) <- guniverse
  str(vals)

  # get a gene set
  if(file.exists(as.character(gsea_list[[1]][1]))){
    if(verbose) cat("Getting gene sets file\n")
    gset_list <- readfile(gsea_list, stringsAsFactors = F, header = TRUE, verbose = verbose)
  }else if(any(unlist(gsea_list) %in% guniverse)){
    if(verbose) cat("Given", class(gsea_list), "of genes\n")
    gset_list <- if(!is.list(gsea_list)) list(gene_list = gsea_list) else gsea_list
  }else if(is.numeric(gsea_list)){
    if(verbose) cat("Human Gene Ontology to HUGO Symbol list from 'liger'\n",
      "http://geneontology.org/page/download-go-annotations\n")
    gset_list <- org.Hs.GO2Symbol.list[seq(gsea_list)]
  }else{
    stop(paste0("Vector or list of feature names to test. It supports a",
                "path to a file where each column is a gene list or a",
                "numeric vector for org.Hs.GO2Symbol.list (liger)"))
  }; gset_list <- lapply(gset_list, function(x){ # removing NAs and empty elements
    y <- x[!is.na(x)]; gsub("'| ", "", y[y != ""])
  }); if(is.null(names(gset_list))) names(gset_list) <- paste0("g", 1:length(gset_list))
  if(verbose) cat("Number of sets:", length(gset_list), "\n")
  if(verbose) str(gset_list[head(1:length(gset_list))])
  if(all(sapply(gset_list, function(x) sum(x %in% names(vals)) ) == 0)) warning("No feature was found")

  ## Correcting for multiple-hypothesis
  if(method == "liger"){
    fgseaRes <- bulk.gsea(vals, gset_list)
    colnames(fgseaRes) <- c("pval", "padj", "ES", "edge")
    fgseaRes <- cbind(pathway = rownames(fgseaRes), fgseaRes)
    fgseaRes$size = sapply(gset_list, function(x) sum(x %in% guniverse) )
    fgseaRes <- data.table(fgseaRes)
    myheight <- 9
  }
  if(method == "fgsea"){
    fgseaRes <- fgsea(
      pathways = gset_list,
      stats = vals,
      minSize = 3,
      maxSize = 500,
      nproc = ncores,
      nperm = 10000
    )
    myheight <- 7
  }
  if(!plot_it) return(fgseaRes)

  topPathways <- c(
    fgseaRes[ES > 0, ][head(order(pval), n = 10), pathway],
    rev(fgseaRes[ES < 0, ][head(order(pval), n = 10), pathway])
  )
  suffix <- ifelse(is.null(metric), "", paste0("_", metric))
  if(length(topPathways) > 2){
    if(verbose) cat("Plotting top enrichements\n")
    fname <- paste0(path, "a1_", method, '_top_up_and_down_pathways', suffix, '.pdf')
    pdf(sub("\\/a1$", "_a1", fname), width = 8, ifelse(length(topPathways) >= 16, 16, 8))
    try(plotGseaTable(gset_list[topPathways], vals, fgseaRes, gseaParam = 1))
    graphics.off()
    tvar <- if("leadingEdge" %in% colnames(fgseaRes)) ncol(fgseaRes) - 1 else ncol(fgseaRes)
    write.csv(data.frame(fgseaRes[, 1:(tvar)]), file = sub("pdf", "csv", fname), row.names = FALSE)
  }

  go_names <- suppressWarnings(mdgsa::getGOnames(names(gset_list), verbose = TRUE));
  names(go_names) <- names(gset_list)
  go_names <- gsub(" activity|activity", "", sub("(^[A-z]{41,}).*", "\\1", go_names))
  if(verbose) cat("\n")
  for(indx in names(gset_list)){
    cat("----------------\n", indx, "\n")
    allgs <- sub("'", "", gset_list[[indx]]); allgs <- unique(allgs[allgs != ""])
    gs <- allgs[allgs %in% guniverse]
    if(verbose) cat(length(gs), "of", length(allgs), "found\n")
    if(length(gs) < 2) next
    suffy <- paste0(
      suffix, "_", method, "_",
      ifelse(length(gs) != length(allgs), paste0(length(gs), "of"), ""),
      length(allgs)
    )

    if(!is.na(go_names[indx])){
      if(verbose) cat(go_names[indx], "\n")
      indxname <- gsub("[[:punct:]]| {1,}", "_", paste0(indx, "_", go_names[indx]))
    }else{ indxname <- indx }
    indxname <- gsub("_{2,}", "_", indxname)

    fname <- paste0(path, indxname, suffy, '.pdf')
    pdf(fname, height = myheight, width = 7, onefile = TRUE)
    if(verbose) cat('Testing', metric, '\n'); set.seed(myseed)
    pval <- paste0("Adjusted P-value: ", round(fgseaRes[pathway == indx, "padj"], 5))
    if(method == "liger"){
      print(gsea(values = vals, geneset = gs, mc.cores = ncores, return.details = TRUE))
      if(!is.na(go_names[indx])) graphics::legend("bottomleft", legend = go_names[indx], bty = 'n', horiz = TRUE)
      graphics::legend("bottomright", legend = pval, bty = 'n', horiz = TRUE)
    }
    if(method == "fgsea"){
      p <- if(classical_plot == "fgsea"){
        plotEnrichment(pathway = gs, stats = vals)
      }else{
        gsea_plotEnrichment(
          pathway = gs,
          stats = vals,
          gsea.metric = metric_names[metric],
          gsea.enrichment.score = round(unlist(fgseaRes[pathway == indx, "ES"]), 2),
          gsea.normalized.enrichment.score = round(unlist(fgseaRes[pathway == indx, "NES"]), 2),
          padj = signif(unlist(fgseaRes[pathway == indx, "padj"]), 2),
          pval = signif(unlist(fgseaRes[pathway == indx, "pval"]), 2),
          classical_plot = classical_plot
        )
      }
      if(any(unlist(class(p)) == "ggplot")) print(p)
    }
    graphics.off()
  }
  return(fgseaRes)
}

metric_names = c('Signal-to-Noise ratio', 'T-Test', 'Ratio of Classes', 'Difference of Classes', 'log2 Ratio of Classes')
names(metric_names) = c('Signal2Noise', 'tTest', 'Ratio_of_Classes', 'Diff_of_Classes', 'log2_Ratio_of_Classes')

#' Metrics for Ranking Genes
#'
#' @description `gsea_metric` creates values useful for ranking in a GSEA.
#'
#' @details This is a function that takes a matrix of measurements and uses it to
#' calculate metric values from mean, sd, or a test that will be used by `gsea_tests`.
#'
#' @param mat Count matrix.
#' @param groups Column name containing the groups.
#' @param metric Type of metric.
#' @param rnames Feature names to use.
#' @param verbose Show progress.
#' @keywords mean sd
#'
#' @return A vector of values.
#'
#' @author Ciro Ramírez-Suástegui, \email{ksuasteguic@gmail.com}
#' @references \url{https://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html}
#' @seealso \code{\link{gsea_tests}}
#'
#' @importFrom stats stats_summary_table
#'
#' @examples
#' gsea_ready = gsea_metric(mat = edata, groups = c(rep("G1", 200), rep("G2", 200)))
#'
#' @export
#'

gsea_metric <- function(
  mat,
  groups,
  metric = c('Signal2Noise', 'tTest', 'Ratio_of_Classes', 'Diff_of_Classes', 'log2_Ratio_of_Classes'),
  verbose = FALSE
) {
  metric <- match.arg(metric)
  tvar <- unique(groups)
  if(verbose) cat(paste0(tvar, collapse = ", "), "\n")
  groups <- sub(paste0("^", tvar[1], "$"), "X1", groups)
  if(verbose) cat("X2 =", tvar[2], "\n")
  groups <- sub(paste0("^", tvar[2], "$"), "X2", groups)
  stattab <- stats_summary_table(
    mat = mat,
    groups = groups,
    moments = c('mn', 'sd'),
    v = verbose
  )
  if(verbose){ print(head(stattab)); print(tail(stattab)) }
  stattab2 <- data.frame(sapply(colnames(stattab), function(x){
    y <- stattab[, x]
    if(grepl("_sd", x)){ # check these lines... might explain the spikes when no filter is applied
      ttab <- as.matrix(data.frame(minsd = .2 * abs(stattab[, sub("_sd", "_mean", x)]), sd = y))
      y <- matrixStats::rowMaxs( ttab )
    }
    if(grepl("_mean", x)) y <- ifelse(y == 0, 1, y)
    y
  }), row.names = rownames(stattab)); stattab <- stattab2; rm(stattab2)
  dmean <- apply( stattab[, rev(grep("_mean", colnames(stattab)))], 1, diff )
  if(verbose) cat(metric, "\n")
  ymetric <- switch(metric,
    Signal2Noise = {
      if(verbose) cat("X1mean - X2mean\n")
      if(verbose) cat("---------------\n")
      if(verbose) cat("  X1sd + X2sd\n")
      dmean / apply( stattab[, grep("_sd", colnames(stattab))], 1, sum )
    }, tTest = {
      if(verbose) cat("          X1mean - X2mean\n")
      if(verbose) cat("  ------------------------------\n")
      if(verbose) cat("sqrt ((X1sd^2 / X1n) +( X2sd^2 / X2n))\n")
      ttab <- sweep(stattab[, grep("_sd", colnames(stattab))] ** 2, 2, table(groups), "/")
      dmean / sqrt( base::rowSums(ttab) )
    }, Ratio_of_Classes = {
      if(verbose) cat("X1mean / X2mean\n")
      stattab[, "X1_mean"] / stattab[, "X2_mean"]
    }, Diff_of_Classes = {
      if(verbose) cat("X1mean - X2mean\n"); dmean
    }, log2_Ratio_of_Classes = {
      if(verbose) cat("log2(X1mean / X2mean)\n")
      log2(stattab[, "X1_mean"] / stattab[, "X2_mean"])
    }
  )
  ymetric
}

#' GSEA from measurements
#'
#' @description `gsea_matrix` Performs the GSEA from a table.
#'
#' @details This is a function that takes a matrix of measurements like
#' expression or a table of metrics and performs GSEA on them.
#'
#' @param mat Count matrix or a table of metrics per comparison as columns.
#' @param groups Column name containing the groups. Accepts a vector where
#' the first element is the groups' column in the metadata.
#' @param metadata data.frame of the metadata.
#' @param metric Type of metric. It is adviced you specify it even if 'mat' is
#' the metrics themselves.
#' @param gsea_list List of vectors of the gene signatures.
#' @param method Method to use.
#' @param path Output path or prefix.
#' @param plot_it Plot every test.
#' @param pct_thr Percentage threshold for filtering before metric calculation.
#' @param verbose Show progress.
#' @keywords mean sd GSEA
#'
#' @return A vector of values.
#'
#' @author Ciro Ramírez-Suástegui, \email{ksuasteguic@gmail.com}
#' @references \url{https://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html}
#' @seealso \code{\link{gsea_tests}}
#'
#' @importFrom stats stats_summary_table
#'
#' @examples
#' gsea_res = gsea_matrix(mat = edata, groups = "Cell_type", metadata = mdata)
#'
#' @export
#'
gsea_matrix <- function(
  mat, # count matrix or a table of metrics per comparison as columns
  groups = NULL, # just the column and it'll do vs REST
  metadata = NULL,
  metric = 'custom',
  gsea_list = 5,
  method = c("fgsea", "liger"),
  path = "./",
  plot_it = TRUE,
  pct_thr = 1,
  classical_plot = FALSE,
  verbose = TRUE
){
  # Digesting
  method <- match.arg(method)
  if(metric == 'custom') warning("You might want to specify 'metric'")

  if(!grepl("/$", path) && dir.exists(path)) path <- paste0(path, "/")
  if(!grepl("_$", path) && !grepl("/$", path)) path <- paste0(path, "_")
  if(grepl("\\/$", path) && !dir.exists(path)) dir.create(path)
  if(verbose) cat("Output at:", path, "\n")

  rownames(mat) <- parse_ens_name(rownames(mat))

  if(any(colnames(mat) %in% rownames(metadata))){
    mat <- mat[, colnames(mat) %in% rownames(metadata)]
    metadata <- metadata[colnames(mat), ]
  }
  if(verbose) cat("Filtering matrix\n")
  # mat <- mat[Matrix::rowMeans(mat > 0) > pct_thr, ]
  mat <- mat[Matrix::rowMeans(mat) > pct_thr, ]
  features <- rownames(mat); group_column <- ""
  if(verbose) cat("Features:", length(features), head(features), tail(features), sep = "\n")
  if(!is.null(metadata)){
    if(groups[1] %in% colnames(metadata) && length(groups) == 1){
      if(verbose) cat("Given column\n")
      mygroups <- levels(factor(metadata[, groups]))
      mygroups <- data.frame(G1 = mygroups, G2 = rep("REST", length(mygroups)))
      mygroups$column <- groups
    }else if(!is.data.frame(groups)){ # these groups vs rest
      if(verbose) cat("Given column and a subset in it\n")
      mygroups <- if(!"REST" %in% groups){
        groups_comb <- gtools::combinations(length(groups[-1]), r = 2, v = groups[-1], set = FALSE)
        data.frame(G1 = groups_comb[, 1], G2 = groups_comb[, 2])
      }else{
        groups_comb <- groups[-1]; groups_comb <- groups_comb[!groups_comb %in% "REST"]
        data.frame(G1 = groups_comb[-1], G2 = rep("REST", length(groups_comb) - 1))
      }
      mygroups$column <- groups[1]
    }else{
      if(verbose) cat("Object class", class(mygroups), " given\n")
      mygroups <- groups
    }; group_column <- paste0(unique(mygroups$column), collapse = "-")
  }else{ # when the metrics are given
    mygroups <- data.frame(t(data.frame(
      strsplit(x = colnames(mat), split = "vs")
    )), row.names = NULL, stringsAsFactors = FALSE)
  }
  if(!'name' %in% colnames(mygroups)) mygroups$name <- do.call(paste, c(mygroups[, 1:2], sep = "vs"))
  mygroups$name <- gsub(pattern = "vsREST", replacement = "", mygroups$name)
  mygroups <- remove.factors(mygroups)
  if(verbose) cat("Tests:\n"); if(verbose) print(head(mygroups, 10))

  resnamef <- paste0(path, "metrics_per_", group_column, ".csv")
  res <- if(file.exists(resnamef)){
    if(verbose) cat("Previously found metrics table\n")
    read.csv(resnamef, row.names = 1, check.names = FALSE)
  }else if(file.exists(metric)){
    if(verbose) cat("Given metrics table file\n")
    read.csv(metric, row.names = 1, check.names = FALSE)
  }else if(!is.data.frame(mat)){
    if(verbose) cat("Calculating metrics\n")
    data.frame(row.names = features)
  }else{
    if(verbose) cat("Given metrics table data.frame\n")
    mat
  }
  res$gene_name <- rownames(res)
  if(!metric %in% colnames(res)){
    res$metric_t <- 0
    colnames(res) <- gsub("^metric_t$", metric, colnames(res))
  }

  # Performing tests
  void <- list()
  for(i in 1:nrow(mygroups)){
    if(verbose) cat("Comparison:", mygroups[i, ]$name, "\n")
    if(!mygroups[i, ]$name %in% colnames(res)){
      metadata$tmp <- as.character(metadata[, mygroups[i, ]$column])
      metadata_t <- if(as.character(mygroups[i, 2]) == "REST"){
        metadata$tmp <- ifelse(!metadata$tmp %in% mygroups[i, 1], "REST", metadata$tmp)
        metadata
      }else{
        metadata[metadata$tmp %in% unlist(mygroups[i, 1:2]), ]
      }
      tvar <- as.character(metadata_t$tmp)
      names(tvar) <- rownames(metadata_t)
      tvar <- c(tvar[tvar %in% mygroups[i, 1]], tvar[tvar %in% mygroups[i, 2]])

      res$comparison_name <- gsea_metric(
        groups = tvar,
        mat = mat[rownames(res), ],
        metric = metric,
        verbose = verbose
      )
      colnames(res) <- gsub("^comparison_name$", mygroups[i, ]$name, colnames(res))
      write.csv(res, file = resnamef)
    }
    if(mygroups[i, ]$name %in% colnames(res)){
      res[, metric] <- res[, mygroups[i, ]$name]
    }
    void[[mygroups[i, ]$name]] <- gsea_tests(
      res = res,
      feature_names = "gene_name",
      metric = metric,
      gsea_list = gsea_list,
      method = method,
      path = paste0(path, mygroups[i, ]$name, "/"),
      plot_it = plot_it,
      classical_plot = classical_plot,
      verbose = verbose
    )
  }
  return(void)
}

#' GSEA summaries
#'
#' @description `gsea_plot_summary` creates plots from GSEA results.
#'
#' @details This is a function that takes a list or an object from `gsea_tests`
#' and plots a heatmap and a radar plot. It also writes a table of the results.
#'
#' @param tests_list List of tests returned by `gsea_tests`.
#' @param type Type of values you want to plot. E. g., Normalized Enrichment Score.
#' @param path Output path or prefix.
#' @param padjthr Adjusted P-value threshold.
#' @param nesthr 'type' threshold.
#' @param na.rm Set NA to 0.
#' @param ... Arguments for pheatmap.
#' @keywords mean sd
#'
#' @return A vector of values.
#'
#' @author Ciro Ramírez-Suástegui, \email{ksuasteguic@gmail.com}
#' @references \url{https://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html}
#' @seealso \code{\link{gsea_tests}}
#'
#' @importFrom pheatmap pheatmap
#' @importFrom data.table rbindlist fwrite dcast.data.table
#' @importFrom rescale scales
#'
#' @examples
#' gsea_tab = gsea_plot_summary(tests_list = gsea_results)
#'
#' @export
#'

gsea_plot_summary <- function(
  tests_list,
  type = c("NES", "padj", "NLE", "ES", "LE"),
  path = "./",
  padjthr = 0.05,
  nesthr = 1,
  pathways = NULL,
  na.rm = FALSE,
  axes = list(columns = NULL, rows = NULL),
  radar_transpose = FALSE,
  ... # arguments for pheatmap
){
  type <- match.arg(type)
  round2 <- function(x) round(x, digits = 2)
  nlog10 <- function(x) -log10(x)
  plot_var <- list(
    NES = list(title = "Normalized Enrichment Score", transf = round2),
    padj = list(title = "Adjuste P-value", transf = nlog10),
    LE = list(title = "Normalized Leading Edge", transf = round2),
    ES = list(title = "Enrichment Score", transf = round2),
    LE = list(title = "Leading Edge", transf = round2)
  )
  plot_var <- plot_var[names(plot_var) %in% type]

  mygseas <- data.table::rbindlist(lapply(names(tests_list), function(x){
    y <- tests_list[[x]]
    z <- cbind(
      comparison = strsplit(basename(x), "_a1")[[1]][1],
      LE = sapply(y$leadingEdge, length),
      y
    ) # z$NLE = z$LE / z$size
    return(z)
  }))
  if(!is.null(pathways)) mygseas <- mygseas[as.character(mygseas$pathway) %in% pathways, ]
  thistab <- mygseas[mygseas$padj <= padjthr & (mygseas[[type]] >= nesthr | mygseas[[type]] <= -nesthr), ]

  fname0 <- paste0(path, "a1_gsea_summary_", nrow(thistab), "sets_", gsub(" .*", "", Sys.time()))
  data.table::fwrite(x = mygseas, file = paste0(fname0, ".txt"), sep = "\t")

  if(nrow(thistab) == 0){ warning("None passed p-value <= ", padjthr, " nor ", type, " of ", nesthr); return(0)}

  mysum <- data.table::dcast.data.table(thistab, comparison ~ pathway, value.var = names(plot_var))
  mysum <- plot_var[[1]]$transf(data.frame(mysum[, -1], stringsAsFactors = FALSE, row.names = mysum$comparison))
  tmysum <- mysum <- t(as.matrix(mysum))
  if(!is.null(axes$columns)) mysum <- mysum[, axes$columns[axes$columns %in% colnames(mysum)]]
  if(!is.null(axes$rows)) mysum <- mysum[axes$rows[axes$rows %in% rownames(mysum)], ]
  cc = if(is.null(axes$columns)) TRUE else FALSE
  cr = if(is.null(axes$rows)) TRUE else FALSE
  if(isTRUE(na.rm)) mysum[is.na(mysum)] <- 0
  # annoc <- data.frame(
  #   Cluster = names(identnames), Identity = unname(identnames), row.names = names(identnames)
  # )
  # annoc <- annoc[, sapply(annoc, function(x) !any(is.na(x)) ), drop = FALSE]
  # anncolist <- lapply(annoc, function(x) v2cols(x, grcols, v = TRUE) )

  x <- pheatmap::pheatmap(
    mat               = mysum,
    cluster_rows      = cr,
    cluster_cols      = cc,
    scale             = 'none',
    border_color      = NA,
    show_colnames     = TRUE,
    show_rownames     = TRUE,
    main              = plot_var[[1]]$title,
    # annotation_col    = annoc,
    # annotation_colors = anncolist,
    na_col            = "#BEBEBE",
    annotation_legend = TRUE,
    annotation_names_col = TRUE,
    annotation_names_row = TRUE,
    drop_levels       = TRUE,
    filename          = paste0(fname0, "_", names(plot_var), ".pdf"),
    width = 10, height = 7,
    ...
  );

  # # in case it failed when not setting NAs to 0
  # mysum <- tmysum[x$tree_row$labels[x$tree_row$order], x$tree_col$labels[x$tree_col$order]]
  # cr <- cc <- FALSE

  # Radar plot
  library(ggradar) # devtools::install_github("ricardo-bion/ggradar")
  library(dplyr)
  if(isTRUE(radar_transpose)) mysum <- data.frame(t(mysum))
  mysum[is.na(mysum)] <- 0
  p <- mysum %>%
    tibble::as_tibble(rownames = "group") %>%
    mutate_at(vars(-group), scales::rescale) %>% ggradar()
  pdf(paste0(fname0, "_", names(plot_var), "_radar.pdf"), width = 15, height = 15)
  print(p)
  dev.off()

  return(mysum)
}

remove.factors <- function (df) {
  for (varnum in 1:length(df)) {
    if ("factor" %in% class(df[, varnum])) {
      df[varnum] = as.character(df[, varnum])
    }
  }
  return(df)
}

parse_ens_name <- function(x, keepens = FALSE){
  tvar <- "X1234" # removing the ensembl name
  if(grepl("ENS", x[1]) && grepl("[0-9]", x[1])){
    ensloc <- all(grepl("^ENS", x[1:10])) # find which side the names are
    if(isTRUE(keepens)) ensloc <- !ensloc
    tvar <- ifelse(ensloc, ".*_", "_.*")
    tvar <- paste0(c(tvar, ifelse(ensloc, ".*\\|", "\\|.*")), collapse = "|")
  }
  tvar <- gsub("\\-", "_", gsub(tvar, "", x)) # so we can keep the dashes!
  tvar <- make.names(tvar, unique = TRUE, allow_ = TRUE)
  gsub("_", "-", tvar) # returning dashes ;)
}

gsea_plotEnrichment <- function(
  pathway,
  stats,
  gseaParam = 1,
  ticksSize = 0.2,
  name = "Gene Set Enrichment Analysis",
  gsea.metric = NA,
  gsea.enrichment.score = NULL,
  gsea.normalized.enrichment.score = NULL,
  pval = NULL,
  padj = NULL,
  classical_plot = FALSE
){
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
  statsAdj <- statsAdj / max(abs(statsAdj))
  pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway, returnAllExtremes = TRUE)
  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8
  x = y = NULL
  ## Create GSEA plot
  if(!isTRUE(classical_plot)){
    g <- ggplot(data = toPlot, mapping = aes(x = x, y = y)) +
      geom_point(color = "darkblue", size = 0.1) +
      geom_hline(yintercept = max(tops), colour = "#ff7070", linetype = "dashed") +
      geom_hline(yintercept = min(bottoms), colour = "#ff7070", linetype = "dashed") +
      geom_hline(yintercept = 0, colour = "black", linetype = "dashed") +
      geom_line(color = "blue") +
      geom_segment(
        data = data.frame(x = pathway),
        mapping = aes(x = x, y = -diff/4, xend = x, yend = diff/4),
        size = ticksSize
      ) +
      theme_bw() +
      theme(
        panel.border = element_blank(), panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(), axis.text.y = element_text(colour = "#BEBEBE")
      ) +
      labs(x = "Rank", y = "Enrichment Score", caption = padj)
    return(g)
  }else{
    # Classic GSEA plot from https://github.com/PeeperLab/Rtoolbox/blob/master/R/ReplotGSEA.R
    gsea.hit.indices = toPlot$x
    gsea.rnk_metric = stats
    gsea.es.profile = toPlot$y
    gsea.gene.set = name
    enrichment.score.range = NULL
    metric.range = NULL
    class.name = ""
    gsea.fdr = padj
    gsea.p.value = pval

    if (is.null(metric.range)) {
      metric.range <- c(min(gsea.rnk_metric), max(gsea.rnk_metric))
    }

    # Set enrichment score range
    if (is.null(enrichment.score.range)) {
      enrichment.score.range <- range(gsea.es.profile[gsea.es.profile!=0])
    }

    # Save default for resetting
    def.par <- par(no.readonly = TRUE)

    # # Create a new device of appropriate size
    # dev.new(width = 3, height = 3)

    # Create a division of the device
    gsea.layout <- layout(matrix(c(1, 2, 3, 4)), heights = c(1.7, 0.5, 0.2, 2))
    # layout.show(gsea.layout)

    # Create plots
    par(mar = c(0, 5, 2, 2))
    plot(
      c(1, gsea.hit.indices, length(gsea.rnk_metric)),
      c(0, gsea.es.profile, 0), type = "l", col = "red", lwd = 1.5, xaxt = "n",
      xaxs = "i", xlab = "", ylab = "Enrichment score (ES)",
      ylim = enrichment.score.range,
      main = list(gsea.gene.set, font = 1, cex = 1),
      panel.first = {
        abline(h = seq(round(enrichment.score.range[1], digits = 1), enrichment.score.range[2], 0.1), col = "gray95", lty = 2)
        abline(h = 0, col = "gray50", lty = 2)
      }
    )
    plot.coordinates <- par("usr")
    text_description <-  paste(
      "Nominal p-value:", gsea.p.value, "\nFDR:", gsea.fdr, "\nES:",
      gsea.enrichment.score, "\nNormalized ES:", gsea.normalized.enrichment.score, "\n"
    )
    if(gsea.enrichment.score < 0) {
      text(length(gsea.rnk_metric) * 0.01, plot.coordinates[3] * 0.98, text_description, adj = c(0, 0))
    } else {
      text(length(gsea.rnk_metric) * 0.99, plot.coordinates[4] - ((diff(plot.coordinates[4:3])) * 0.03),
       text_description, adj = c(1, 1))
    }

    par(mar = c(0, 5, 0, 2))
    plot(0, type = "n", xaxt = "n", xaxs = "i", xlab = "", yaxt = "n",
         ylab = "", xlim = c(1, length(gsea.rnk_metric)))
    abline(v = gsea.hit.indices, lwd = 0.75)

    par(mar = c(0, 5, 0, 2))
    rank.colors <- gsea.rnk_metric - metric.range[1]
    rank.colors <- rank.colors / (abs(diff(metric.range[2:1])))
    rank.colors <- ceiling(rank.colors * 255 + 1)
    tryCatch({
      rank.colors <- colorRampPalette(c("blue", "white", "red"))(256)[rank.colors]
    }, error = function(e) {
      warning("Please use the metric.range argument to provide a metric range that",
           "includes all metric values")
    })
    # Use rle to prevent too many objects
    rank.colors <- rle(rank.colors)
    barplot(
      matrix(rank.colors$lengths), col = rank.colors$values, border = NA, horiz = TRUE, xaxt = "n",
      xlim = c(1, length(gsea.rnk_metric))
    )
    box()
    text(length(gsea.rnk_metric) / 2, 0.7, labels = class.name)
    text(length(gsea.rnk_metric) * 0.01, 0.7, "Positive", adj = c(0, NA))
    text(length(gsea.rnk_metric) * 0.99, 0.7, "Negative", adj = c(1, NA))

    par(mar = c(5, 5, 0, 2))
    rank.metric <- rle(round(gsea.rnk_metric, digits = 2))
    is_there_a_metric = isTRUE((length(gsea.metric) > 0) & !is.na(gsea.metric))
    plot(
      gsea.rnk_metric, type = "n", xaxs = "i",
      xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk_metric)),
      ylim = metric.range, yaxs = "i",
      ylab = if(!is_there_a_metric){ "Ranking metric" }else{ gsea.metric },
      panel.first = abline(
        h = seq(metric.range[1] / 2, diff(metric.range[2:1]) / 4, metric.range[2] / 2),
        col = "gray95", lty = 2)
    )
    barplot(
      unname(rank.metric$values), col = "lightgrey", lwd = 0.1, xaxs = "i",
      xlab = "Rank in ordered gene list", xlim = c(0, length(gsea.rnk_metric)),
      ylim = c(-1, 1), yaxs = "i", width = rank.metric$lengths, border = NA,# las = 2,
      ylab = ifelse(!is_there_a_metric, "Ranking metric", gsea.metric), space = 0, add = TRUE
    )
    box()

    # Reset to default
    par(def.par)
  }
}
