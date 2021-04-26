#!/usr/bin/R

#' From table sumarisation
#'
#' This function is to summarise aspects of a matrix.
#' @param mat matrix or data.frame of numeric values
#' @param groups named (column names from mat) vector
#' @param moments stats to calculate, mean (mn), median 'md', percentage (p)
#' mean or median of percentage (pmn or pmd) and quantiles (q)
#' @param expr_cutoff threshold for a value to be considered 'expressed'
#' @param cnames columns to use
#' @param rnames rows to use
#' @param do_log add shited log 2 of the values, default is FALSE
#' @keywords stat_summary
#'
#' @return Returns a data.frame with statistics of your data
#'
#' @importFrom handy cbindList
#' @importFrom gtools mixedsort
#'
#' @export
#'
#' @examples
#' addinfo <- stats_summary_table(mat = edata, groups = annot_group)
#'

moments = c(
  mn = "_mean", md = "_median", p = "_percentage",
  pmn = "_meanOfPositive", pmd = "_medianOfPositive",
  sd = "_sd", q = "_q"
)
stats_summary_table <- get_stat_report <- function(
    mat,
    groups = NULL, # named [with mat colnames] vector
    cnames = NULL, # column names
    rnames = NULL, # row names
    moments = c('mn', 'p', 'pmn'),
    expr_cutoff = 0,
    do_log = FALSE, # add log2 data
    datatype = "",
    sep_str = "_",
    verbose = FALSE
  ){
  if(verbose) cat('------- Stat summary -------\n')

  if(ncol(mat) == 0 || nrow(mat) == 0){
    stop("No data given")
  }

  if(is.null(rownames(mat))){
    warning('\'matrix\' should have row names')
    rownames(mat) <- 1:nrow(mat)
  }
  if(is.null(colnames(mat))){
    warning('\'matrix\' needs to have column names')
    colnames(mat) <- 1:ncol(mat)
  }

  if(verbose) cat('Type of data:', datatype, '\n')
  if(!is.null(groups)){
    if(is.unsorted(groups)){
      if(verbose) cat('Not ordered')
      if(is.null(cnames)){
        if(verbose) cat(' - ordering')
        groups <- gtools::mixedsort(groups)
        if(is.null(names(groups))) stop('\'groups\' needs to be a named vector')
      }
      if(verbose) cat('\n')
    }
  }else{
    if(verbose) cat('No groups given - treating all as one\n')
    groups <- rep("Identity", ncol(mat))
    names(groups) <- colnames(mat)
  }

  if(is.null(cnames)) cnames <- names(groups) else groups <- groups[names(groups) %in% cnames]
  if(is.null(rnames)) rnames <- rownames(mat)
  if(verbose) cat('Cells:', length(cnames), head(cnames), '\n')
  if(verbose) cat('Genes:', length(rnames), head(rnames), '\n')
  mat <- try(as.matrix(mat[rnames, cnames, drop = FALSE]))
  if(class(mat) == 'try-error'){
    cat("Missing genes:", str(rnames[!rnames %in% rownames(mat)]), "\n")
    cat("Missing cells:", str(cnames[!cnames %in% colnames(mat)]), "\n")
    str(mat)
  }

  groups <- factor(groups, levels = unique(groups))
  grps <- levels(groups)
  if(length(grps) == 0) grps <- unique(groups)
  if(verbose) cat('Groups (n=',length(grps),'): ', str(grps), '\n', sep = "")

  if(verbose) cat('Calculating groups stats\n') ######################################
  sum_stat <- data.frame(row.names = rnames, stringsAsFactors = FALSE, check.names = FALSE)

  if(requireNamespace("Rfast", quietly = TRUE)){
    if(verbose) cat(' * Using Rfast\n')
    medians_fun = Rfast::rowMedians; means_fun = Rfast::rowmeans
  }else{
    medians_fun = matrixStats::rowMedians; means_fun = matrixStats::rowMeans2
  }

  if('b' %in% moments && verbose) cat(' - Bases\n')

  if('md' %in% moments){
    if(verbose) cat(' - Median\n')
    if('b' %in% moments) sum_stat[, paste0('Bmedian', datatype)] <- medians_fun(mat, na.rm = TRUE)
    medians <- apply(mat, 1, function(vec) tapply(vec, groups, median, na.rm = TRUE) )
    if(length(grps) != 1) medians <- t(medians)
    medians <- data.frame(medians, check.names = FALSE);
    colnames(medians) <- paste0(colnames(medians), sep_str, 'median', datatype)
    sum_stat <- cbind(sum_stat, medians)
  }

  if('mn' %in% moments){
    if(verbose) cat(' - Mean\n')
    if('b' %in% moments) sum_stat[, paste0('Bmean', datatype)] <- means_fun(mat, na.rm = TRUE)
    means <- apply(mat, 1, function(vec) tapply(vec, groups, matrixStats::mean2, na.rm = TRUE) )
    if(length(grps) != 1) means <- t(means)
    means <- data.frame(means, check.names = FALSE); colnames(means) <- grps
    colnames(means) <- paste0(colnames(means), sep_str, 'mean', datatype)
    sum_stat <- cbind(sum_stat, means)
  }

  if('sd' %in% moments){
    if(verbose) cat(' - Standard deviation\n')
    if('b' %in% moments) sum_stat[, paste0('Bsd', datatype)] <- matrixStats::rowSds(mat, na.rm = TRUE)
    means <- apply(mat, 1, function(vec) tapply(vec, groups, stats::sd, na.rm = TRUE) )
    if(length(grps) != 1) means <- t(means)
    means <- data.frame(means, check.names = FALSE); colnames(means) <- grps
    colnames(means) <- paste0(colnames(means), sep_str, 'sd', datatype)
    sum_stat <- cbind(sum_stat, means)
  }

  if('p' %in% moments){
    if(verbose) cat(' - Percentages\n')
    if('b' %in% moments) sum_stat[, paste0('BexprFrac', datatype)] <- rowSums(mat > expr_cutoff, na.rm = TRUE) / ncol(mat) * 100
    poscells <- apply(mat, 1, function(vec) tapply(vec, groups, function(cc) sum(cc > expr_cutoff, na.rm = TRUE) / length(cc) ) )
    if(length(grps) != 1) poscells <- t(poscells)
    poscells <- data.frame(poscells * 100, check.names = FALSE);
    colnames(poscells) <- grps
    colnames(poscells) <- paste0(colnames(poscells), sep_str, 'percentage', datatype)
    sum_stat <- cbind(sum_stat, poscells)
  }

  if('pmd' %in% moments){
    if(verbose) cat(' - Median of positive values\n')
    poscent <- apply(mat, 1, function(vec) tapply(vec, groups, function(cc) median(cc[cc > expr_cutoff], na.rm = TRUE) ) )
    if(length(grps) != 1) poscent <- t(poscent)
    poscent <- data.frame(poscent, check.names = FALSE);
    colnames(poscent) <- paste0(colnames(poscent), sep_str, 'medianOfPositive', datatype)
    sum_stat <- cbind(sum_stat, poscent)
  }

  if('pmn' %in% moments){
    if(verbose) cat(' - Mean of positive values\n')
    if('b' %in% moments){
      poscells <- apply(mat, 1, function(vec) tapply(vec, rep("ALL", ncol(mat)), function(cc) mean(cc[cc > 0], na.rm = TRUE) ) )
      sum_stat[, paste0('BmeanOfPositive', datatype)] <- poscells
    }
    poscent <- apply(mat, 1,function(vec) tapply(vec, groups, function(cc) mean(cc[cc > expr_cutoff], na.rm = TRUE) ) )
    if(length(grps) != 1) poscent <- t(poscent)
    poscent <- data.frame(poscent, check.names = FALSE);
    colnames(poscent) <- paste0(colnames(poscent), sep_str, 'meanOfPositive', datatype)
    sum_stat <- cbind(sum_stat, poscent)
  }

  if('q' %in% moments){
    if(verbose) cat(' - Quantiles\n')
    quants <- dplyr::bind_cols(lapply(c(25, 50, 75), function(myq){
      quant <- apply(mat, 1, function(vec) tapply(vec, groups, quantile, probs = myq / 100, na.rm = TRUE) )
      if(length(grps) != 1) quant <- t(quant)
      quant <- data.frame(quant, check.names = FALSE)
      colnames(quant) <- paste0(colnames(quant), sep_str, 'q', myq, datatype)
      return(quant)
    }))
    quants <- quants[, order(colnames(quants))]
    sum_stat <- cbind(sum_stat, quants)
  }

  if(verbose) cat('------- ---- ------- -------\n')
  return(sum_stat)
}
