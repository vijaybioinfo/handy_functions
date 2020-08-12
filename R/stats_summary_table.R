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
#' @param indiv_values add ordered matrix at the end, default is FALSE
#' @keywords stat_summary
#'
#' @return Returns a data.frame with statistics of your data
#'
#' @importFrom handy cbindList
#'
#' @export
#'
#' @examples
#' addinfo <- get_stat_report(mat = edata, groups = annot_group)
#'

get_stat_report <- function(
    mat,
    groups = NULL, # named [with mat colnames] vector
    cnames = NULL, # column names
    rnames = NULL, # row names
    moments = c('mn', 'p', 'pmn'),
    expr_cutoff = 0,
    do_log = FALSE, # add log2 data
    datatype = "",
    indiv_values = FALSE, # add matrix at the end
    with_rownams = FALSE,
    v = FALSE
  ){
  # moments = c('mn', 'md', 'p', 'pmn', 'pmd', 'q'),
  if(v) cat('------- Stat summary -------\n')

  if(is.null(rownames(mat))){
    warning('\'matrix\' should have row names')
    rownames(mat) <- 1:nrow(mat)
  }
  if(is.null(colnames(mat))){
    warning('\'matrix\' needs to have column names')
    colnames(mat) <- 1:ncol(mat)
  }

  if(v) cat('Type of data:', datatype, '\n')
  if(!is.null(groups)){
    if(is.unsorted(groups)){
      if(v) cat('Not ordered')
      if(is.null(cnames)){
        if(v) cat(' - ordering')
        groups <- gtools::mixedsort(groups)
        if(is.null(names(groups))) stop('\'groups\' needs to be a named vector')
      }
      if(v) cat('\n')
    }
  }else{
    if(v) cat('No groups given - treating all as one\n')
    groups <- rep("Identity", ncol(mat))
    names(groups) <- colnames(mat)
  }

  if(is.null(cnames)) cnames <- names(groups) else groups <- groups[names(groups) %in% cnames]
  if(is.null(rnames)) rnames <- rownames(mat)
  if(v) cat('Cells:', length(cnames), '\nGenes: ', length(rnames),'\n')
  mat <- try(as.matrix(mat[rnames, cnames, drop = FALSE]))
  if(class(mat) == 'try-error'){
    cat("Missing genes:", commas(rnames[!rnames %in% rownames(mat)]), "\n")
    cat("Missing cells:", commas(cnames[!cnames %in% colnames(mat)]), "\n")
    str(mat)
  }

  groups <- factor(groups, levels = unique(groups))
  grps <- levels(groups)
  if(length(grps) == 0) grps <- unique(groups)
  if(v) cat('Groups (n=',length(grps),'): ', commas(grps, hn = 10), '\n', sep = "")

  if(v) cat('Calculating groups stats\n') ######################################
  sum_stat <- data.frame(rownams = rnames, stringsAsFactors = F, check.names = F)
  # res <- microbenchmark(
  #     rowSums(mat > 0, na.rm = T) / ncol(mat) * 100, # looks like this is faster
  #     rowMeans(mat > 0, na.rm = T) * 100
  # )

  if('b' %in% moments && v) cat(' - Bases\n')

  if('md' %in% moments){
    if(v) cat(' - Median\n')
    # if('b' %in% moments) sum_stat[, paste0('Bmedian', datatype)] <- Biobase::rowMedians(mat, na.rm = T)
    medians <- apply(mat, 1, function(vec) tapply(vec, groups, median, na.rm = T) )
    if(length(grps) != 1) medians <- t(medians)
    medians <- data.frame(medians, check.names = F); # colnames(medians) <- grps
    colnames(medians) <- paste0(colnames(medians), '_median', datatype)
    sum_stat <- cbind(sum_stat, medians)
  }

  if('mn' %in% moments){
    if(v) cat(' - Means\n')
    if('b' %in% moments) sum_stat[, paste0('Bmean', datatype)] <- Matrix::rowMeans(mat, na.rm = T)
    means <- apply(mat, 1, function(vec) tapply(vec, groups, Matrix::mean, na.rm = T) )
    if(length(grps) != 1) means <- t(means)
    means <- data.frame(means, check.names = F); colnames(means) <- grps
    colnames(means) <- paste0(colnames(means), '_mean', datatype)
    sum_stat <- cbind(sum_stat, means)
  }

  if('sd' %in% moments){
    if(v) cat(' - SD\n')
    if('b' %in% moments) sum_stat[, paste0('Bmean', datatype)] <- matrixStats::rowSds(mat, na.rm = T)
    means <- apply(mat, 1, function(vec) tapply(vec, groups, sd, na.rm = T) )
    if(length(grps) != 1) means <- t(means)
    means <- data.frame(means, check.names = F); colnames(means) <- grps
    colnames(means) <- paste0(colnames(means), '_sd', datatype)
    sum_stat <- cbind(sum_stat, means)
  }

  if('p' %in% moments){
    if(v) cat(' - Percentages\n')
    if('b' %in% moments) sum_stat[, paste0('BexprFrac', datatype)] <- rowSums(mat > expr_cutoff, na.rm = T) / ncol(mat) * 100
    poscells <- apply(mat, 1, function(vec) tapply(vec, groups, function(cc) sum(cc > expr_cutoff, na.rm = T) / length(cc) ) )
    if(length(grps) != 1) poscells <- t(poscells)
    poscells <- data.frame(poscells * 100, check.names = F); # colnames(poscells) <- grps
    # poscells <- data.frame(poscells/ncol(mat)*100, check.names = F);
    colnames(poscells) <- grps
    colnames(poscells) <- paste0(colnames(poscells), '_exprFrac', datatype)
    sum_stat <- cbind(sum_stat, poscells)
  }

  if('pmd' %in% moments){
    if(v) cat(' - Median of positive values\n')
    poscent <- apply(mat, 1, function(vec) tapply(vec, groups, function(cc) median(cc[cc > expr_cutoff], na.rm = T) ) )
    if(length(grps) != 1) poscent <- t(poscent)
    poscent <- data.frame(poscent, check.names = F); # colnames(poscent) <- grps
    colnames(poscent) <- paste0(colnames(poscent), '_exprMedian', datatype)
    sum_stat <- cbind(sum_stat, poscent)
  }

  if('pmn' %in% moments){
    if(v) cat(' - Mean of positive values\n')
    if('b' %in% moments){
      poscells <- apply(mat, 1, function(vec) tapply(vec, rep("ALL", ncol(mat)), function(cc) mean(cc[cc > 0], na.rm = T) ) )
      sum_stat[, paste0('Bmean_exprFrac', datatype)] <- poscells
    }
    poscent <- apply(mat, 1,function(vec) tapply(vec, groups, function(cc) mean(cc[cc > expr_cutoff], na.rm = T) ) )
    if(length(grps) != 1) poscent <- t(poscent)
    poscent <- data.frame(poscent, check.names = F); # colnames(poscent) <- grps
    colnames(poscent) <- paste0(colnames(poscent), '_exprMean', datatype)
    sum_stat <- cbind(sum_stat, poscent)
  }

  if('q' %in% moments){
    if(v) cat(' - Quantiles\n')
    quants <- cbindList(lapply(c(25, 50, 75), function(myq){
      quant <- apply(mat, 1, function(vec) tapply(vec, groups, quantile, probs = myq / 100, na.rm = T) )
      if(length(grps) != 1) quant <- t(quant)
      quant <- data.frame(quant, check.names = F); # colnames(quant) <- grps
      colnames(quant) <- paste0(colnames(quant), '_Q', myq, datatype)
      return(quant)
    }))
    quants <- quants[, order(colnames(quants))]
    sum_stat <- cbind(sum_stat, quants)
  }

  if(isTRUE(do_log)){
    if(v) cat(' - Mean of shifted log values\n')
    means_log <- apply(log2(mat + 1), 1, function(vec) tapply(vec, groups, mean, na.rm = T) )
    if(length(grps) != 1) means_log <- t(means_log)
    means_log <- data.frame(means_log, check.names = F); # colnames(means_log) <- grps
    colnames(means_log) <- paste0('log2(', colnames(means_log), datatype, '+1)')
    sum_stat <- cbind(sum_stat, means_log)
  }

  if(isTRUE(indiv_values)){
    if(v) cat(' - Adding individual sample values\n   ')
    tvar <- data.frame(rownams = rnames, row.names = rnames, stringsAsFactors = FALSE)
    mat <- as.matrix(mat)
    for(i in grps){ # order sample based on group order specified
      if(v) cat('.')
      tvar <- cbind(tvar, mat[rnames, names(groups[groups %in% i])])
    }; if(v) cat('\n')
    sum_stat <- cbind(sum_stat, tvar[, -1])
  }

  if(v) cat('------- ---- ------- -------\n')
  if(with_rownams) return(sum_stat) else return(sum_stat[, -1])
}
