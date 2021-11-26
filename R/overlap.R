
### Multiple intersections ### -------------------------------------------------
overlap_interset <- function (x){
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], overlap_interset(x[-1]))
  }
}
overlap_unset <- function (x) {
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], overlap_unset(x[-1]))
  }
}
overlap_diffset <- function (x, y) {
  # Remove the union of the y's from the common x's.
  # x and y are lists of characters.
  xx <- overlap_interset(x)
  yy <- overlap_unset(y)
  setdiff(xx, yy)
}

overlap_list <- function(
  groups_list,
  combinations = NULL,
  sep = NULL,
  sharedmax = Inf,
  v = FALSE
){
  if(!is.list(groups_list)){ stop('No list given') }
  if(!is.null(names(groups_list))){
    grps <- names(groups_list)
  }else{
    grps <- paste0('X', 1:length(groups_list))
    names(groups_list) <- grps
  }
  if(v) cat('Names:', paste0(grps, collapse = ', '),'\n')
  if(length(grps) == 1) return(groups_list)
  combinations <- overlap_combn(grps = grps, combinations = combinations, sep = sep, sharedmax = sharedmax, v = v)
  if(v) cat('Calculating overlaping elements\n')
  tvar <- lapply(combinations, function(i){ # Get overlapping elements
            overlap_diffset( groups_list[i], groups_list[setdiff(names(groups_list), i)] )
          })
  names(tvar) <- names(combinations)
  return(tvar)
}

# get a list of group combinations
overlap_combn <- function(
  grps,
  combinations = NULL,
  sep = NULL,
  logchar = "<=",
  sharedmax = Inf,
  nameorder = NULL,
  simplify = FALSE,
  v = FALSE
){
  if(is.null(sep)){
    tvar <- c('n', '~', 'u', '_', '-', 'xx')
    sep <- head(tvar[!tvar %in% unlist(strsplit(grps, split = ""))], 1)
  }
  if(is.null(combinations)){
    combinations <- unlist(lapply(1:length(grps),function(j){ # Getting vector of groups per combination
        combn(grps, j, simplify = FALSE)
    }), recursive = FALSE)
    combinations <- lapply(combinations, function(i){
      if(is.null(nameorder)) sort(i) else this_order(i, ref = nameorder)
    })
    names(combinations) <- sapply(combinations, function(i) paste0(i, collapse = sep) )
  }
  if(v) cat('Combinations:', length(combinations),'\n')
  if(length(grps) > sharedmax){
    cat('Maximum number of combinations', sharedmax,'\n')
    expr <- paste("combinations <- combinations[sapply(combinations, length)", logchar, sharedmax, "]")
    eval(expr = parse(text = expr))
  }
  if(isTRUE(simplify)) return(data.frame(t(vlist2df(combinations)), stringsAsFactors = FALSE))
  return(combinations)
}

# Fisher exact text overlap
overlap_test <- function(
  x,
  y = NULL,
  total = NULL,
  names = c("List1", "List2"),
  row_str = c('in_', 'out_'),
  main = NULL,
  return_plot = TRUE,
  verbose = FALSE
){
  suppressPackageStartupMessages(library(gridExtra))
  if(is.list(x)){
    if(!is.null(names(x))) mynames = names(x)[1:2]
    y = x[[2]]
    total = if(length(x) > 2){
      if(is.numeric(x[[3]]) && length(x[[3]]) == 1) x[[3]] else length(x[[3]])
    }else if(is.null(total)){
      warning("Taking list as universe")
      length(unlist(x))
    }else{ total }
    x = x[[1]]
  }
  if(is.null(y)) stop("Provide second list of features")
  if(is.null(total)) stop("Provide number of features in the universe or the list itself")
  mynames[is.na(mynames)] <- paste0("List", 1:2)[which(is.na(mynames))]

  if(verbose) cat("Building object\n")
  go_object <- GeneOverlap::newGeneOverlap(
    listA = x,
    listB = y,
    genome.size = total
  )
  if(verbose) cat("Performing the test\n")
  go_object <- GeneOverlap::testGeneOverlap(go_object)
  if(verbose) cat("Contingency table\n")
  ddf <- data.frame(GeneOverlap::getContbl(go_object)[2:1, 2:1])
  colnames(ddf) <- paste0(row_str, mynames[1])
  rownames(ddf) <- paste0(row_str, mynames[2])
  if(is.null(main)) main = paste0(mynames, collapse = " vs ")
  ft <- fisher.test(ddf); if(verbose) print(ft) # chisq.test(ddf)
  ddfsum <- ddf; ddfsum$total <- rowSums(ddfsum); ddfsum <- rbind(ddfsum, colSums(ddfsum)); rownames(ddfsum)[3] <- 'total'
  dimnames(ddfsum) <- lapply(dimnames(ddfsum), function(x) casefold(gsub("_", " ", x), T) )
  ddf; lout <- list(object = go_object, conttable = ddfsum)

  if(verbose) cat("Grob table\n")
  tg <- gridExtra::tableGrob(ddfsum)
  x_arrgorb <- gridExtra::arrangeGrob(
    grobs =  list(
      grid::textGrob(ft$method),
      grid::textGrob(casefold(gsub("_", " ", main), T)),
      gridExtra::tableGrob(ddfsum),
      grid::textGrob(paste("95 % CI =", paste0(round(ft$conf.int[1:2], 4), collapse = ', '))),
      grid::textGrob(paste("Odds ratio =", formatC(ft$estimate))),
      grid::textGrob(paste0("Alt: ", ft$alternative, ', P-value =', formatC(ft$p.value)))
    ), ncol = 1, heights = c(1/8, 1/8, 2, 1/8, 1/8, 1/8)
  )
  lout$table_grob <- x_arrgorb
  if(!isTRUE(return_plot)){
    grid.draw(x_arrgorb)
  }else if(verbose) cat("You can plot it with 'grid.draw(out$table_grob)'. Add 'plot.new()' for a new page\n")
  lout
}
