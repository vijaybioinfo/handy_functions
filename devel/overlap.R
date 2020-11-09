
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
