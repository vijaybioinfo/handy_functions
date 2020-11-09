#!/usr/bin/R

# Directory check
dircheck <- function(dname, ...){
  if(!grepl("/$", dname) && dir.exists(dname)) dname <- paste0(dname, "/")
  tvar <- !grepl("_$", dname) && !grepl("/$", dname) && !file.exists(dname)
  if(tvar) dname <- paste0(dname, "_")
  if(grepl("\\/$", dname)) dir.create(dname, ...)
  return(dname)
}

# change directory and if doesn't exists create it
setwdc <- function(dname, ...){
  if(is.null(dname)) stop("'dirname' is null")
  if(dname == '') stop("'dirname' needs to be characters")
  dname <- dircheck(dname)
  tvar <- basename(dirname(dname))
  if(tvar == ".") tvar <- basename(getwd())
  if(!tvar == basename(dname)){
    dir.create(dname, ...)
    setwd(dname)
  }else{ warning("parent of 'dirname' is the same: ", dname, " not created") }
}

# get parent directory N
dirnamen <- function(x, n = 1){
  for(i in 1:n) x <- dirname(x)
  return(x)
}

remove.factors <- function (df) {
  for (varnum in 1:length(df)) {
    if ("factor" %in% class(df[, varnum])) {
      df[varnum] = as.character(df[, varnum])
    }
  }
  return(df)
}

# ggplot type colours
gg_color_hue <- function(n) {
  if(length(n) > 1){
    mynames <- n
    n <- length(n)
  }else{ mynames <- NULL }
  hues = seq(15, 375, length = n + 1)
  tmp <- hcl(h = hues, l = 65, c = 100)[1:n]
  if(!is.null(mynames)) names(tmp) <- mynames
  return(tmp)
}

# vector or matrix to colours
# Check https://www.schemecolor.com/
v2cols <- function(
  select,
  sour = NULL,
  cname = NULL,
  uniq = TRUE,
  fw = 'gg', # force wgcna colours
  myseed = 27,
  in_hex = FALSE,
  verbose = FALSE
){
  colorsnogray <- colors()[grep('white|gray|grey', grDevices::colors(), invert = T)]
  if(is.factor(select)) select <- levels(select)
  if(isTRUE(uniq)) select <- unique(select)
  if(is.null(cname)) cname <- 1
  select <- select[!is.na(select)]
  ordselect <- select
  tmp <- select %in% rownames(sour) | select %in% names(sour)
  if(!all(tmp) && !is.null(sour) && is.null(dim(select))){ # if they are not all in (row)names and you get a source
    dimname <- unique(c(names(sour), unlist(dimnames(sour))))
    unkno <- select[!select %in% dimname]
    select <- select[select %in% dimname]
    if(verbose) warning('Names: ', show_commas(unkno),' not in object')
    if(length(select)){
      select <- v2cols(select, sour = sour, cname = cname, uniq = uniq, fw = fw, verbose = verbose)
    }
    if(isTRUE(fw)){
      tvar <- WGCNA::labels2colors(unkno)
    }else if(fw == 'gg'){
      tvar <- gg_color_hue(length(unkno))
    }else{
      colorsnogray <- colorsnogray[!colorsnogray %in% select]
      set.seed(myseed); tvar <- sample(colorsnogray, length(unkno))
    }
    names(tvar) <- unkno; tmp <- grepl("white", tvar)
    set.seed(myseed); tvar[tmp] <- sample(colorsnogray, sum(tmp))
    select <- c(select, tvar)
    return(select[ordselect])
  }
  if(length(dim(select)) > 1 && !is.null(sour)){
    warning('Selection is data.frame while source is not NULL')
  }
  if(verbose) cat('Colours from ')
  if((is.null(sour) && length(dim(select)) > 1) || isTRUE(fw)){
    if(verbose) cat('labels2colors\n')
    thiscols <- WGCNA::labels2colors(select)
    tmp <- grepl("white", thiscols)
    set.seed(myseed); thiscols[tmp] <- sample(colorsnogray, sum(tmp))
  }else if(is.null(sour)){
    if(verbose) cat('ggplot\n')
    thiscols <- gg_color_hue(length(select))
  }else if(class(sour) == 'data.frame'){
    if(verbose) cat('data.frame\n')
    thiscols <- sour[select, cname]
  }else if(class(sour) == 'character'){
    if(verbose) cat('vector\n')
    thiscols <- sour[select]
  }
  if(length(dim(thiscols)) < 2){
    if(in_hex) select <- col2hex(select)
    names(thiscols) <- select
  }else{
    if(in_hex) select <- apply(select, 2, col2hex)
    colnames(thiscols) <- colnames(select)
    rownames(thiscols) <- rownames(select)
  }
  return(thiscols)
}


### List/Tables manipulation ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# create matrix with given names
mat_names <- function(rnames, cnames, arrname = NULL){
  x <- list(rnames, cnames)
  tvar <- matrix(nrow = length(x[[1]]), ncol = length(x[[2]]), dimnames = x)
  if(!is.null(arrname)){
    tvar <- replicate(length(arrname), tvar)
    dimnames(tvar)[[3]] <- arrname
  }
  return(tvar)
}

joindf <- function(
  x,
  y,
  keep_from_y = "NULL123",
  type = c("left", "right", "full", "none"),
  verbose = FALSE
){
  type <- match.arg(type)
  # str(x); str(y)
  if(type != "none"){
    x$tmpcol123 <- rownames(x)
    y$tmpcol123 <- rownames(y)
    yvars <- colnames(y)[(!colnames(y) %in% colnames(x)) | colnames(y) %in% keep_from_y]
    xvars <- colnames(x)[!colnames(x) %in% yvars] # exclude them from x if exist
    if(verbose) cat('Keeping in y:', yvars, '\n')
    if(verbose) cat('Droping in x:', ifelse(isTRUE(keep_from_y == 'NULL123'), "None", keep_from_y), '\n')
    if(verbose) cat("Using:", type, "\n")
    eval(expr = parse(text = paste0("functy <- dplyr::", type, "_join")))
    z <- functy(
      x = x[, xvars, drop = FALSE],
      y = y[, c(yvars, 'tmpcol123'), drop = FALSE],
      by = "tmpcol123"
    )
    rownames(z) <- z$tmpcol123
    z <- z[, -which(colnames(z) == "tmpcol123")]
  }else{
    mycnames <- unique(c(colnames(x), colnames(y)))
    z <- data.frame(mat_names(rnames = rownames(x), cnames = mycnames))
    z[rownames(x), colnames(x)] <- x
    icells <- intersect(rownames(x), rownames(y))
    z[icells, colnames(y)] <- y[icells, ]
  }
  z
}

# Collapse column summarising the others
column_collapse <- function(metab, rname = 1, verbose = FALSE){
  metab[, rname] <- as.character(metab[, rname])
  if(!any(duplicated(metab[, rname]))){
    if(verbose > 0) cat("."); return(metab) # return if already
  }; if(verbose > 0) cat("s")
  gnames <- unique(metab[, rname][duplicated(metab[, rname])])
  colsum <- colnames(metab)[!colnames(metab) %in% rname] # columns to summarise
  # check if column is same length
  csize <- sapply(colsum, function(co){ ttab <- table(metab[, c(co, rname)]); nrow(ttab) == ncol(ttab) })
  cclass <- sapply(metab[, colsum, drop = FALSE], class) == "character" # keep numerics
  colsum <- names(which(!(csize & cclass))) # eliminate repeated
  if(verbose > 1) cat("\nSummarising:", paste0(colsum, collapse = ", "), "\n")
  # myformula <- paste(rname, "~", paste(colsum, collapse = " + "))
  # y <- dcast.data.table(setDT(metab), myformula)
  if(verbose > 0) cat(length(gnames))
  if(verbose > 1) cat(" elements are repeated\n")
  sx <- data.frame(data.table::rbindlist(lapply(gnames, function(g){ # g <- gnames[1]; g <- "'CISH"
    per_g <- lapply(colsum, function(cname){ # cname <- colsum[1]
      y <- metab[which(metab[, rname] == g), cname]
      y <- if(is.numeric(y)) round(y, 4) else y
      return(paste0(sort(unique(y)), collapse = " & "))
    })
    per_g <- data.frame(per_g, check.names = FALSE, stringsAsFactors = FALSE)
    colnames(per_g) <- colsum
    per_g
  })), check.names = FALSE, stringsAsFactors = FALSE, row.names = gnames)
  sx[, rname] <- gnames
  sx <- sx[, c(rname, colnames(sx)[-ncol(sx)])]
  sx <- list(sx, metab[!metab[, rname] %in% gnames, colnames(sx)])
  sx <- data.frame(data.table::rbindlist(sx, use.names = TRUE, fill = TRUE), check.names = FALSE, stringsAsFactors = FALSE)
  sx
}

# combine tables by a column
summ_tables <- function(
  ltabs,
  column2row = "gene_name",
  tab4roworder = NULL,
  colgroup_trans = NULL,
  colexclude = NULL,
  verbose = FALSE
) {
  # ltabs <- list(res, statseffclusters, res_seu)
  # colgroup_trans = list(cluster = newlabs)
  if(verbose) cat("Columns:", show_commas(sapply(ltabs, ncol)),
  "\nRows:", show_commas(sapply(ltabs, nrow)), "\n")
  if(!is.null(colexclude)){
    colexclude <- paste0(colexclude, collapse = "|")
    if(verbose) cat("Excluding patterns in columns", colexclude, "\n")
  }
  ltabs <- lapply(ltabs, remove.factors)
  if(verbose) cat("Adding", column2row, "from rows when needed: ")
  ltabs <- lapply(ltabs, function(x){
    if(!column2row %in% colnames(x)){
      if(verbose) cat("r"); x[, column2row] <- rownames(x)
    }else if(verbose) cat(".")
    x[, column2row] <- sub("'", "", x[, column2row])
    if(!is.null(colexclude)) x <- x[, colnames(x)[!grepl(colexclude, colnames(x))], drop = FALSE]
    x
  }); if(verbose) cat("\n")
  if(!is.null(colgroup_trans)){
    if(verbose) cat("Changing groups in column(s):", show_commas(names(colgroup_trans)), "\n")
    ltabs <- lapply(ltabs, function(x){
      cnames <- colnames(x)[colnames(x) %in% names(colgroup_trans)]
      for(cname in cnames){
        x[, cname] <- unname(colgroup_trans[[cname]][as.character(x[, cname])])
      }
      x
    })
  }
  if(verbose) cat("Removing repeated column names\n")
  for(i in 2:length(ltabs)){
    cnames <- unique(unlist(sapply(ltabs[1:(i - 1)], colnames)))
    tvar <- colnames(ltabs[[i]]) %in% cnames & !colnames(ltabs[[i]]) %in% column2row
    cat("Table", i, "- repeated:", sum(tvar), "\n")
    ltabs[[i]] <- ltabs[[i]][, !tvar]
  }
  if(verbose) cat("Summarise tables with repeated", column2row, ": ")
  ltabs_sum <- lapply(ltabs, function(mytab){ # mytab <- ltabs[[1]]
    column_collapse(metab = mytab, rname = column2row, verbose = verbose)
  }); if(verbose) cat("\n")
  # Fix the rownames
  # sapply(ltabs_sum, function(x) head(rownames(x)) )
  if(verbose) cat("Setting", column2row, "as rows\n")
  ltabs_sum <- lapply(ltabs_sum, function(x){ rownames(x) <- x[, column2row]; x })
  if(verbose) cat("Final columns:", show_commas(sapply(ltabs_sum, ncol)),
  "\nFinal rows:", show_commas(sapply(ltabs_sum, nrow)), "\n")
  # bind the tables
  cnames <- unique(c(column2row, unlist(sapply(ltabs_sum, colnames))))
  tab4roworder <- unique(c(tab4roworder, 1:length(ltabs_sum)))
  rnames <- unique(unlist(lapply(ltabs_sum[tab4roworder], rownames)))
  if(verbose) cat("Columns, rows: ", length(cnames), ", ", length(rnames), "\n", sep = "")
  summtab <- data.frame( mat_names(rnames = rnames, cnames = cnames),
  stringsAsFactors = FALSE, check.names = FALSE)
  for(x in ltabs_sum){
    ctmp <- colnames(x)[colnames(x) %in% cnames]
    summtab[rownames(x), ctmp] <- x[, ctmp]
  }
  if('columns' %in% names(colgroup_trans)){
    tvar <- unname(colgroup_trans[['columns']][cnames])
    colnames(summtab) <- ifelse(is.na(tvar), cnames, colgroup_trans[['columns']][cnames])
  }
  summtab
}
vlist2df <- function(x, maxlen = NULL){
  if(!is.list(x)) stop('No list given')
  if(is.data.frame(x)) return(x)
  if(is.null(maxlen)) maxlen <- seq_len(max(sapply(x, length)))
  tvar <- data.frame(sapply(x, "[", i = maxlen), stringsAsFactors = F)
  colnames(tvar) <- names(x)
  return(tvar)
}


### Cosmetics ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show_found <- function(x, y, element = 'Elements', verbose = FALSE){
  if(verbose && sum(duplicated(x))) cat('Duplicated:', show_commas(x[duplicated(x)]), '\n')
  if(verbose) cat('Finding ', element, ': ', show_commas(x), '\n', sep = '')
  tmp <- !x %in% y
  if(sum(tmp) && verbose) cat('Missing:', show_commas(x[tmp]), '\n')
  x <- x[!tmp]
  if(verbose) cat('Returning', length(x), casefold(element), '\n')
  return(x)
}

show_commas <- function(x, hn = 3){
  tvar <- length(x)
  if(tvar > 1){
    tmp <- paste0(head(x[-tvar], hn), collapse = ', ')
    connect <- ifelse(tvar-1 > hn, ' ... and ', ' and ')
    tmp <- paste0(tmp, connect, x[tvar])
    if(tvar-1 > hn) tmp <- paste0(tmp, ' (', tvar, ')')
  }else{
    tmp <- x
  }
  return(tmp)
}

show_parameters <- function(
  parms,
  flaggy = NULL, # Logical vector indicating whethe the parameter is default
  fast_track = FALSE # Print variable assignments to be copied in case you want
  # to run the code manually
){
  cat('Parameters:\n')
  if(is.null(names(parms))) stop('No parameter\'s names in vector\n')
  if(is.null(flaggy) && length(flaggy) == 0){
    flaggy = rep(TRUE, length(parms)); names(flaggy) <- names(parms)
    defv <- c('NA', 'NA')
  }else{ defv <- c('T', 'F') }
  maxchar <- max(nchar(parms))
  cat(paste0(rep(' ',maxchar),collapse = ''),'Default\tPARAMETERS\n')
  void <- lapply(1:length(parms),function(i){
    tvar <- parms[i]
    addspaces <- paste0(rep(' ', maxchar-nchar(tvar)),collapse = '')
    lol <- suppressWarnings(try(evaluation <- paste0(eval(as.name(names(tvar))), collapse=', ')))
    if(class(lol) == 'try-error'){
      cat('lol')
      assign(names(tvar),eval.parent(as.name(names(tvar)),n=3),pos=-1)
      evaluation <- paste0(eval(as.name(names(tvar))), collapse=', ')
    }
    ifdef <- ifelse(flaggy[names(tvar)],  defv[1], defv[2])
    cat(tvar, addspaces, ifdef,'\t--->',evaluation,'\n')
    if(is.character(eval(as.name(names(tvar))))[1]){
      if(length(eval(as.name(names(tvar)))) > 1){
        return(paste0(names(tvar)," = c(\'",paste0(eval(as.name(names(tvar))),collapse='\', \''),'\')\n'))
      }else{
        return(paste0(names(tvar)," = \'",paste0(eval(as.name(names(tvar))),collapse='\', \''),'\'\n'))
      }
    }else{
      if(length(eval(as.name(names(tvar)))) > 1){
        return(paste0(names(tvar)," = c(",paste0(eval(as.name(names(tvar))),collapse=', '),')\n'))
      }else{
        return(paste0(names(tvar)," = ",eval(as.name(names(tvar))),'\n'))
      }
    }
  })
  if(fast_track){ cat('\nFast track:\n'); void <- lapply(void, function(i){ cat(' ', i) }) }
  cat('\n*******************************************************************\n')
}

# modify string to create new lines for fitting long strings
newlines <- function(stringy, ln = 6, sepchar = ' '){
  stringy <- as.character(stringy); if(!grepl(sepchar, stringy)) return(stringy)
  if(Hmisc::stringDims(stringy)$width > ln){
    vec <- unlist(strsplit(stringy, sepchar))
    ln <- ifelse((ln + 1) > length(vec), length(vec) - 1, ln)
    newvec <- vector(); ivec <- seq(ln + 1, length(vec), ln)
    for(i in ivec){
      newelem <- paste(vec[(i-ln):(i-1)], collapse = sepchar)
      addline <- sepchar == '' && !grepl('.* $', newelem) && !grepl('^ .*', vec[i])
      if(addline) newelem <- paste0(newelem, '-')
      newelem <- sub('^ +', '', sub(' +$', '', newelem))
      newvec <- c(newvec, newelem)
    }
    if(i <= length(vec)){
      newelem <- paste0(vec[i:length(vec)], collapse = sepchar)
      newelem <- sub('^ +', '', sub(' +$', '', newelem))
      newvec <- c(newvec, newelem)
    }
    stringy <- paste(newvec, collapse = '\n')
  }
  return(stringy)
}

# add needed spaces to reach a certain number of characters
addspaces <- function(x, m) paste0(x, paste0(rep(' ', m-nchar(x)), collapse = ''))

# Remove 'boring' interesting features
cosmy <- function(x, patties = "^rps|^rpl|^mt-|rp[0-9]{1,}\\-"){
  tvar <- grepl(patties, casefold(features_parse_ensembl(x)))
  if(sum(tvar)) x <- x[!tvar]
  return(x)
}

# head for matrices
headmat <- function(x, n = 5) head(x[, head(colnames(x), n)], n)
tailmat <- function(x, n = 5) tail(x[, tail(colnames(x), n)], n)
headtail <- function(x, n = 5){
  if(is.null(dim(x))) return(c(head(x, n), tail(x, n)))
  rbind(head(x, n), tail(x, n))
}
bordering <- function(dtab, cnames = 1, n = 5, func = headtail){
  y <- unique(unlist(lapply(cnames, function(x){
    rownames(dtab[func(order(dtab[, x]), n), ])
  })))
  y[y %in% rownames(dtab)]
}

is.file.finished <- function(x, size = 3620) file.exists(x) && isTRUE(file.size(x) > size)

# create factors with levels in the right order
factormix <- function(x){
  if(!is.character(x)) return(x)
  y <- as.character(x)
  factor(y, levels = gtools::mixedsort(unique(y)))
}

### Order ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# order vector (when repeated values) following a reference order
order_vector <- function(x, ref, which = FALSE){
  y <- unique(x)
  order_all <- c(ref[ref %in% y], y[!y %in% ref])
  cat(paste0(order_all, collapse = ", "), "\n")
  y <- c()
  if(isTRUE(which)){
    for(i in order_all) y <- c(y, which(x %in% i))
  }else{
    for(i in order_all) y <- c(y, x[x %in% i])
  }
  return(y)
}
# order based on column
order_df <- function(
  annot,
  order_by = 'corder',
  cname = NULL,
  grupos = NULL,
  mat = NULL,
  lgenes = NULL,
  samplit = FALSE,
  verbose = FALSE
){
  annot$corder <- 'ALL'
  if(!order_by %in% c('pca', 'mean', 'hc', colnames(annot))) order_by <- "corder"
  if(verbose) cat('Setting columns order by:', show_commas(order_by), '\n')
  if(order_by == 'pca' && is.null(mat)){
    stop("You need an expression matrix, 'mat'")
    return(rownames(annot))
  }
  if(is.null(cname)) cname <- order_by[1]
  if(verbose) cat('Covariates:', show_commas(cname, Inf), '\n')
  minsiz <- min(table(annot[, cname]))
  if(isTRUE(samplit)){
    if(verbose) cat("Sampling\n")
    if(!is.numeric(samplit)) samplit <- NULL else samplit <- -samplit
    annot <- annot[sample_grp(annot, cname = cname, maxln = samplit), ]
  }
  grupos <- grupos[grupos %in% as.character(annot[, cname])]
  if(length(grupos) == 0){ warning("Given groups not in data"); grupos <- NULL }
  if(is.null(grupos)) grupos <- gtools::mixedsort(unique(as.character(annot[, cname])))
  if(verbose) cat('Groups:', show_commas(grupos, Inf), '\n')
  if(is.null(lgenes) && !is.null(mat)){
    lgenes <- rep(list(rownames(mat)), length(grupos))
    names(lgenes) <- grupos
  }
  list_order <- lapply(grupos, function(x){
    mynames <- filters_subset_df(c(cname, x), annot, verbose = verbose)
    if(verbose) cat(x, "")
    myorder <- ifelse(isTRUE(order_by %in% colnames(annot) && order_by != "corder"), 'acol', order_by)
    if(verbose && (!myorder %in% c("acol", "corder"))) cat(myorder); if(verbose) cat("\n")
    if(cname == order_by) return(mynames)
    if(!is.null(mat)) mygenes <- show_found(lgenes[[x]], rownames(mat), v = FALSE)
    mynames <- switch(myorder,
      acol = {
        if('corder' %in% order_by){
          torder_by <- colnames(annot)[grepl('^orig\\.|^res', colnames(annot))]
        }else{ torder_by <- order_by }
        tannot <- order(annot[mynames, torder_by], decreasing = TRUE)
        str(mynames[tannot])
        mynames[tannot]
      }, pca = {
        smat <- mat[mygenes, mynames]
        pc1 <- irlba::prcomp_irlba(t(smat), n = min(c(3, dim(smat) - 1)) )$x[, 1]
        names(pc1) <- mynames; rev(names(sort(pc1, decreasing = T)))
      }, hc = {
        smat <- t(mat[mygenes, mynames])
        hc1 <- hclust(dist(smat, method = 'euclidean'), method = 'average')
        hc1$labels[hc1$order]
      }, mean = {
        names(sort(colMeans(as.matrix(mat[mygenes, mynames])), decreasing = T))
      },
      corder = mynames
    )
    return(mynames)
  })
  list_order
}
# to get a df in the orther of a specified vector
order_df1 <- function(ids, df, cname){
  y <- unlist(sapply(ids, function(x) which(df[, cname] == x) ))
  nrows <- 1:nrow(df)
  if(length(y) != nrow(df)) y <- c(y, nrows[!nrows %in% y])
  return(y)
}
# now order a vector
order_group <- function(x, groups){
  groups <- as.character(groups)
  x <- as.character(x)
  c(groups[groups %in% x], x[!x %in% groups])
}

# get columns with a max N of groups
extract_grouping_cols <- function(
  metadat,
  onames = NULL,
  maxn = 27,
  ckeep = NULL, cignore = NULL,
  types = "character",
  na_rm = FALSE,
  verbose = FALSE
){
  metadat <- remove.factors(metadat)
  if(is.null(onames)) onames <- colnames(metadat)
  if('numeric' %in% types) maxn <- Inf
  onames <- rev(onames)
  if(verbose) cat("Selecting types:", paste0(types, collapse = ", "), '\n')
  tvar <- sapply(metadat[, onames, drop = FALSE], class) %in% types # if it's the type and has no NA's
  if(isTRUE(na_rm)) tvar <- tvar & !sapply(metadat[, onames, drop = FALSE], function(x) any(is.na(x)) )
  onames <- onames[unname(tvar)]
  if(!is.null(ckeep)){
    if(length(ckeep) > 1){
      onames <- onames[onames %in% ckeep]
    }else{ onames <- onames[grepl(ckeep, onames, ignore.case = TRUE)] }
  }
  if(!is.null(cignore)) onames <- onames[!grepl(cignore, onames, ignore.case = TRUE)]
  if(length(onames) == 0){ if(verbose) cat("Taking all columns\n"); onames <- colnames(metadat) }
  if(verbose) cat("Filtering", length(onames), "\n")
  nnames <- sapply(onames, function(x) length(unique(metadat[, x])) )
  if(verbose) cat(sum(duplicated(nnames)), "same N across columns\n")
  onames <- names(sort(nnames))
  onames <- onames[nnames[onames] != nrow(metadat)] # excluding  N = total rows
  onames <- onames[nnames[onames] > 1] # excluding N = 1
  onames <- onames[nnames[onames] <= maxn] # N <= 27 because ggpairs doesn't allow more than this
  if(length(onames) == 1) return(onames)
  if(verbose) cat("Match 2-group combination:", length(onames), "\n")
  pairs <- try(gtools::combinations(n = length(onames), r = 2, v = onames))
  if(class(pairs) != 'try-error'){
    kha <- apply(pairs, 1, function(x){
      y <- table(metadat[, c(x)]) # as.data.frame.matrix
      diag(y) <- 0; any(rowSums(y) > 0) # take when none but diagonal are 0, aka, same N groups
    })
    if(verbose && any(!kha)) cat(sum(!kha), "matching 1-to-1 groups:", paste0(pairs[!kha, 2], collapse = ", "), "\n")
    onames <- onames[!onames %in% pairs[!kha, 2]]
  }else{
    if(verbose) cat("n =", length(onames), "/ unique =", length(unique(onames)), "-", paste0(unique(onames), collapse = ', '), "\n")
  }
  kha <- apply(t(apply(metadat[, onames, drop = FALSE], 1, duplicated)), 2, sum) # get dups across columns
  onames <- onames[(kha != nrow(metadat))]
  if(verbose) cat("Returning:", length(onames), "\n")
  onames
}

# make a list of rownames out a column's categories or
# return categories for each group in another column
make_list <- function(x, colname = colnames(x)[1], col_objects = NULL, grouping = FALSE){
  x <- remove.factors(x)
  tvar <- lapply(unique(x[, colname]), function(y){
    rnames <- filters_subset_df(list(c(colname, y)), x, v = F)
    if(is.null(col_objects)) return(rnames)
    x[rnames, col_objects]
  })
  names(tvar) <- unique(x[, colname])
  if(isTRUE(grouping)){
    tmp <- rep(names(tvar), sapply(tvar, length))
    names(tmp) <- unlist(tvar); tvar <- tmp
  }
  return(tvar)
}

### Features handling ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Eliminating ensmbl number
features_parse_ensembl <- parse_ens_name <- function(x, keepens = FALSE){
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

# find genes even if they have ensembl ids
features_matching <- findgenes <- function(x, y){
  if(all(x %in% y)) return(x[x %in% y])
  tvar <- grepl(paste(paste0("^", parse_ens_name(x), "$"), collapse = "|"), parse_ens_name(y), ignore.case = TRUE)
  y[tvar]
}

### Counts transformations ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count_transformation <- function(
  cts,
  transf = "log2",
  pcount = 1,
  norm_fact = 1e+6,
  verbose = FALSE
){
  cts <- as.matrix(cts)
  dtype <- casefold(sub("log2", "", transf), upper = TRUE) # extract transformation type
  if(grepl("cpm", casefold(transf))){
    if(verbose) cat('Counts to CPM\n')
    sweep(cts, 2, colSums(cts), '/') * norm_fact
  }
  if(dtype == ""){ # if no info, try to identify the type of data
    if(verbose) cat('Detecting type\n')
    cts_mean <- mean(colSums(cts))
    if(verbose) cat('Mean total sum:', cts_mean, '\n')
    csums <- which(c(1e6, 1e4) %in% cts_mean) # per million or 10,000 (used by Seurat)
    if(length(csums) == 0){ csums = "sad"; dtype = "CTS" }
    dtype <- switch(csums, "1" = "CPM", "2" = "CP10K", dtype)
  }

  if(isTRUE(grepl('log', transf))){
    bs <- as.numeric(sub('log', '', transf))
    dtype <- paste0("log", bs, "(", dtype, " + ", pcount, ")")
    if(verbose) cat(dtype, 'transformation\n')
    if(is.na(bs)) bs <- exp(1)
    cts <- log(cts + pcount, bs)
  }else if(isTRUE(transf == 'vst')){
    if(verbose) cat('Variance Stabilizing transformation\n')
    cts <- varianceStabilizingTransformation(mat)
  }else{
    if(verbose) cat('No log transformation performed\n')
  }

  if(verbose) cat('Type of counts:', dtype, '\n')
  return(list(data = cts, type = dtype))
}

### Plotting ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theme_axes <- function (
  base_size = 11,
  base_family = "",
  base_line_size = base_size/22,
  base_rect_size = base_size/22
){
  half_line <- base_size/2
  t <- theme(
    line = element_line(
      colour = "black", size = base_line_size,
      linetype = 1, lineend = "butt"), rect = element_rect(fill = "white",
      colour = "black", size = base_rect_size, linetype = 1
    ),
    text = element_text(
      family = base_family, face = "plain",
      colour = "black", size = base_size, lineheight = 0.9,
      hjust = 0.5, vjust = 0.5, angle = 0, margin = margin(),
      debug = FALSE
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "none",
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    strip.text = element_blank(),
    strip.background = element_blank()
  )
  # ggplot_global$theme_all_null %+replace% t
}
plot_rm_layer <- function(gp, lays = "ext", verbose = FALSE){
  if(verbose) cat("Layers:", sapply(gp$layers, function(x) class(x$geom)[1] ), "\n")
  tvar <- sapply(gp$layers, function(x) !grepl(lays, class(x$geom)[1], ignore.case = TRUE) )
  gp$layers <- gp$layers[tvar]
  gp
}
plot_blank <- function(gp, lays = "ext", ...){
  plot_rm_layer(gp, lays = lays, ...) + theme_axes()
}

plot_flush <- function(x){ graphics.off(); invisible(file.remove('Rplots.pdf')) }
