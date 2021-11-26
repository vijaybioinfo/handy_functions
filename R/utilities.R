#!/usr/bin/R

# Directory check
dircheck = prefix_check = function(
  path, sufix = "_", ...
){
  if(!grepl("\\/$", path) && dir.exists(path)) path <- paste0(path, "/")
  if(!is.null(sufix)){
    tvar <- !grepl(paste0(sufix, "$"), path) &&
      !grepl("\\/$", path) && !file.exists(path)
    if(tvar) path <- paste0(path, sufix)
  }
  if(grepl("\\/$", path)){
    dir.create(path = path, ...)
  }else if(!dir.exists(dirname(path))) dir.create(path = dirname(path), ...)
  return(path)
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

# Archive files/folders
file.archive = function(pattern, name = "archive/", exclude = NULL, type = "mv"){
  # format(Sys.time(), '%Y_%m_%d/')
  # type = match.arg(type)
  outdir = paste0(name, Sys.Date(), "/")
  dir.create(outdir, recursive = TRUE)
  y <- if(any(file.exists(pattern))){
    list.files(
      path = dirname(pattern),
      pattern = basename(pattern),
      full.names = TRUE, all.files = TRUE
    )
  }else{ pattern }
  if(!is.null(exclude)) y <- y[!grepl(pattern = exclude, x = y)]
  if(length(y)){ # wil try to mirror the folder tree
    command = paste0(type, " ", y, " ", outdir, dirname(y), "/"); cat(command, sep = "\n")
    # ask <- if(interactive()) readline("Press 1[ENTER] to move file(s): ") else 1
    cat("Press 1[ENTER] to move file(s): ")
    ask <- if(interactive()) readline("") else readLines("stdin", n = 1)
    if(ask == 1){
      cat("\n"); lapply(1:length(y), function(x){
        outdir_i = paste0(outdir, dirname(y[x]), "/")
        if(!dir.exists(outdir_i)) dir.create(outdir_i)
        command_i = paste(type, y[x], outdir_i); cat(command_i, sep = "\n")
        system(command_i)
      })
    }
  }; invisible(x = NULL)
}

remove.factors <- function (df) {
  for (varnum in 1:length(df)) {
    if ("factor" %in% class(df[, varnum])) {
      df[varnum] = as.character(df[, varnum])
    }
  }
  return(df)
}

# https://stackoverflow.com/questions/3452086/getting-path-of-an-r-script
script_rscript <- function() {
  rscript_file <- gsub(".*=(.*)", "\\1", grep("--file=", base::commandArgs(), value = TRUE))
  return(rscript_file)
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

# Summarises a table by the categories of a column
summarise_table = function(
  x,
  column = colnames(x)[1],
  sep = "|",
  numeric_fun = mean,
  expand = FALSE,
  verbose = FALSE
){
  if(verbose) str(x)
  if(verbose) cat("Column:", column, "\n")
  x[, column] <- factor(x[, column])
  x <- x[, c(column, setdiff(colnames(x), column))]
  summarised <- if(ncol(x) == 1){
    return(x)
  }else{
    if(verbose) cat("Summarising:", nlevels(x = x[, column]), "elements\n")
    y <- lapply(
      X = setNames(nm = levels(x = x[, column])),
      FUN = function(ident) {
        slice <- x[which(x[, column] %in% ident), , drop = FALSE];
        if(verbose > 1) cat(" -", ident, "\n")
        # Collapsed into the same columns
        z <- if(!expand){
          temp <- try(setNames(2:ncol(slice), colnames(slice)[-1]))
          if(class(temp) == 'try-error'){ str(slice); stop(ident, " failed") }
          data.frame(lapply(
            X = temp,
            FUN = function(i){
              y <- slice[, i]; if(is.factor(y)) y <- droplevels(y)
              if(!is.numeric(y)) paste0(levels(factormix(y)), collapse = sep) else numeric_fun(y)
          }));
        }else{ # Expand columns
          props <- lapply(slice[, -1, drop = FALSE], table)
          z <- data.frame(t(reshape2::melt(props)), stringsAsFactors = FALSE)
          colnames(z) <- unlist(lapply(props, names))
          z$total <- nrow(slice); z[2, ]
        }; z[, column] <- ident; z
    }); if(verbose) cat("\n");
    y <- as.data.frame(data.table::rbindlist(y, fill = TRUE))
    y <- y[, c(column, setdiff(colnames(y), column))];
    rownames(y) <- as.character(y[, column])
    if(expand) for(i in colnames(y)[-1]) y[, i] <- as.numeric(y[, i])
    for(i in colnames(x))
      if(is.factor(x[, i]) && (i %in% colnames(y)))
        y[, i] <- factor(as.character(y[, i]), levels(x[, i]))
    y
  }; return(summarised)
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
vlist2df_diff <- function(x, y, delim = "Not found"){
  final <- lapply(names(x), function(z){
    x1 <- !y[[z]] %in% x[[z]] # values missing in the subset
    if(sum(x1) > 0) c(x[[z]], delim, y[[z]][x1]) else x[[z]]
  }); ddf <- vlist2df(final)
  sizes <- paste0(sapply(x, length), "of", sapply(y, length))
  names(ddf) <- paste0(names(x), ";", sizes)
  ddf[is.na(ddf)] <- ""
  return(ddf)
}

# process list to get data.frames
list2dfs <- function(mylist){
  if(all(sapply(mylist, is.vector))){
    mylist <- list(vlist2df(mylist))
  }else{
    mylist <- lapply(mylist, function(x){
      if(is.data.frame(x)) return(x)
      if(is.vector(x)) return(data.frame(x, stringsAsFactors = F))
      if(is.matrix(x)){
        tvar <- data.frame(x, stringsAsFactors = F)
        colnames(tvar) <- colnames(x)
        return(tvar)
      }
    })
  }
  mylist
}

# make list of data.frames or list of lists of vectors the same number of rows
list2evendf <- function(mylist){
  mylist <- list2dfs(mylist)
  maxlen <- max(sapply(mylist, nrow))
  mylist <- lapply(mylist, function(z){
      tvar <- mat_names(rnames = nrow(z):(maxlen-1), cnames = colnames(z))
      if(nrow(z) != maxlen) tvar <- rbind(z, tvar) else tvar <- z
      colnames(tvar) <- colnames(z)
      return(tvar)
  })
}

# merge list of data.frames or matrices with differing nrows
cbindList <- function(mylist){
  mylist <- list2evendf(mylist)
  listdf <- do.call(cbind, list2evendf(mylist))
  return(listdf)
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
    add_nspaces <- paste0(rep(' ', maxchar-nchar(tvar)),collapse = '')
    lol <- suppressWarnings(try(evaluation <- paste0(eval(as.name(names(tvar))), collapse=', ')))
    if(class(lol) == 'try-error'){
      cat('lol')
      assign(names(tvar),eval.parent(as.name(names(tvar)),n=3),pos=-1)
      evaluation <- paste0(eval(as.name(names(tvar))), collapse=', ')
    }
    ifdef <- ifelse(flaggy[names(tvar)],  defv[1], defv[2])
    cat(tvar, add_nspaces, ifdef,'\t--->',evaluation,'\n')
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

is.file.finished <- function(x, size = 3620) file.exists(x) & file.size(x) > size

# create factors with levels in alphanumeric order
factormix <- function(x){
  if(!is.character(x)) return(x)
  y <- as.character(x)
  factor(y, levels = gtools::mixedsort(unique(y)))
}
ident_combine = function(mdata, columns = colnames(mdata), sep = "-"){
  for(i in columns) mdata[, i] <- factormix(mdata[, i])
  if(length(columns) > 1){
    tvar <- interaction(mdata[, columns], sep = sep, lex.order = TRUE)
    factor(tvar, levels(tvar)[levels(tvar) %in% as.character(tvar)])
  }else{ mdata[, columns] }
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
        hc1 <- hclust(dist(smat))
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

# make a list of rownames out a column's categories or
# return categories for each group in another column
# check tibble::deframe
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
  pattern <- "X1234" # removing the ensembl name
  if(grepl("ENS", x[1]) && grepl("[0-9]", x[1])){
    ensloc <- all(grepl("^ENS", x[1:10])) # find which side the names are
    if(isTRUE(keepens)) ensloc <- !ensloc
    pattern <- ifelse(ensloc, ".*_", "_.*")
    pattern <- paste0(c(pattern, ifelse(ensloc, ".*\\|", "\\|.*")), collapse = "|")
  }; newnames <- gsub(pattern, "", x); tmp <- grepl("\\-", newnames)
  newnames[tmp] <- gsub("\\-", "_", newnames[tmp]) # so we can keep the dashes!
  newnames <- make.names(newnames, unique = TRUE, allow_ = TRUE)
  newnames[tmp] <- gsub("_", "-", newnames[tmp]) # returning dashes ;)
  newnames
}

# find genes even if they have ensembl ids
features_matching <- findgenes <- function(x, y){
  if(all(x %in% y)) return(x[x %in% y])
  tvar <- grepl(paste(paste0("^", parse_ens_name(x), "$"), collapse = "|"), parse_ens_name(y), ignore.case = TRUE)
  y[tvar]
}

features_find_symbols <- function(
  features_init,
  annotation = NULL,
  gene_name = c("gene_name", "gene"),
  verbose = TRUE
){
  annotation <- if(is.character(annotation)){
    annotation <- annotation[file.exists(annotation) & !dir.exists(annotation)]
    if(verbose) cat("Using", annotation, "\n")
    if(grepl("json$", annotation)){
      config <- rjson::fromJSON(paste(readLines(annotation), collapse=""))
      config <- sapply(config, function(x) if(is.list(x)) x else as.list(x) )
      config$config$annotation_file
    }
  }else{ annotation }
  gannot <- if(is.character(annotation)){
    if(verbose) cat("Reading", annotation, "\n")
    readfile(annotation, stringsAsFactors = FALSE)
  }else{ annotation }
  gene_name <- which(gene_name %in% colnames(gannot))[1]
  if(is.data.frame(gannot)){
    if(verbose) cat("Matching features column:")
    rnames_col = which(sapply(gannot, function(x) all(features_init %in% x) ))[1]
    if(!is.na(rnames_col)){
      if(verbose) cat(rnames_col)
      rownames(gannot) <- gannot[, rnames_col]
      gannot <- gannot[features_init, ]
    }; if(verbose) cat("\n")
  }
  newnames <- if(isTRUE(gene_name %in% colnames(gannot))){
    if(all(features_init == rownames(gannot))){
      if(features_init[1] != gannot[1, gene_name]){
        if(verbose) cat("Appending", gene_name, "\n")
        gannot$feature_id = paste0(features_init, "_", gsub("'", "", gannot[, gene_name]))
        gannot$feature_id
      }else{ features_init }
    }else{ #if(mean(gannot[, gene_name] %in% features_init) > .9){
      gannot = NULL; features_init
    }
  }else{ features_init }
  if(verbose) str(newnames)
  return(list(annotation = gannot, features = features_init))
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
    cts <- sweep(cts, 2, colSums(cts), '/') * norm_fact
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
    if(is.na(bs)) bs <- exp(1)
    dtype <- paste0("log", bs, "(", dtype, " + ", pcount, ")")
    if(verbose) cat('log', bs, 'transformation\n')
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

MASS_kde2d <- function(x, y, ...) {
  band.nrd = MASS::bandwidth.nrd
  dens <- MASS::kde2d(x, y,
    h = c(ifelse(band.nrd(x) == 0, 0.1, band.nrd(x)),
          ifelse(band.nrd(y) == 0, 0.1, band.nrd(y))),
    ...
  ); ix <- findInterval(x, dens$x); iy <- findInterval(y, dens$y)
  return(dens$z[cbind(ix, iy)])
}

values_capped <- function(
  x, value_min = NULL, value_max = NULL,
  quantiles = 0.95, verbose = FALSE
){
  value_ran <- round(range(x, na.rm = TRUE), 3)
  if(is.null(value_min)) value_min <- value_ran[1]
  if(is.null(value_max)) value_max <- quantile(x = x, probs = quantiles, na.rm=TRUE)
  y <- x; y[y < value_min] <- value_min
  y[y > value_max] <- value_max
  caption_this = paste0("Values ranged from ", value_ran[1],
    " to ", value_ran[2], "; capped at ", value_min, " and ",
    quantiles, " (quantile, ", round(value_max, 3), ")")
  cat(caption_this, "\n"); y
}
