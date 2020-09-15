#!/usr/bin/R

###############################
# Handy functions for n-tasks #
###############################

# This functions are designed to be use across all scripts from Vijay's Lab
# NOTE: highly sensitive, dont change input parameters names so confidently

# source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')

# When, for some reason, Ctrl+L stops working
clr <- function() system("clear")

## from counts to counts per factor (1M)
cts2cpm <- function(x, fact = 1e+6){
  x <- as.matrix(x)
  sweep(x, 2, colSums(x), '/') * fact
}
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
    cts <- cts2cpm(x = cts, fact = norm_fact)
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

# Converting RPKM to TPM
fpkmToTpm <- function(fpkm){
  # https://haroldpimentel.wordpress.com/2014/05/08/what-the-fpkm-a-review-rna-seq-expression-units/
  # exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
  sweep(fpkm, 2, colSums(fpkm), '/' ) * 1e6 # from Ariel
}

# Call when a process is taking too long and you wan to be notified when it finishes
process_has_ended <- function(x = NULL, email = "${USER}"){
  x <- ifelse(is.null(x), "Process", x)
  command <- paste0('echo "', x,  ' finished" | mail -s "Job status" "', email, '@lji.org"')
  cat(command, '\n')
  tvar <- system(command, intern = TRUE)
}
# process_has_ended()

# Directory check
dircheck <- function(dname, appe = NULL){
  if(sum(dname == '') || is.null(dname)){
    warning('No path given'); return(dname)
  }
  if(!grepl("\\/$", dname)) dname <- paste0(dname, "/")
  if(is.character(appe)) dname <- sub(".$", appe, dname)
  sub("_{2,}", "_", sub("\\/{2,}", "/", dname))
}
ll <- function(x = getwd()) system(paste('ls -hlt', x))

# get a name or a "blank" character ('')
newname <- function(x, default = '') ifelse(x != default, x, '')

remove.factors <- function (df) {
  for (varnum in 1:length(df)) {
    if ("factor" %in% class(df[, varnum])) {
      df[varnum] = as.character(df[, varnum])
    }
  }
  return(df)
}

# String comparison
compareStr <- function(list1, list2, threshold = 0.9, v = FALSE){
  if(v) cat('--- String comparison in progress ---\n')
  if(v) cat('Using grepl()\n')
  tvar <- sapply(list1, function(x){
    tmp <- which(grepl(x, list2))
    if(length(tmp) > 1) cat(x, ' duplicated\n')
    list2[tmp][1]
  })
  sppDF <- data.frame(list1 = names(tvar), list2 = tvar, score = 1,
    stringsAsFactors = F)
  sppDF <- sppDF[complete.cases(sppDF), ]
  if(nrow(sppDF) != length(list1)){
    if(v) cat('Trying Levenshtein Similarity\n')
    sppDF <- data.table::rbindlist(lapply(list1, function(x){
      data.table::rbindlist(lapply(list2, function(sp){
        score <- RecordLinkage::levenshteinSim(as.character(x), as.character(sp))
        if(is.na(score)) score <- -Inf
        if(score > threshold){
          data.frame('list1' = x, 'list2' = sp, 'score' = score)
        }
      }))
    }))
  }
  tvar <- sum(duplicated(sppDF$list1))
  if(tvar > 0){
    if(v) cat('Correcting for duplicated -',tvar,'\n')
    sppDF <- sppDF[!duplicated(sppDF$list1, fromLast = T), ]
  }
  if(v) cat(nrow(sppDF), 'matches found.\n')
  return(sppDF)
}

# Load packages and install them if necessary #
load_packs <- function(
    packs,
    install = FALSE,
    silence = FALSE,
    lib = .libPaths(),
    v = TRUE, vpb = FALSE
  ){
  if(v) cat('*******************************************************************\n')
  .libPaths(new = lib[1])
  if(v) cat('Library:', .libPaths()[1], '\n')
  if(v){ cat('Loading', length(packs),'package(s)\n** '); cat(packs, sep='\n** ') }
  if(v && vpb) pb <- txtProgressBar(min = 0, max = length(packs), width = NA,  style = 3)
  deps <- vector()
  for(i in 1:length(packs)){
    if(v) cat('.'); if(v && vpb) setTxtProgressBar(pb, i)
    deps <- append(deps, suppressMessages(require(packs[i], character.only = T)))# lib.loc = lib[1])))
  }; if(v) cat('\n'); if(v && vpb) close(pb)
  if(v) cat("  |=======================================================| 100%\n")
  if(v) cat(length(packs[deps]), 'package(s) loaded\n')
  if(length(packs[deps]) < length(packs) && install){ # To check missing packages
    packs <- packs[!deps]; if(v) cat(length(packs), 'not loaded:\n** ')
    if(v){ cat(packs, sep = '\n** '); cat('\nInstalling...\n') }
    deps <- lapply(packs,function(x){
        if(v) cat('My turn:', x, '\n')
        install.packages(x, repos = 'https://cloud.r-project.org/', quiet = silence)#, lib = lib[1])
        instp <- try(library(x, character.only = T, quietly = silence))# lib.loc = lib[1]))
        if(class(instp) == 'try-error'){
          if(v) cat('---------------\n')
          if(silence){
            suppressMessages(source('https://bioconductor.org/biocLite.R'))
            suppressMessages(instp <- biocLite(x, suppressUpdates = TRUE))
          }else{
            source('https://bioconductor.org/biocLite.R')
            instp <- biocLite(x, suppressUpdates = TRUE)
          }
        };if(v) cat('\n')
    })
    deps<-unlist(lapply(packs,function(x){ suppressMessages(require(x, character.only=T)) }))
    if(v){ cat(length(packs[deps]),'package(s) installed and charged:\n** ')
           cat(packs[deps], sep='\n** ') }
  }
  if(v) cat('*******************************************************************\n')
}

# Present parameters
present_params <- function(
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

# Get the object saved in a RData file
theObjectSavedIn <- function(
  saveFile,
  ob_name = NULL,
  v = FALSE
) {
  env <- new.env()
  if(v) cat('Reading object ')
  if(grepl('rds$', saveFile)){
    if(v) cat("RDS\n")
    env$rds <- readRDS(saveFile)
  }else{ if(v) cat("RDATA\n"); load(saveFile, envir = env) }
  loadedObjects <- objects(env, all = TRUE)
  if(length(loadedObjects) == 1){
    ob_name <- loadedObjects
  }else if(is.null(ob_name)){
    if(length(loadedObjects) > 1) stop(length(loadedObjects), ' objects in data, specify one')
    ob_name <- loadedObjects
  }else if(is.numeric(ob_name)){ # mmda para hacerlo bonito
    if(ob_name > length(loadedObjects)){ # if there fewer objects than the asked number
      if(grepl('1$', as.character(ob_name))) tvar <- paste0(ob_name, 'st')
      if(grepl('2$', as.character(ob_name))) tvar <- paste0(ob_name, 'nd')
      if(grepl('3$', as.character(ob_name))) tvar <- paste0(ob_name, 'rd')
      if(sum(c('12', '13', '11') %in% as.character(ob_name))) tvar <- paste0(ob_name, 'th')
      if(!exists('tvar')) tvar <- paste0(ob_name, 'th')
      stop(length(loadedObjects), ' objects in data, ', tvar, ' not found')
    }
    ob_name <- loadedObjects[ob_name]
  }
  if(v) cat('Taking:', ob_name, '\n')
  if(!sum(ob_name %in% loadedObjects)) stop('Object not found')
  env[[ob_name]]
}

# Save history automatically
.Last <- function() {
  graphics.off()
  # system("screen -ls | grep Attached | grep -o '\\.[a-z]*' | grep -o '[a-z]*'", intern = TRUE)
  fname <- gsub("_{2,}", "_", paste0("~/.Rhistory_", gsub(" |:", "_", date())))
  cat("Saving history to:", fname, "\n")
  cat("You were at:", getwd(), "\n")
  savehistory(file = fname)
  cat("Science rocks!\n")
}

# check group made strings are equivalent
comp_comb <- function(x, y, sep = 'n', v = TRUE){
  if(!all(c(length(x), length(y)))) return(FALSE)
  tvar <- sort(unlist(strsplit(as.character(x), sep)))
  tmp <- sort(unlist(strsplit(as.character(y), sep)))
  if(v){
    cat("x:", paste0(x, collapse = sep), '\n')
    cat("y:", paste0(y, collapse = sep), '\n')
  }
  if(length(tvar) != length(tmp)) return(FALSE)
  all(tvar == tmp)
}

# find equivalent string
find_eqs <- function(findx, vec, sep = 'n', v = FALSE){
  ii <- unlist(strsplit(as.character(findx), sep))
  if(sum(duplicated(ii), na.rm=T) > 0){
    cat('Warning - duplication:', ii[duplicated(ii)], 'in', findx,'\n')
  }
  tvar <- sapply(vec, function(x) comp_comb(findx, x, sep, v) ) # get logic vectors per group
  vec[tvar]
}

# Number of character occurences in a string
countChars <- function(char, string) {
    tvar <- gsub(char, "", string)
    return (nchar(string) - nchar(tvar))
}

# add needed spaces to reach a certain number of characters
addspaces <- function(x, m) paste0(x, paste0(rep(' ', m-nchar(x)), collapse = ''))

# add capitals
stringtrans <- function(
  x,
  sep = " ",
  npos = 1,
  up = TRUE
) { # just use str_to_sentence
  # x <- c("Roses are red, violets are blue", "My favourite colour is green")
  # str_replace_all(x, colours, col2hex)
  if(length(npos) == 1) npos <- rep(npos, 2)
  s <- strsplit(x, sep)
  sapply(s, function(y){
    s2trans <- substring(text = y, first = npos[1], last = npos[2])
    paste(casefold(x = s2trans, upper = up), substring(y, 2),
    sep="", collapse=" ")
  })
}

# print separated by commas
commas <- function(x, hn = 3){
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

# get only found
getfound <- function(x, sour, element = 'Elements', v = FALSE){
  if(v && sum(duplicated(x))) cat('Duplicated:', commas(x[duplicated(x)]), '\n')
  if(v) cat('Finding ', element, ': ', commas(x), '\n', sep = '')
  tmp <- !x %in% sour
  if(sum(tmp) && v) cat('Missing:', commas(x[tmp]), '\n')
  x <- x[!tmp]
  if(v) cat('Returning', length(x), casefold(element), '\n')
  return(x)
}

# limit length of a vector
vlimit <- function(x, mylim = Inf){
  if(length(x) > mylim) return(x[1:mylim]) else return(x)
}

# Erase all excluding a set of objects
cleanEnv <- cleanenv <- function(keep = NULL){
  if(is.null(keep)) keep <- 'no variables'
  myvars <- ls(pos = 1L)
  cat('Keeping:',paste0(keep, collapse = ', '),'\n')
  rm(list = myvars[!myvars %in% keep], pos = 1L)
  eval.parent(gc(), n = 1); gc()
  source('/mnt/BioHome/ciro/scripts/functions/handy_functions.R')
}

# repeat char/number in a data.frame with x rows
dfrep <- function(x, n) data.frame(Identity=rep(as.character(n), length(x)), row.names = x, stringsAsFactors = F)

# to get a df in the orther of a specified vector
getorder <- function(ids, df, cname){
  y <- unlist(sapply(ids, function(x) which(df[, cname] == x) ))
  nrows <- 1:nrow(df)
  if(length(y) != nrow(df)) y <- c(y, nrows[!nrows %in% y])
  return(y)
}
# now order a vector
group_order <- function(x, groups){
  groups <- as.character(groups)
  x <- as.character(x)
  c(groups[groups %in% x], x[!x %in% groups])
}

# transform string/vector/list to list for subset
getsubset_pattern <- translist <- function(pat){
  if(is.null(pat)) return(character(0))
  if(grepl(":|;", pat[[1]][1])){
    pat <- strsplit(unlist(strsplit(pat, ";")), ":")
    pat <- lapply(pat, function(x) unlist(strsplit(gsub(x, pattern = " ", replacement = ""), ",")) )
    return(pat)
  }
  if(sum(grepl("list", pat))) pat <- eval(expr = parse(text = pat))
  if(!is.list(pat)) pat <- list(pat)
  tvar <- any(sapply(pat, function(x) any(grepl('~', x)) && length(x) > 1 ))
  if(tvar) return(pat) # if it has the separator and is a vector > 1, return it
  pat <- lapply(pat, function(x){
    if(grepl('c\\(.*)', x[1])){
      eval(expr = parse(text = x))
    }else if(grepl('~', x[1])){
      unlist(strsplit(x, '~'))
    }else{ return(x) }
  })
  pat
}
pattern2list <- function(x){
  sapply(x, function(y){
    if(!grepl(":|;", y)) return(y)
    z <- gsub(";", "'), c\\('", gsub(":", "', '", y))
    paste0("\"list(c('", z, "'))\"")
  }, USE.NAMES = FALSE)
}

# subset data.frame with c('col','values') or a list(c('col1','values1'), c('col2','values2'))
getsubset <- function(
  x = NULL,
  df,
  op = 'and',
  v = FALSE
){
  df <- remove.factors(df)
  if(is.null(x)) return(rownames(df))
  x <- translist(x)
  # if(is.null(x)) x <- c(colnames(df)[1], unique(df[, 1]))
  if(length(x) == 1 && class(x) != 'list') x <- c(x, unique(df[, 1]))
  x <- x[!is.na(x)]
  if(class(x) != 'list') x <- list(x)
  if(v){
    maxchar <- max(c(14, max(nchar(sapply(x, head, 1)))))
    tvar <- paste(addspaces('Category', maxchar), ':=>\t Size \t Classes\n')
    cat(tvar)
    cat(paste0(rep("-", nchar(tvar) + 2), collapse = ''), '\n')
  }
  x <- x[sapply(x, function(s){
    if(!s[1] %in% colnames(df)){ warning("No '", s[1], "' in ", commas(colnames(df))); FALSE }
    TRUE
  })]
  x <- lapply(x, function(s){
    stmp <- sub("^\\-", "", s[-1])
    tvar <- stmp %in% unique(df[, s[1]])
    if(!all(tvar)) warning("No '", commas(stmp[!tvar]), "' in ", s[1])
    c(s[1], s[-1][tvar])
  })
  tmp <- sapply(x, function(s){
    tvar <- grepl('^\\-', s)
    if(sum(tvar)){ # in case we have a negative selection
      tmp <- unique(df[, s[1]])
      if(!is.null(levels(df[, s[1]]))) tmp <- as.character(unique(df[, s[1]]))
      tmp <- tmp[!tmp %in% sub('^\\-', '', s[tvar])] # reject negative ones
      s <- c(s[1], unique(c(tmp, s[-1][!tvar[-1]]))) # select targetted ones
      if(length(s) == 1) warning('no selection in ', s)
    }
    thisrows <- df[, s[1]] %in% s[-1]
    if(v){
      cat(addspaces(s[1], maxchar), ':=>\t', sum(thisrows), '\t', commas(s[-1]), '\n')
    }
    return(thisrows)
  })
  if(op == 'and') tmp <- rowSums(tmp) == length(x) # all x elements are in the row
  if(op == 'or') tmp <- apply(tmp, 1, any) # only some elements are in the row
  tvar <- rownames(df[tmp, , drop = FALSE])
  if(length(tvar) == 0) warning('In total no selection')
  if(v && length(x) > 1) cat(addspaces("Total elements", maxchar), ':=>\t', length(tvar), '\t', length(x), '\n')
  return(tvar)
}

summary_subset <- function(x){
  tvar <- sapply(x, function(x) paste0(paste0(x[-1], collapse = "AND"), "_from_", x[1]) )
  if(length(tvar) > 0) paste0(tvar, collapse = "_and_") else ""
}

# how to fit a number of objects in a 2D dimention
fitgrid <- function(x, nCol = 'x', v = FALSE){
  if(is.numeric(x) && length(x) == 1) x <- 1:x
  vecL <- length(x)
  if(v) cat('Length:', vecL, '\n')
  if (is.null(x = nCol)) {
    if(v) cat('Adjust for nCol:', nCol, '\n')
    nCol <- 2
    if (length(vecL) == 1) nCol <- 1
    if (length(vecL) > 6) nCol <- 3
    if (length(vecL) > 9) nCol <- 4
    nRow <- floor(x = length(vecL) / nCol - 1e-5) + 1
  }else{
    if(vecL %in% 1:2) nCol <- 1 else if(!is.numeric(nCol)) nCol <- round(sqrt(vecL) + .5)
    if(vecL == 1) nRow <- 1 else nRow <- round(vecL / nCol)
    if(nRow * nCol < vecL) nRow <- nRow + 1
    if(!sqrt(vecL) %% 1) nRow <- nCol <- sqrt(vecL) # if sqrt is integer
  }
  if(v) cat('nRow:', nRow, '\nnCol', nCol, '\n')
  return(c(nRow, nCol, ifelse(vecL <= 2,  5, 4)))
}

# get the legend from a plot
g_legend<-function(myplot){
  mylegend <- get_legend(myplot)
  if(file.exists('Rplots.pdf')) file.remove('Rplots.pdf')
  dev.off()
  return(mylegend)
}

# get data.frame with set operations
set_ops <- function(x, y, cname = NULL, addlen = FALSE, rename = NULL, tab = TRUE, v = F){
  if(!is.null(cname[1]) && sum(c(class(x), class(y)) == 'data.frame') == 2){
    if(length(cname) == 1) cname <- c(cname, cname)
    if(!cname[1] %in% colnames(x)){
      stop(cname[2], ' not in X column names'); return(0)
    }else if(cname[2] %in% colnames(y)){
      x <- x[, cname[1]]
      y <- y[, cname[2]]
    }else{ stop(cname[2], ' not in Y column names'); return(0) }
  }else if(!is.null(cname)){
    warning('No data.frame given: ', class(x), ', ', class(y), ' - ignoring ', cname)
  }
  if(class(x) != class(y)){
    warning('No matching classes: ', class(x), ', ', class(y))
    return(0)
  }
  if(sum(c(class(x), class(y)) == 'data.frame') == 2){
    warning('Data frames given - using vectors is recommended or give cname parameter')
  }
  setlist <- list(
    inter = intersect(x, y),
    union = union(x, y),
    x = setdiff(x, y),
    y = setdiff(y, x)
  )
  if(length(rename) == 2 && class(rename) == 'character' && !sum(rename %in% c('x', 'y'))){
    if(v) cat('Renaming\n')
    names(setlist) <- c('inter', 'union', rename[1], rename[2])
  }else{ rename <- c('x', 'y') }
  if(v){
    cat(paste0(rename[1], ': ', length(x)), '\n')
    cat(paste0(rename[2], ': ', length(y)), '\n')
    cat('Elements in operations:\n')
    print(sapply(setlist, length))
  }
  if(addlen) setlist <- sapply(setlist, function(x) c(length(x), x) )
  if(tab){
    data.frame(vlist2df(setlist), stringsAsFactors = F)
  }else{ return(setlist) }
}

# make a list of rownames out a column's categories or
# return categories for each group in another column
make_list <- function(x, colname = colnames(x)[1], col_objects = NULL, grouping = FALSE){
  x <- remove.factors(x)
  tvar <- lapply(unique(x[, colname]), function(y){
    rnames <- getsubset(c(colname, y), x, v = F)
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

# list of vectors to data.frame
vlist2df <- function(x, maxlen = NULL){
  if(!is.list(x)) stop('No list given')
  if(is.data.frame(x)) return(x)
  if(is.null(maxlen)) maxlen <- seq_len(max(sapply(x, length)))
  tvar <- data.frame(sapply(x, "[", i = maxlen), stringsAsFactors = F)
  colnames(tvar) <- names(x)
  return(tvar)
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

# Merge data.frames by column indexing
cbindList_i <- function(x, i = NULL, type = 'left'){
  suppressPackageStartupMessages(library("dplyr"))
  # if(!is.null(names(x))){
  #   cnames <- names(x)
  # }else{
  #   cnames <- paste0('X', 1:length(x))
  #   names(x) <- cnames
  # }
  if(is.null(i)){
    x <- list2evendf(x)
    x <- lapply(x, function(df) cbind(df, rancol487 = 1:nrow(df)))
    i <- "rancol487"
  }
  x <- switch(type,
    left = x %>% Reduce(function(dtf1, dtf2) left_join(dtf1, dtf2, by = i), .),
    full = x %>% Reduce(function(dtf1, dtf2) full_join(dtf1, dtf2, by = i), .),
    inner = x %>% Reduce(function(dtf1, dtf2) inner_join(dtf1, dtf2,by = i), .)
  )
  x <- x[, c(i, colnames(x)[colnames(x) != i])]
  # colnames(x) <- c(i, cnames)
  if('rancol487' %in% i) x <- x[, -which(colnames(x) == i)]
  x
}

# bind eliminating repeated new columns
cbind_repcol <- function(x, y, k = '1random302', v = FALSE){
  # new columns AND selected to keep
  tvar <- (!colnames(y) %in% colnames(x)) | colnames(y) %in% k
  tmp <- colnames(x) %in% colnames(y)[tvar] # exclude them from x if exist
  if(v){
    cat('Keeping in y:', colnames(y)[tvar], '\n')
    cat('Droping in x:', colnames(x)[tmp], '\n')
  }
  cbind(x[, !tmp, drop = F], y[rownames(x), tvar, drop = F])
}

# bind updating repeated new columns with second data.frame
# x <- metadata[1:(nrow(metadata) - 10002), 1:14]
# y <- metadata[1:(nrow(metadata) - 1020), 10:18]
# x$origlib <- factor(x$origlib)
# y$orig.peptide <- factor(y$orig.peptide)

joindf <- function(
  x,
  y,
  keep_from_y = "NULL123",
  type = c("left", "right", "full"),
  by_column = NULL
){
  type <- match.arg(type)
  if(is.null(by_column)){
    if(all(rownames(y) %in% rownames(x))){
      x$tmpcol123 <- rownames(x)
      y$tmpcol123 <- rownames(y)
      by_column = 'tmpcol123'
    }else{
      tvar <- apply(y, 2, function(mycolumn){ # find which column contains barcodes
        sum(mycolumn %in% rownames(x)) # then check if they're in the rownames
      });
      if(!any(tvar > 0)) stop('Non-matching names')
      rownames(y) <- y[, which.max(tvar[-length(tvar)])]
    }
  }
  yvars <- colnames(y)[(!colnames(y) %in% colnames(x)) | colnames(y) %in% keep_from_y]
  xvars <- colnames(x)[!colnames(x) %in% yvars] # exclude them from x if exist
  eval(expr = parse(text = paste0("functy <- dplyr::", type, "_join")))
  z <- functy(
    x = x[, xvars, drop = FALSE],
    y = y[, c(yvars, by_column), drop = FALSE],
    by = by_column
  )
  rownames(z) <- z[, by_column]
  z <- z[, -which(colnames(z) == by_column)]
  z
}
joindf <- function(
  x,
  y,
  keep_from_y = "NULL123",
  type = c("left", "right", "full", "none"),
  v = FALSE
){
  type <- match.arg(type)
  # str(x); str(y)
  if(type != "none"){
    x$tmpcol123 <- rownames(x)
    y$tmpcol123 <- rownames(y)
    yvars <- colnames(y)[(!colnames(y) %in% colnames(x)) | colnames(y) %in% keep_from_y]
    xvars <- colnames(x)[!colnames(x) %in% yvars] # exclude them from x if exist
    if(v) cat('Keeping in y:', yvars, '\n')
    if(v) cat('Droping in x:', ifelse(isTRUE(keep_from_y == 'NULL123'), "None", keep_from_y), '\n')
    if(v) cat("Using:", type, "\n")
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
    # classes1 <- c(sapply(x, class), sapply(y, class))
    # classes1 <- classes1[!duplicated(names(classes1))]
    z[rownames(x), colnames(x)] <- x
    icells <- intersect(rownames(x), rownames(y))
    z[icells, colnames(y)] <- y[icells, ]
  }
  # str(z)
  z
}

# get set operations from columns in two data.frames
df_overlap <- function(x, y, select_group = c('column', 'all'), compare = NULL,
                       per_group = FALSE, tab = FALSE, rename = NULL, v = F){
  if(v) cat('-------------------------\n')
  if(is.list(select_group)){
    if(v) cat('- It\'s a list!\n')
    void <- lapply(select_group, function(select){
      print(head(x))
      print(head(y))
      df_overlap(x, y, select_group = select, compare = compare,
        per_group = per_group, tab = tab, rename = rename, v = v)
    })
    if(tab){
      group_summary <- data.frame(data.table::rbindlist(
        lapply(void, function(x) data.frame(x[['group_summary']]))
      ), stringsAsFactors = T)
      myops <- lapply(void, function(x) x[['myops']])
      if(per_group){
        tvar <- make.names(unlist(sapply(select_group, function(x){ x[-1] })), unique = T)
      }else{
        tvar <- make.names(unlist(sapply(select_group, function(x){
          paste0(x[-1], collapse = '~')
        })), unique = T)
      }
      rownames(group_summary) <- tvar
    }else{ return(void) }
  }else{ # getting the groups
    tmp <- unique(c(x[, select_group[1]], y[, select_group[1]])); tmp <- tmp[!is.na(tmp)]
    if(per_group){
      if(v) cat('Groups in', select_group[1], '=>',commas(tmp), '\n')
      select_group <- c(select_group[1], tmp)
    }
    if(v) cat('Selection =>', commas(select_group[-1]), '\n')
  }
  if(per_group && !is.list(select_group)){
    if(v) cat('- Per group\n')
    myops <- list()
    for(gr in select_group[-1]){
      if(gr == 'all') gr <- tmp
      tvar <- df_overlap(x, y, select_group = c(select_group[1], gr), compare = compare,
                 per_group = FALSE, tab = tab, rename = rename, v = v)
      if(!exists('group_summary')){
        group_summary <- tvar[['group_summary']]
      }else{ group_summary <- rbind(group_summary, tvar[['group_summary']]) }
      if(length(gr) > 1) gr <- 'all'
      myops[[gr]] <- tvar[['myops']]
    }
  }else if(!is.list(select_group)){
    if(v) cat('---------- CALL ----------\n')
    if('all' %in% select_group[-1]) select_group <- unique(c(select_group[select_group!='all'], tmp))
    x <- x[getsubset(select_group, x, v), ]; tmp2 <- tmp
    y <- y[getsubset(select_group, y, v), ]
    myops <- set_ops(x, y, cname = compare, rename = rename, tab = tab, v = v)
    if(sum(tmp2 %in% select_group[-1], na.rm=T) == length(tmp2)) select_group <- c(1, 'all')
    if(is.data.frame(myops)){
      group_summary <- data.frame(t(sapply(myops, function(x){
        length(x[x != '' & !is.na(x)])
      })), row.names = paste0(select_group[-1], collapse = '~'))
    }else{ group_summary <- sapply(myops, length) }
  }
  if(tab){
    if(length(myops)>1 && is.data.frame(myops[[1]])){
      gaps <- TRUE
      myops <- sapply(myops, function(x) cbind(x, gap = rep('###', nrow(x))),
        simplify = FALSE, USE.NAMES = TRUE)
    }else{ gaps <- FALSE }
    myops <- cbindList(myops)
    if(gaps) myops <- myops[, -ncol(myops)]
  }
  return(list(group_summary = group_summary, myops = myops))
}


## Tests ##
# x <- data_frame(i = c("a","b","c"), j = 1:3)
# y <- data_frame(i = c("b","c","d"), k = 4:6)
# z <- data_frame(i = c("c","d","a"), l = 7:9)
# z <- data_frame(i = c("c","d","a", "n", "4"), l = 7:11)
# cbindList_i(list(x, y, z), type = 'full')

# find the r^2
lm_eqn <- function(df, x, y){
  formla <- as.formula(y ~ x)
  m <- lm(formla, df);
  eq <- substitute(italic(y) == a + b %.% italic(x)*", "~~italic(r)^2~"="~r2,
       list(a = format(coef(m)[1], digits = 2),
            b = format(coef(m)[2], digits = 2),
           r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

#### for groups intersections #### ---------------------------------------------
interset <- function (x){
  # Multiple set version of intersect
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    intersect(x[[1]], x[[2]])
  } else if (length(x) > 2){
    intersect(x[[1]], interset(x[-1]))
  }
}
unset <- function (x) {
  # Multiple set version of union
  # x is a list
  if (length(x) == 1) {
    unlist(x)
  } else if (length(x) == 2) {
    union(x[[1]], x[[2]])
  } else if (length(x) > 2) {
    union(x[[1]], unset(x[-1]))
  }
}
diffset <- function (x, y) {
  # Remove the union of the y's from the common x's.
  # x and y are lists of characters.
  xx <- interset(x)
  yy <- unset(y)
  setdiff(xx, yy)
}
overlap_calc <- function(
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
  combinations <- getcombn(grps = grps, combinations = combinations, sep = sep, sharedmax = sharedmax, v = v)
  if(v) cat('Calculating overlaping elements\n')
  tvar <- lapply(combinations, function(i){ # Get overlapping elements
            diffset( groups_list[i], groups_list[setdiff(names(groups_list), i)] )
          })
  names(tvar) <- names(combinations)
  return(tvar)
}
#### ------------------------ #### ---------------------------------------------

# get a list of group combinations
getcombn <- function(
  grps,
  combinations = NULL,
  sep = NULL,
  logchar = "<=",
  sharedmax = Inf,
  nameorder = NULL,
  simplify = FALSE,
  v = FALSE
){
  if(is.null(sep)) sep <- getsepchar(grps)
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

# when you have a vector of filters you want to place separated by
further_filters <- function(
  x,
  cnames = NULL,
  th = "",
  op = "&",
  v = FALSE
){
  if(is.list(x)){
    if(v) cat("It is a list of", length(x), "elements\n")
    # if it is a list, it will take the name of each element as column
    # and the content as the condition
    x <- paste(sapply(1:length(x), function(cond){
      numlim <- x[cond]
      newexpr <- paste0("x[, '", names(x[cond]), "']")
      if(!grepl(">|<|=", numlim)){
        numlim <- paste0(">", th, " ", numlim)
        print(newexpr)
        eval(expr = parse(text = paste0("tvar <- any(-", newexpr, ">0)")))
        if(tvar) newexpr <- paste0("abs(", newexpr, ")")
      }
      paste0(newexpr, numlim)
    }), collapse = paste0(" ", op, " "))
  }else if(!is.null(cnames)){ # if it is a expression
    if(v) cat("Digesting:", x, "\n")
    xcolumns <- paste0("(", paste0(cnames, collapse = "|"), ")")
    x <- gsub(xcolumns, paste0("x[, '", "\\1", "']"), x)
  }else{
    warning("Nothing to filter\n")
  }
  x
}

fix_gene_excel_err <- function(x){
  if(any(grepl("1-Sep", x))) warning("SEP-1 or SEPT1 conflict for 1-Sep")
  y <- gsub("([0-9]{1,})\\-Sep", "SEPT\\1", x)
  y <- gsub("([0-9]{1,})\\-Mar", "MARCH\\1", y)
  y <- gsub("([0-9]{1,})\\-Nov", "NOV\\1", y)
  y <- gsub("December ([0-9]{1,})", "DEC\\1", y)
  y <- gsub("October ([0-9]{1,})", "OCT\\1", y)
  y <- gsub("2006/09/02", "SEPT2", y)
  return(y)
}

getDEGenes <- function(
  x,
  pv = 0.2,
  fc = 0.0,
  upreg = NA,
  catch = NULL,
  pvtype = 'padj',
  lfc.type = 'log2FoldChange',
  th = "=", # add equal to > and <
  gene_name = NULL,
  further = NULL, # list or expression of columns which you want to filter with
  v = FALSE
){
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if(sum(gene_name %in% colnames(x))){
    rownames(x) <- x[, gene_name]
  }
  if(!is.null(catch)){
    pv <- max(x[catch, pvtype], na.rm = TRUE)
    # fc <- max(abs(x[catch, lfc.type]))
  }
  if(v) cat(pvtype,': ', pv, '\n', sep = "")
  if(v) cat(lfc.type,': ', fc, '\n', sep = "")

  assiggy <- "g <- rownames(x[which(x[, pvtype] <"
  mid <- "x[, lfc.type]"
  suf <- "fc"
  if(is.na(upreg)){
    if(v) cat('Both\n')
    express <- paste0(assiggy, th, "pv & abs(", mid, ") >", th, suf)
  }else if(isTRUE(upreg)){
    if(v) cat('Upregulated\n')
    express <- paste0(assiggy, th, "pv & ", mid, ">", th, suf)
  }else{
    if(v) cat('Downregulated\n')
    express <- paste0(assiggy, th, "pv & ", mid, "<", th, " -", suf)
  }
  if(!is.null(further)){ # Extra condition for which filter
    tvar <- further_filters(x = further, th = th, cnames = colnames(x), v = v)
    if(tvar != further) express <- paste(express, "&", tvar)
  }
  express <- paste0(express, "), ])")
  if(v) cat("------------\n", express, "\n------------\n")
  eval(expr = parse(text = express))
  n.gs <- length(g) # Get number retrieved genes
  if(v) cat('Retrieved genes:', commas(g), '\n')
  if(n.gs == 0 || is.null(g) || is.na(g)){
    return(character(0))
  }else{
    return(g)
  }
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
v2cols <- function(
  select,
  sour = NULL,
  cname = NULL,
  uniq = TRUE,
  fw = 'gg', # force wgcna colours
  myseed = 27,
  in_hex = FALSE,
  v = FALSE
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
    if(v) warning('Names: ', commas(unkno),' not in object')
    if(length(select)){
      select <- v2cols(select, sour = sour, cname = cname, uniq = uniq, fw = fw, v = v)
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
  if(v) cat('Colours from ')
  if((is.null(sour) && length(dim(select)) > 1) || isTRUE(fw)){
    if(v) cat('labels2colors\n')
    thiscols <- WGCNA::labels2colors(select)
    tmp <- grepl("white", thiscols)
    set.seed(myseed); thiscols[tmp] <- sample(colorsnogray, sum(tmp))
  }else if(is.null(sour)){
    if(v) cat('ggplot\n')
    thiscols <- gg_color_hue(length(select))
  }else if(class(sour) == 'data.frame'){
    if(v) cat('data.frame\n')
    thiscols <- sour[select, cname]
  }else if(class(sour) == 'character'){
    if(v) cat('vector\n')
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

# get cells that express specific genes
gen_sharing <- function(
    mydata,
    genes,
    neg = TRUE,
    thiscells = NULL,
    df = FALSE,
    v = FALSE,
    ...
  ){
  if(is.null(thiscells)){
    if(v) cat('Using all cells\n')
    thiscells <- colnames(mydata)
  }
  if(v) cat('Genes:\n')
  tmp <- sapply(genes, function(g){
    if(v) cat(' -', g,'\n')
    clls <- mydata[g, thiscells] > 0
    clls <- names(clls[clls])
    return(clls)
  })
  if(v) cat('Calculating overlaps\n')
  tmp <- overlap_calc(tmp, ...)
  if(isTRUE(neg)){
    if(v) cat('Adding negative population\n')
    tmp$NEG <- setdiff(thiscells, as.vector(unlist(tmp)))
  }
  if(length(genes) > 1){
    if(isTRUE(df)){
      if(v) cat('As table\n')
      tmp <- vlist2df(tmp)
    }
  }
  return(tmp)
}

# proportion of cells expressing genes
gen_prop <- function(
    mydata,
    thiscells,
    cname = 'res.0.6',
    norm_type = 'e',
    v = FALSE
  ){
  if(norm_type == 'e' && v) cat('Normalising to the number of expressing cells\n')
  if(norm_type == 'c' && v) cat('Normalising to the number of cells in cluster\n')
  print(names(thiscells))
  props <- cbindList_i(
    lapply(names(thiscells), function(myname){
      if(v) cat(' - ',myname,'\n')
      thismeta <- mydata[rownames(mydata) %in% thiscells[[myname]], ]
      xtreg <- table(thismeta[, cname])
      all <- table(mydata[, cname])
      # if a group is missing
      xtreg[names(all)[!names(all) %in% names(xtreg)]] <- 0
      xtreg <- xtreg[sort(names(xtreg))]
      all <- all[sort(names(all))]
      if(v) cat('.')
      if(norm_type == 'e'){ # might not make sense
        tvar <- table(rep(names(xtreg), sum(xtreg))) # all expressing cells
        mydata <- data.table::data.table((xtreg / tvar) * 100)
      }else if(norm_type == 'c'){ # cluster size
        mydata <- data.table::data.table((xtreg / all) * 100)
      }else{
        mydata <- data.table::data.table(xtreg)
      }
      colnames(mydata) <- c(cname, myname)
      return(mydata)
    }), i = cname)
  if(v) cat('Getting proportions\n')
  # return(props)
  void <- get_props(props = props, group = "Gene", cluster = cname,
    norm_type = 'n', return_plot = TRUE, v = T)
  return(list(props = props, plots = void))
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

# set directory with names
setdir <- function(x, pats = NULL, v = FALSE){
  pats <- unlist(strsplit(pats, "/"))
  pats <- unique(pats)
  x <- dircheck(x)
  if(v) cat("Setting folder with root:", x, "\n")
  if(v) cat("Patterns:", commas(pats), "\n")
  n <- length(pats)
  for(i in 1:n){
    tmp <- list.files(x)
    tmp <- tmp[tmp %in% pats]
    if(length(tmp) == 1){
      if(v) cat(rep(' ', i - 1), "|--", tmp, "\n")
      x <- paste0(x, tmp, '/')
    }
  }
  return(x)
}

set_dir_pats <- function(x, pj = '10XData', v = FALSE){
  mypats <- c('outs', 'filtered_feature_bc_matrix', 'filtered_gene_bc_matrices_mex',
    pj, 'hg19', 'filtered_gene_bc_matrices_h5.h5', 'filtered_feature_bc_matrix.h5')[-c(6:7)]
  # if(!dir.exists(x)) return(x)
  # list.files(x, recursive = T, pattern = "filtered_", full.name = T) # takes too long
  setdir(x, pats = mypats, v = v)
}

# find a file upstream in a path
findfile <- function(name, path = getwd(), read = FALSE, nup = 4, v = FALSE, ...){
  path <- dircheck(path)
  tmp <- ''
  for(i in 1:nup){
    myfile <- paste0(path, tmp, name)
    tmyfile <- system(paste('ls', myfile), intern = T)
    if(length(tmyfile) > 0) myfile <- tmyfile
    if(file.exists(myfile)){
      if(v) cat(basename(myfile), 'found!\n')
      if(read) return(readfile(myfile, v = v, ...)) else return(myfile)
    }
    tmp <- paste0(tmp, '../'); # print(myfile)
  }
  if(v) cat(name, 'not found\n')
  # if(!read) stop(myfile, ' not found')
  return(data.frame(lol = 'void', stringsAsFactors = F))
}

# get arguments
matchArgs <- function(x, func){
  x <- x[names(x) %in% methods::formalArgs(func)]
  if(!length(x)) return(NULL)
  x[!duplicated(names(x))]
}

# detect file format and read it
readfile <- function(
  myfile,
  usedt = FALSE,
  v = FALSE,
  futil = list(stringsAsFactors = FALSE, check.names = FALSE),
  ...
){
  futil = c(list(...), futil); futil <- futil[!duplicated(names(futil))]
  if(length(futil) && v) cat('Extra parameters:', commas(names(futil)), '\n')
  if(isTRUE(usedt) && grepl("tsv$|csv$|txt$", myfile)){
    file_content <- data.frame(data.table::fread(myfile), stringsAsFactors = FALSE, check.names = FALSE)
    rownames(file_content) <- file_content[, 1]; file_content <- file_content[, -1, drop = FALSE]
  }else if(grepl(paste0(c('rda', 'rta', 'rdata', 'robj'), '$', collapse = "|"), myfile, ignore.case = TRUE)){
    if(v) cat(' - From RData file\n'); #futil <- matchArgs(futil, theObjectSavedIn)
    file_content <- theObjectSavedIn(myfile)
  }else if(grepl('csv$', myfile, ignore.case = TRUE)){
    if(v) cat(' - From CSV file\n')
    file_content <- read.csv(myfile, ...)
  }else if(grepl('h5$|h5$', myfile, ignore.case = TRUE)){
    if(v) cat(' - 10X H5\n')
    file_content <- Seurat::Read10X_h5(myfile)
  }else if(grepl('txt$|tsv$', myfile, ignore.case = TRUE)){
    if(v) cat(' - From TXT file\n')
    file_content <- read.table(myfile, ...)
  }else if(grepl('rds$', myfile, ignore.case = TRUE)){
    if(v) cat(' - From RDS file\n'); futil <- matchArgs(futil, readRDS)
    file_content <- readRDS(myfile, ...)
  }else if(grepl('loom$', myfile, ignore.case = TRUE)){
    file_content <- connect(filename = myfile, mode = "r")
  }else{ cat('Unknown format:', basename(myfile), '\n') }
  file_content
}
############ MORE COMPLEX ONES ############

# find power for set of values in a matrix
findPower <- function(x, limit = 0.9, type = 'min', do_plot = F, outplot = '', v = F){
  powers = c(c(1:10), seq(from = 12, to=20, by=2))
  sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = v)
  if(outplot!='') outplot <- dircheck(outplot) # Fix for directory
  if(do_plot){
    if(v) cat('Plotting')
    pdf(paste0(outplot,'soft_threshold_power.pdf'), width = 8, height = 5 )
    par(mfrow = c(1,2))
    plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
         xlab='Soft Threshold (power)', ylab='Scale Free Topology Model Fit,signed R^2',type='n',
        main = paste('Scale independence'));
    text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
        labels=powers, cex=0.9, col='red');
    abline(h=0.90,col='red')
    plot(sft$fitIndices[,1], sft$fitIndices[,5],
        xlab='Soft Threshold (power)', ylab='Mean Connectivity', type='n',
        main = paste('Mean connectivity'))
    text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9, col='red')
    dev.off()
  }

  corrs <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
  if(v) cat('Correlations:\n'); print(corrs)
  if(!type %in% c('min', 'max')) cat('MIN after cutoff chosen\n')
  if(type == 'min') tvar <- min(which(corrs > limit))
  if(type == 'max') tvar <- max(which(corrs > limit))
  if(tvar == Inf){ tvar <- min(which(corrs == max(corrs))); cat('Threshold not passed, taking maximum correlation:',tvar,'\n') }
  cat('Tentative soft power:',tvar,'\n')
  return(tvar)
}

batch_correction <- function(
  mat,
  coldata = NULL,
  batch = "batch",
  covar = NULL,
  type = 'combat',
  cut.mean = 0,
  transf = 'log10',
  v = FALSE
){
  mat2deseq <- function(x){
    DESeqTransform(SummarizedExperiment(x))
  }
  if(class(mat) == 'DESeqDataSet'){
    coldata <- colData(mat)
    mat <- counts(mat)
  }else if(!is.null(coldata)){
    coldata <- DataFrame(coldata)
  }else{
    stop("You need to give 'coldata'")
  }
  mat <- as.matrix(mat)

  coldata$batch <- factor(coldata[, batch])
  if(!is.null(covar)){
    if(v) cat('Using covariate\n')
    coldata$covar <- factor(coldata[, covar])
  }

  # transformation
  idx <- rowMeans(mat) > cut.mean # filtering
  if(v) cat(sum(!idx),' genes cut by average...\n')
  mat <- mat[idx, ]

  # transformation
  if(grepl('log', transf)){
    if(v) cat('Shifted log10 transformation\n')
    bs <- as.numeric(sub('log', '', transf))
    if(is.na(bs)) bs <- exp(1)
    print(min(mat))
    print(max(mat))
    vsd <- mat2deseq(log(mat + 1, bs))
    cat('cries\n')
  }else if(transf == 'vst'){
    if(v) cat('Variance Stabilizing transformation\n')
    vsd <- mat2deseq(varianceStabilizingTransformation(mat))
  }else{
    if(v) cat('No transformation performed\n')
    vsd <- mat2deseq(mat)
  }
  print(ncol(vsd))
  colData(vsd) <- coldata
  colnames(vsd) <- 1:ncol(vsd)
  # print(vsd)

  if(v) cat('Correction: ')
  if(type == 'limma'){
    if(v) cat('Limma\n')
    if(!is.null(covar)) covar <- as.numeric(vsd$covar)
    assay(vsd) <- limma::removeBatchEffect(x = assay(vsd), batch =  vsd$batch, covariates = covar)
  }else if(type == 'combat'){
    if(v) cat('ComBat\n')
    if(!is.null(covar)){
      modcombat <- model.matrix(as.formula(paste0( '~ as.factor(',covar,')')), colData(vsd))
    }else{
      modcombat <- model.matrix(~1, data = colData(vsd))
    }
    if(v){
      assay(vsd) <- sva::ComBat(dat = assay(vsd), batch = vsd$batch, mod = modcombat)
    }else{
      assay(vsd) <- suppressMessages(ComBat(dat = assay(vsd), batch = vsd$batch, mod = modcombat))
    }
  }else{
    warning("No batch correction performed, select 'limma' or 'combat'")
  }
  if(v) cat('Finished\n\n')
  return(vsd)
}

see_batch_effect <- function(
  mat,
  coldata = NULL,
  type = 'limma',
  batch = "batch",
  covar = NULL,
  cut.mean = 0,
  transf = 'vst',
  return_plots = TRUE,
  plot_type = 'suas',
  grcols = NULL,
  v = FALSE,
  pcs = 3,
  genes = NULL
){
  tmp <- list() # plot
  intg <- batch
  if(!is.null(covar)){
    intg <- c(batch, covar)
  }
  if(pcs > 3) pcs = 3

  if(v) cat('------ transformed data\n')
  vsd <- batch_correction(mat, coldata = coldata, batch = batch, covar = covar,
              type = 'type', cut.mean = cut.mean, transf = transf, v = v)
  if(is.null(genes)) genes <- rownames(assay(vsd)) else cat('using',length(genes),'genes\n')
  if(plot_type == 'suas'){
    if(is.null(covar)) tvar <- NULL else tvar <- batch
    tvar <- pcaplot(t(assay(vsd)[genes, ]), main.f = "Batch effect",
                  coldata = colData(vsd), batchCol = batch, cd = covar,
                        ptsize = 3, printplot = FALSE, logp = FALSE,
                        cols = c(covar, batch), shapes = tvar,
                        grcols = grcols, hidleg = TRUE)
    tmp <- tvar[1:pcs]
  }else{
    tmp[[1]] <- plotPCa(vsd, intgroup = intg, titl = "Batch effect", legend = TRUE, ntop = 5000)
  }
  if(v) cat('------ batch correction\n')
  vsd <- batch_correction(assay(vsd), coldata = colData(vsd), batch = batch, cut.mean = cut.mean,
                covar = covar, type = type, transf = 'none', v = F)
  if(plot_type == 'suas'){
    if(is.null(covar)) tvar <- NULL else tvar <- batch
    tvar <- pcaplot(t(assay(vsd)[genes, ]), main.f = paste(type, "removed batch effect"),
                  coldata = colData(vsd), batchCol = batch, cd = covar,
                        ptsize = 3, printplot = FALSE, logp = FALSE,
                        cols = c(covar, batch), shapes = tvar,
                        grcols = grcols)
    tmp <- c(tmp, tvar[1:pcs])
  }else{
    tmp[[2]] <- plotPCa(vsd, intgroup = intg, legend = TRUE, titl = paste(type, "removed batch effect"))
  }
  if(return_plots) return(tmp)
  tvar <- fitgrid(tmp)
  do.call(grid.arrange, c(tmp, ncol = tvar[1], nrow = tvar[2]))
}

# creates a list of names based on a list's of DF's colum names
getnames <- function(
  stats_names,
  padjs = 'padj',
  fcs = 'log2FoldChange'
){
  if(length(padjs) < length(stats_names)){
    padjs <- as.list(rep(padjs, length(stats_names)))
    names(padjs) <- stats_names
  }
  if(length(fcs) < length(stats_names)){
    fcs <- as.list(rep(fcs, length(stats_names)))
    names(fcs) <- stats_names
  }
  return(list(padjs = padjs, fcs = fcs))
}

# get a separator character
getsepchar <- function(
  groups, # group names
  opts = c('n', '~', 'u', '_', '-', 'xx')
){
  if(length(opts) == 0) opts = c('n', '~', 'u', '_', '-', 'xx')
  if(all(sapply(groups, function(x) x == casefold(x) ))){
    tvar <- sapply(opts, function(x) any(unlist(strsplit(x, split = "")) %in% letters) )
    opts <- c(opts[!tvar], opts[tvar])
  }
  tvar <- sapply(opts, function(x) !grepl(x, groups) )
  newsep <- head(opts[colSums(tvar) == length(groups)], 1)
  if(length(newsep) == 0){
    newsep <- head(opts, 1)
    warning("Ambiguity might be found with '", newsep, "' as separator\n", sep = "")
  }
  return(newsep)
}

# combine stats and expression
cbind_dea <- function(
  comps_stats,
  expr_mat = NULL,
  fpvtype = 'padj',
  ffctype = 'log2FoldChange',
  v = FALSE
){
  degslist <- unique(unlist(lapply(comps_stats, rownames)))
  tvar <- getnames(names(comps_stats), fpvtype, ffctype)
  fpvtype <- tvar[[1]]
  ffctype <- tvar[[2]]
  if(v) message('Merging ', length(degslist), ' genes')
  suppressPackageStartupMessages(library("doParallel"))
  registerDoParallel(cores=parallel::detectCores()) # For stat table for comparisons
  comprs.stats <- foreach(cp = names(comps_stats), .combine = cbind) %dopar% {
    tmp <- comps_stats[[cp]]
    if(isTRUE('gene_name' %in% colnames(tmp))) rownames(tmp) <- sub("'","",tmp$gene_name)
    tvar <- !degslist %in% rownames(tmp)
    if(sum(tvar)) cat('Missing genes in', cp, ':', sum(tvar), '\n')
    comptab <- data.frame(tmp[degslist[!tvar], c(ffctype[[cp]], fpvtype[[cp]])])
    tmp <- mat_names(degslist[tvar], c(ffctype[[cp]], fpvtype[[cp]])) # adding missing genes
    comptab <- rbind(comptab, tmp)
    comptab <- comptab[degslist, ] # recovering order
    colnames(comptab) <- c(paste0(ffctype[[cp]],'(',cp,')'),paste0(fpvtype[[cp]],'(',cp,')'))
    return(comptab)
  }
  rownames(comprs.stats) <- degslist
  # tvar <- lapply(comps_stats, function(x) as.data.frame(t(as.matrix(x))) )
  # comprs.stats_dt <- t(data.table::rbindlist(tvar, fill = T, use.names = TRUE))

  # Adding extra stats of the stats: mean, max and min of p-values and fold changes
  degslist[!degslist %in% rownames(comprs.stats)]
  rownames(comprs.stats)[!rownames(comprs.stats) %in% degslist]
  tvar <- sort(unique(sub("^(.*)\\(.*", "\\1", colnames(comprs.stats))))
  if(grepl('adj', tvar[1])) tvar <- sort(tvar, decreasing = T) # is nomenclature is putting padj first
  # adding boundaries
  if(v) cat('boundaries\n')
  boundaries <- cbindList(lapply(tvar, function(x){
    mat <- as.matrix(comprs.stats[, grepl(x, colnames(comprs.stats))])
    tmp <- data.frame(
      min = matrixStats::rowMins(mat, na.rm = TRUE),
      mean = rowMeans(mat, na.rm = TRUE),
      max = matrixStats::rowMaxs(mat, na.rm = TRUE))
    colnames(tmp) <- paste0(x, "_", colnames(tmp))
    return(tmp)
  }))
  head(boundaries)
  tvar <- unlist(sapply(tvar, function(x) which(grepl(x, colnames(comprs.stats))), simplify = F))
  comprs.stats <- cbind(comprs.stats[, tvar], boundaries)
  head(comprs.stats); tail(comprs.stats)
  comprs.stats <- cbind(gene_name = paste0("'", degslist), comprs.stats)

  # Using metap methods for combination of p-values

  if(is.null(expr_mat)) return(comprs.stats)
  if(v) message('Merging stats with expression a.k.a useless table')
  useless_table <- cbind(
    comprs.stats,
    expr_mat[rownames(comprs.stats), ]
  )
  return(useless_table)
}

# subset all strings in vector #
cutstrings <- function(x, ssize = 4, v = FALSE){
  x <- as.character(x)
  if(v) cat(length(x), 'elements\n')
  if(v) cat('Mean size', mean(nchar(x)), '\n')
  if(v) cat('Down to', ssize*2, '\n')
  ifelse(nchar(x) > ssize*2, paste0(substring(x, 1, ssize), '-', substring(x, nchar(x)-ssize, nchar(x))), x)
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
# found in column
colfound <- function(x, tab) sum(sapply(x, function(y) grepl(paste0('^', y, '$'), colnames(tab))))

# Eliminating ensmbl number
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

# find genes even if they have ensembl ids
findgenes <- function(x, y){
  if(all(x %in% y)) return(x[x %in% y])
  tvar <- grepl(paste(paste0("^", parse_ens_name(x), "$"), collapse = "|"), parse_ens_name(y), ignore.case = TRUE)
  y[tvar]
}

found_partial <- function(x, y){
  # tvar <- paste0(x, collapse = "|")
  z <- unlist(sapply(x, function(tvar) y[grepl(tvar, y)] ), use.names = FALSE)
  return(z)
}

cosmy <- function(x, patties = "^rps|^rpl|^mt-|rp[0-9]{1,}\\-"){
  tvar <- grepl(patties, casefold(parse_ens_name(x)))
  if(sum(tvar)) x <- x[!tvar]
  return(x)
}

# set rownames from column
setrows <- function(tab, x = 'gene_name', v = FALSE){
  x <- getfound(x, colnames(tab))
  if(sum(x %in% colnames(tab))){
    if(v) cat('Setting', commas(x), 'as rows\n')
    newrows <- do.call('paste', c(tab[x], sep=''))
    newrows <- sub("'", "", newrows)
    if(v) cat(commas(newrows), '\n')
    rownames(tab) <- newrows
  }
  return(tab)
}

# Add totals to a numeric table
calc_tots <- function(x){
  x <- rbind(x, total = colSums(x))
  cbind(x, total = rowSums(x)) # adding totals
}

table_pct <- function(df, cnames, total = NULL, subsamp = NULL){
  fdf <- table(df[, cnames])
  if(!is.null(total)){
    tdf <- data.table::melt(table(df[, c(cnames[1], total)]))
    tdf <- tdf[tdf[, 3] > 0, ]
    tvar <- stats::aggregate(tdf[, 3], by = list(Category = tdf[, 2]), FUN = sum)
    rownames(tvar) <- tvar[, 1]
    tdf[, 3] <- tvar[tdf[, 2], 2]
    rownames(tdf) <- tdf[, 1]
    fdf <- rbind(fdf, total = colSums(fdf))
    fdf <- cbind(fdf, total = tdf[rownames(fdf), 3])
    fdf['total', 'total'] <- sum(fdf['total', ], na.rm = T)
    colnames(fdf)[3] <- paste0(colnames(fdf)[3], "_", sub("orig\\.", "", total))
  }else{
    fdf <- calc_tots(fdf)
  }
  if(!is.null(subsamp)){
    fdf <- fdf[c(rownames(table(df[getsubset(subsamp, df, v = TRUE), cnames])), 'total'), ]
    rownames(fdf)[nrow(fdf)] <- paste0(rownames(fdf)[nrow(fdf)], "_", sub("orig\\.", "", cnames[2]))
  }
  return(fdf)
}

# To get the significative PCs given a threshold for JackStraw data
choosePCs <- function (object, PCs = 1:5, mythres = 0.05, score.thresh = 1e-05){
  vrs <- as.numeric(sub("(^...).*", "\\1", object@version))
  if(vrs < 3){
    pAll <- object@dr$pca@jackstraw@emperical.p.value
  }else{
    pAll <- object@reductions$pca@jackstraw@empirical.p.values
  }
  pAll <- as.data.frame(pAll[, PCs, drop = FALSE]); pAll$Contig <- rownames(x = pAll)
  pAll.l <- data.table::melt(data = pAll, id.vars = "Contig")
  colnames(x = pAll.l) <- c("Contig", "PC", "Value")
  qq.df <- NULL; score.df <- NULL
  for (i in PCs) {
      pc.score <- suppressWarnings(prop.test(x = c(length(x = which(x = pAll[,
          i] <= score.thresh)), floor(x = nrow(x = pAll) *
          score.thresh)), n = c(nrow(pAll), nrow(pAll)))$p.val)
      if (length(x = which(x = pAll[, i] <= score.thresh)) ==
          0) {
          pc.score <- 1
      }
      if (is.null(x = score.df)) {
          score.df <- data.frame(PC = paste0("PC", i), Score = pc.score)
      }
      else {
          score.df <- rbind(score.df, data.frame(PC = paste0("PC",
              i), Score = pc.score))
      }
  }; cat('JackStraw scoring:', commas(round(score.df$Score, 5)), '\n')
  tvar <- !score.df$Score <= mythres
  tvar <- ifelse(sum(tvar) != 0, min( which(tvar == TRUE) ) - 1, length(tvar))
  cat('Chosen:', tvar, '\n')
  tvar
}

# to decide what pattern's case identifies some strings
choose_case <- function(x, vec){
  tvar <- c(sum(grepl(x, vec)),
    sum(grepl(casefold(x), vec)),
    sum(grepl(casefold(x, upper = TRUE), vec))
  )
  tvar <- min(which(tvar == max(tvar)))
  c(x, casefold(x), casefold(x, upper = TRUE))[tvar]
}

# to divide data in quantiles
quantile_breaks <- function(xs, n = 10) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
over_quantile <- function(x, nn = 11, q = 90){
  cutoff <- quantile_breaks(x, n = nn) # over the 90% percentile
  tvar <- round(as.numeric(sub("%", "", names(cutoff))), 3)
  cuto <- cutoff[min(which(tvar >= q))]
}

# calculate module scores ranking genes and sum those ranks for each cell
# then regress the number of genes
genes2sign <- function(
  object,
  genes_use = NULL,
  names = 'signature',
  reg_var = "nFeature_RNA",
  v = FALSE
){
  if(is.null(genes_use)){
    if(v) cat('Random selection of genes\n'); set.seed(27)
    genes_use <- lapply(names, function(x) sample(rownames(object@assays$RNA@data), 20) )
  }
  names(genes_use) <- names
  for(name in names(genes_use)){
    if(v) cat('Name:', name, '\n')
    gg <- getfound(genes_use[[name]], rownames(object@assays$RNA@data), element = 'Genes', v = v)
    dat <- object@assays$RNA@data[gg, rownames(object@meta.data)]
    if(v) cat('Ranking\n')
    rank.mat <- apply(dat, 1, rank, ties.method = "average")
    if(v) cat('Sum\n')
    cell_score <- apply(rank.mat, 1, sum)
    object@meta.data[, name] <- cell_score
    form <- paste0(name, "~", paste0(reg_var, collapse = "+"))
    if(v) cat('Residuals', form,'\n')
    form <- formula(form)
    object@meta.data[, name] <- scale(summary(lm(form, object@meta.data))$residuals)[, 1]
    if(v) print(summary(object@meta.data[, name]))
  }
  return(object)
}

# calculate module scores with z-score
genes2sign <- function(
  object,
  genes_use = NULL,
  names = 'signature',
  reg_var = "nFeature_RNA",
  v = FALSE
){
  if(is.null(genes_use)){
    if(v) cat('Random selection of genes\n'); set.seed(27)
    genes_use <- lapply(names, function(x) sample(rownames(object@assays$RNA@data), 20) )
  }
  names(genes_use) <- names
  for(name in names(genes_use)){
    if(v) cat('Name:', name, '\n')
    str(genes_use)
    thesegenes <- if(is.list(genes_use[[name]])) genes_use[[name]][[1]] else genes_use[[name]]
    gg <- getfound(thesegenes, rownames(object@assays$RNA@data), element = 'Genes', v = v)
    dat <- object@assays$RNA@data[gg, rownames(object@meta.data)]
    if(is.numeric(genes_use[[name]][[2]])){
      if(v) cat("Transforming\n")
      dat <- sweep(x = dat, MARGIN = 1, genes_use[[name]][[2]], '*')
    }
    if(v) cat('Scaling\n')
    scadat <- t(scale(t(dat)));
    scadat[scadat < (-2)] <- -2
    scadat[scadat > (2)] <- 2
    if(v) cat('Sum\n')
    cell_score <- apply(X = scadat, MARGIN = 2, FUN = sum)
    object@meta.data[, name] <- cell_score
    if(v) print(summary(object@meta.data[, name]))
  }
  return(object)
}

# create factors with levels in the right order
factormix <- function(x){
  if(!is.character(x)) return(x)
  y <- as.character(x)
  factor(y, levels = gtools::mixedsort(unique(y)))
}

# midpoints
mids <- function(x) c(x[1]/2, x[-length(x)] + diff(x)/2)

# find elbow
get_elbow <- function(
  x,
  y = NULL,
  threshold = .90,
  decide = FALSE,
  v = FALSE
) {
  if(is.null(y)){ y <- x; x <- 1:(length(y)) }
  # ap <- approx(x, y, n=1000, yleft=min(y), yright=max(y))
  # x <- ap$x; y <- ap$y
  mindexes <- sapply(threshold, function(z) {
    d1 <- diff(y) / diff(x) # first derivative
    d2 <- diff(d1) / diff(x[-1]) # second derivative
    thh <- abs(quantile(d2, probs = z))
    indices <- which(abs(d2) >= thh)
    if(v){
      cat('SDEV thh:', thh, '\n')
      cat("D'':", commas(round(d2[indices], 3)), '\n')
      cat('SDEVs:', commas(round(y[indices], 3)), '\n')
      cat('Chosen:', x[max(indices)], '->', y[max(indices)], '\n')
    }
    max(indices)
  })
  if(length(threshold) == 1)
    return(list(elbows = mindexes, ptable = data.frame(y = max(y), x = mindexes + 0.25, N = 1), elbow = mindexes))
  elfreqs <- table(mindexes)
  if(isTRUE(decide)){
    # if(length(elfreqs) == 2) bestfit <- mean(as.numeric(names(elfreqs))) # if only two points
    # maxes <- as.numeric(names(elfreqs)[elfreqs == elfreqs[which.max(elfreqs)]])
    # if(length(maxes) > 1) bestfit <- mean(maxes) # get between more freq points identified as elbows
    # # maybe hard to decide, so agree on the more frequent and the lowest points
    # if(length(elfreqs) > 3) bestfit <- mean(c(as.numeric(names(elfreqs[1])), max(maxes)))
    bestfit <- round(mean(as.numeric(names(elfreqs))))
  }else{ bestfit <- mean(mindexes) }
  qposition <- data.frame(y = max(y), N = elfreqs[as.character(mindexes)])
  colnames(qposition) <- c("y", "x", "N"); qposition$x <- as.numeric(as.character(qposition$x)) + 0.25
  list(elbows = mindexes, ptable = qposition, elbow = bestfit)
}

## replace categoris of a column in another
creplace <- function(x, cin, cfrom, subs = NULL, newname = NULL){
  if(is.null(subs)) subs <- unique(x[, cin])
  if(is.null(newname)) newname <- 'newname'
  # print(table(x[, cin]))
  x[, newname] <- x[, cin]
  x[x[, cin] %in% subs, newname] <- x[x[, cin] %in% subs, cfrom]
  # print(table(x[, newname]))
  x
}

## match threshold/function length to n parameters/entries
transformations <- function(x, y){
  if(length(x) > 1 && length(y) != length(x))
    stop(length(y), " features/parameters ", length(x), " thresholds")
  if(length(y) > 1 && length(x) == 1){
    # cat("Repeating threshold\n")
    x <- x[1:length(y)]; x[is.na(x)] <- x[1] # min(x, na.rm = TRUE)
  }
  names(x) <- y; return(x)
}

findsample <- function(x, sour, ln = NULL, v = FALSE){
  if(is.null(ln)) ln <- length(sour)
  gg <- getfound(x, sour, v = v)
  if(!length(gg)){
    if(v) cat('Sampling... ')
    set.seed(myseed); gg <- sample(sour, ln)
    if(v) cat('Got:', commas(gg), '\n')
  }
  gg
}

sample_grp <- function(
  annot,
  cname = 1,
  maxln = NULL, # if negative, returns that number per group
  v = FALSE
){
  annot <- remove.factors(annot)
  grsize <- table(annot[, cname])
  factored <- rep(min(grsize), length(grsize)) # take the smallest group size
  if(!is.null(maxln)){
    if(is.character(maxln)){ # sample to a total of cells
      maxln <- as.numeric(gsub("[A-z]", "", maxln))
      maxln <- round(ifelse(maxln / nrow(annot) > 1, 1, maxln / nrow(annot)), 2)
    }; #print(maxln)
    maxln <- min(c(maxln, max(grsize))) # can't take more than the biggest
    if(maxln > 1){
      maxln <- (maxln / max(grsize)) # maximum size per group
    }else if(maxln[1] < 0){
      maxln <- (abs(maxln) / grsize)
    }
    if(maxln[1] == 0) maxln <- 1
    factored <- round(maxln * grsize) # take a percentage
    factored <- ifelse(factored > grsize, grsize, factored)
  }
  names(factored) <- names(grsize); set.seed(27)
  scells <- lapply(names(grsize), function(x){
    sample(rownames(annot[annot[, cname] == x, ]), factored[x])
  }); names(scells) <- names(grsize)
  finalgrsize <- sapply(scells, length)
  scells <- unlist(scells)
  if(v){
    cat("Given group:", commas(names(grsize)), "\n")
    cat("Init. sizes:", commas(grsize), "\n")
    cat("Final sizes:", commas(finalgrsize), "\n")
    cat("% from total:", (sum(finalgrsize) / nrow(annot)) * 100, "\n")
    cat("Returning:", length(scells), "\n")
  }
  scells
}

make_title <- function(x){
  y <- gsub("orig|\\.", "_", casefold(x, upper = TRUE), ignore.case = TRUE)
  y <- gsub("_", " ", gsub("_{2,}", "_", y))
  y <- gsub(" {1,}", " ", y)
  y <- gsub(" $|^ ", "", y)
  gsub("_$|^_", "", y)
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

# find words in file
infile <- function(x, this_file, write_to = FALSE, ...){
  fnames <- lapply(this_file, head, 1); x <- unique(x)
  found <- lapply(this_file, function(fname){
    cat('Fetching:', basename(fname[1]), '\n')
    filetext <- readfile(fname[1], v = T, ...)
    if(is.na(fname[2])) fname[2] <- 'gene_name'
    if(fname[2] %in% colnames(filetext)){
      cat('Column:', fname[2], '\n')
      filtered <- filetext[filetext[, fname[2]] %in% x, ]
    }else{
      rowco <- apply(filetext, 2, function(tt) which(tt %in% paste0("'", genes)) )
      rowco <- data.table::melt(rowco[sapply(rowco, length) > 0])
      cat('Column(s):', commas(unique(rowco[, 2])), '\n')
      filtered <- filetext[rowco[, 1], ]
    }
    #print(headmat(filtered))
    if(isTRUE(write_to)){
      newfname <- sub('\\.csv', '', basename(fname[1]))
      newfname <- paste0(names(fnames[fnames == fname[1]]), '_', length(x),'found_in_', newfname, '.csv')
      cat('Writing to:', newfname, '\n')
      write.csv(filtered, file = newfname)
    }; cat('\n')
    filtered
  })
}

# from vector, get these patterns
getpats <- grepatrn <- function(x, accept, reject = '1random302'){
  x <- x[grepl(pattern = paste0(accept, collapse = "|"), x = x, ignore.case = TRUE)]
  x[!grepl(pattern = paste0(reject, collapse = "|"), x = x, ignore.case = TRUE)]
  # greps <- paste0(paste0("grepl('", accept, "', x)"), collapse = " | ")
  # x <- eval(expr = parse(text = paste0("x[", greps, "]")))
  # greps <- paste0(paste0("grepl('", reject, "', x)"), collapse = " | ")
  # eval(expr = parse(text = paste0("x[!(", greps, ")]")))
}

## select columns with a range of number of elements
select_columns <- function(df, ...){
  cnames <- sapply(df, function(x) length(table(x)) )
  names(cnames)[cnames > 1 & cnames < 200 & grepatrn(names(cnames), ...)]
}

# Create a column for cells expressing a gene
gettag <- function(gg, sig = '+', prefix = FALSE){
  if(is.list(gg)){ cat("'gg' already a list\n"); return(gg) }
  lapply(gg, function(x){
      if(isTRUE(prefix)) return(c(paste0('tag_', x), paste0(sig, x)))
      c(paste0('tag_', x), paste0(x, sig))
  })
}
add_gene_tag <- function(
  lgenes,
  annot,
  mat,
  thresh = 0,
  tag = c('tag', '+', '-'),
  prefix = FALSE,
  v = FALSE
){
  if(length(thresh) != length(lgenes)){ # check if there is a threshold for each gene
    if(length(thresh) == 1) thresh <- rep(thresh, length(lgenes))
    if(!is.null(names(thresh))) thresh <- thresh[lgenes]
    thresh <- thresh[1:length(lgenes)]
    if(sum(is.na(thresh))) warning('No threshold specified for ', commas(lgenes[is.na(thresh)]), '; using 0\'s')
    thresh[is.na(thresh)] <- 0
  }
  names(thresh) <- lgenes
  newcolms <- cbindList(lapply(lgenes, function(thisgene){
    if(v) cat('----\nGene:', thisgene, '- threshold:', thresh[thisgene], '\n')
    if(prefix){
      tags <- paste0(tag, "_", casefold(thisgene, upper = TRUE))
    }else{
      tags <- paste0(tag[1], "_", casefold(thisgene, upper = TRUE))
      tags <- c(tags, paste0(casefold(thisgene, upper = TRUE), tag[-1]))
    }
    annot[, tags[1]] <- ifelse(as.vector(t(mat[thisgene, rownames(annot)] > thresh[thisgene])), tags[2], tags[3])
    tvar <- rownames(annot[annot[, tags[1]] == tags[3], ]) # check values
    tvar <- unlist(mat[thisgene, tvar])
    # if(v){ cat("Negative summary:\n"); print(summary(tvar)) }
    tvar <- rownames(annot[annot[, tags[1]] == tags[2], ]); tvar <- unlist(mat[thisgene, tvar])
    # if(v){ cat("Positive summary:\n"); print(summary(tvar)) }
    if(v){ cat('Proportions'); print(table(annot[, tags[1]])) }
    return(annot[, tags[1], drop = FALSE])
  }))
  newcolms
}

# New Read10X
Read10Xnc <- function (data.dir = NULL)
{
    full.data <- list()
    for (i in seq_along(data.dir)) {
        run <- data.dir[i]
        if (!dir.exists(run)) {
            stop("Directory provided does not exist")
        }
        if (!grepl("\\/$", run)) {
            run <- paste(run, "/", sep = "")
        }
        barcode.loc <- paste0(run, c("barcodes.tsv.gz", "barcodes.tsv"))
        gene.loc <- paste0(run, c("features.tsv.gz", "genes.tsv"))
        matrix.loc <- paste0(run, c("matrix.mtx.gz", "matrix.mtx"))
        barcode.loc <- barcode.loc[file.exists(barcode.loc)]
        gene.loc <- gene.loc[file.exists(gene.loc)]
        matrix.loc <- matrix.loc[file.exists(matrix.loc)]
        if (!sum(file.exists(barcode.loc))) {
            stop("Barcode file missing")
        }
        if (!sum(file.exists(gene.loc))) {
            stop("Gene name file missing")
        }
        if (!sum(file.exists(matrix.loc))) {
            stop("Expression matrix file missing")
        }
        data <- readMM(file = matrix.loc)
        cell.names <- readLines(barcode.loc)
        gene.names <- readLines(gene.loc)
        if (all(grepl(pattern = "\\-1$", x = cell.names))) {
            cell.names <- as.vector(x = as.character(x = sapply(X = cell.names,
                FUN = ExtractField, field = 1, delim = "-")))
        }
        rownames(x = data) <- make.unique(names = as.character(x = sapply(X = gene.names,
            FUN = ExtractField, field = 2, delim = "\\t")))
        if (is.null(x = names(x = data.dir))) {
            if (i < 2) {
                colnames(x = data) <- cell.names
            }
            else {
                colnames(x = data) <- paste0(i, "_", cell.names)
            }
        }
        else {
            colnames(x = data) <- paste0(names(x = data.dir)[i],
                "_", cell.names)
        }
        full.data <- append(x = full.data, values = data)
    }
    full.data <- do.call(cbind, full.data)
    return(full.data)
}

# to check on a table the proportions of ONE column in others
explore_feature <- function(df, cnames, cconst){
  summ <- data.frame(t(cbindList(lapply(cnames, function(x){
    mytab <- table(df[, x], df[, cconst])
    tvar <- dimnames(mytab)
    mytab <- data.frame(as.matrix.data.frame(mytab), row.names = tvar[[1]])
    colnames(mytab) <- tvar[[2]]
    mytab <- cbind(type = rep(sub("orig.", "", x), nrow(mytab)), mytab)
    t(mytab)
  }))), stringsAsFactors = FALSE, check.names = FALSE)
  summ[, -1] <- sapply(summ[, -1], as.numeric)
  summ
}

new_features <- function(elist, cname = "nNew_Features"){
  mytab <- rbind(data.frame(PC = 1, nNew_Features = length(elist[[1]]), Features = commas(elist[[1]], Inf)),
    data.frame(data.table::rbindlist(lapply(2:length(elist), function(x){
      y <- sum(!elist[[x]] %in% unlist(elist[1:(x - 1)]))
      z <- commas(elist[[x]][!elist[[x]] %in% unlist(elist[1:(x - 1)])], Inf)
      if(!y) z <- ''
      data.frame(PC = x, nNew_Features = y, Features = z)
    })))
  )
  colnames(mytab) <- sub("nNew_Features", cname, colnames(mytab))
  colnames(mytab) <- sub("Features", sub(".*_(.*)", "\\1", cname), colnames(mytab))
  mytab
}

# get top genes and highly expressed from a comparison
get_markers <- function(
  fname,
  meang = NULL,
  pval = 0.05,
  fch = 2,
  pvtype = 'padj',
  lfc.type = 'log2FoldChange',
  nmean = 0.8,
  prefix = NULL,
  do.plot = TRUE
){
  genestab <- read.csv(fname, stringsAsFactors = F, check.names = F)
  eline <- paste("setorder(genestab,",  pvtype, ")"); eval(expr = parse(text = eline))
  headmat(genestab, 10)
  rownames(genestab) <- sub("'", "", genestab[, 'gene_name'])
  # print(colnames(genestab))
  genes <- getDEGenes(genestab, pv = pval, fc = fch, upreg = TRUE, pvtype = pvtype, lfc.type = lfc.type, th = "", v = TRUE)
  tvar <- genestab[genes, ]
  meang <- head(colnames(tvar)[grepl(meang, colnames(tvar))], 1)
  cat('Filtering for', meang, '\n')
  tvar <- tvar[tvar[, meang] > nmean, ]
  # eline <- paste("setorder(tvar,",  lfc.type, ")"); eval(expr = parse(text = eline))
  tvar <- tvar[order(tvar[, lfc.type]), ]
  mygenes <- rownames(tvar)
  datp = data.table::melt(data.frame(N = 1:nrow(tvar),
    genes = mygenes,
    FC = tvar[, lfc.type],
    FDR = -log10(tvar[, pvtype]), stringsAsFactors = F),
  id.vars = c('N', 'genes'))
  head(datp)
  if(isTRUE(do.plot)){
    pdf(paste0('threshold_dinamic_', length(mygenes), 'genes', prefix, '.pdf'), height = 5, width = 8)
    print(ggplot(datp, aes(x = N, y = value)) + geom_point() +
      facet_wrap(~ variable, ncol = 2, scales = 'free_y') +
      ggtitle(label = paste('N:', length(mygenes)), subtitle = paste('qval:', pval, '\nFC:', fch, '\nMean:', nmean)) +
      theme(strip.background = element_rect(fill = "transparent", linetype = 0)))
    dev.off()
  }
  cat('My genes', commas(mygenes), '\n')
  mygenes
}

create_pairs <- function(x){
  mycombs <- lapply(x, function(y){
    groups <- unique(y)
    myclass <- if(length(groups) < 30){
      gtools::combinations(length(groups), r = 2, v = groups, set = TRUE, repeats.allowed = FALSE)
    }else{
      matrix(1, ncol = 2)
    }
    data.frame(myclass)
  })
  colname <- rep(names(mycombs), sapply(mycombs, nrow))
  mycombs <- data.frame(data.table::rbindlist(mycombs))
  mycombs$Category <- colname
  mycombs$Name <- "void"
  colnames(mycombs) <- c("Class_1", "Class_2", "Category", "Name")
  mycombs
}

# to compare approaches
# library(microbenchmark)
# res <- microbenchmark(
#     sweep_dense = sweep(mat, 1, vec, '*'),
#     sweep_sparse = sweep(mat_sparse, 1, vec, '*'),
#     mult_dense = mat * vec,
#     mult_sparse = mat_sparse * vec
# ); autoplot(res)

clust_similarity <- function(df){
  df[, 1] <- as.character(df[, 1])
  df[, 2] <- factor(df[, 2])
  df[, 3] <- factor(df[, 3])
  list(
    JI = clusteval::cluster_similarity(labels1 = as.numeric(df[, 2]), labels2 = as.numeric(df[, 3])),
    ARI = mclust::adjustedRandIndex(x = df[, 2], y = df[, 3]),
    NMI = NMI::NMI(X = df[, 1:2], Y = df[, c(1, 3)])$value
  )
}

# Collapse column summarising the others
column_collapse <- function(metab, rname = 1, v = FALSE){
  if(!any(duplicated(metab[, rname]))){
    if(v) cat("."); return(metab) # return if already
  }; if(v) cat("s")
  gnames <- unique(metab[, rname][duplicated(metab[, rname])])
  colsum <- colnames(metab)[!colnames(metab) %in% rname] # columns to summarise
  # check if column is same length
  csize <- sapply(colsum, function(co){ ttab <- table(metab[, c(co, rname)]); nrow(ttab) == ncol(ttab) })
  cclass <- sapply(metab[, colsum, drop = FALSE], class) == "character" # keep numerics
  colsum <- names(which(!(csize & cclass))) # eliminate repeated
  # myformula <- paste(rname, "~", paste(colsum, collapse = " + "))
  # y <- dcast.data.table(setDT(metab), myformula)
  if(v) cat(length(gnames))
  sx <- data.frame(data.table::rbindlist(lapply(gnames, function(g){ # g <- gnames[1]; g <- "'CISH"
    per_g <- lapply(colsum, function(cname){ # cname <- colsum[1]
      y <- metab[which(metab[, rname] == g), cname]
      y <- if(is.numeric(y)) round(y, 4) else y
      return(paste0(unique(y), collapse = " & "))
    })
    per_g <- data.frame(per_g, stringsAsFactors = FALSE)
    colnames(per_g) <- colsum
    per_g
  })), stringsAsFactors = FALSE, row.names = gnames)
  sx[, rname] <- gnames
  sx <- list(sx, metab[!metab[, rname] %in% gnames, colnames(sx)])
  sx <- data.frame(data.table::rbindlist(sx), stringsAsFactors = FALSE)
  sx
}

# summarise tables by a column
summ_tables <- function(
  ltabs,
  column2row = "gene_name",
  tab4roworder = NULL,
  colgroup_trans = NULL,
  colexclude = NULL,
  v = FALSE
) {
  # ltabs <- list(res, statseffclusters, res_seu)
  # colgroup_trans = list(cluster = newlabs)
  if(v) cat("Columns:", commas(sapply(ltabs, ncol)),
            "\nRows:", commas(sapply(ltabs, nrow)), "\n")
  if(!is.null(colexclude)){
    colexclude <- paste0(colexclude, collapse = "|")
    if(v) cat("Excluding patterns in columns", colexclude, "\n")
  }
  ltabs <- lapply(ltabs, remove.factors)
  if(v) cat("Adding", column2row, "from rows when needed: ")
  ltabs <- lapply(ltabs, function(x){
    if(!column2row %in% colnames(x)){
      if(v) cat("r"); x[, column2row] <- paste0("'", rownames(x))
    }else if(v) cat(".")
    if(!is.null(colexclude)) x <- x[, colnames(x)[!grepl(colexclude, colnames(x))], drop = FALSE]
    x
  }); if(v) cat("\n")
  if(!is.null(colgroup_trans)){
    if(v) cat("Changing groups in column(s):", commas(names(colgroup_trans)), "\n")
    ltabs <- lapply(ltabs, function(x){
      cnames <- colnames(x)[colnames(x) %in% names(colgroup_trans)]
      for(cname in cnames){
        x[, cname] <- unname(colgroup_trans[[cname]][as.character(x[, cname])])
      }
      x
    })
  }
  if(v) cat("Removing repeated column names\n")
  for(i in 2:length(ltabs)){
    cnames <- unique(unlist(sapply(ltabs[1:(i - 1)], colnames)))
    tvar <- colnames(ltabs[[i]]) %in% cnames & !colnames(ltabs[[i]]) %in% column2row
    cat("Table", i, "- repeated:", sum(tvar), "\n")
    ltabs[[i]] <- ltabs[[i]][, !tvar]
  }
  if(v) cat("Summarise tables with repeated", column2row, ": ")
  ltabs_sum <- lapply(ltabs, function(mytab){ # mytab <- ltabs[[1]]
    column_collapse(metab = mytab, rname = column2row, v = v)
  }); if(v) cat("\n")
  # Fix the rownames
  # sapply(ltabs_sum, function(x) head(rownames(x)) )
  if(v) cat("Setting", column2row, "as rows\n")
  ltabs_sum <- lapply(ltabs_sum, function(x) setrows(tab = x, x = column2row) )
  if(v) cat("Final columns:", commas(sapply(ltabs_sum, ncol)),
            "\nFinal rows:", commas(sapply(ltabs_sum, nrow)), "\n")
  # bind the tables
  cnames <- unique(c(column2row, unlist(sapply(ltabs_sum, colnames))))
  tab4roworder <- unique(c(tab4roworder, 1:length(ltabs_sum)))
  rnames <- unique(unlist(lapply(ltabs_sum[tab4roworder], rownames)))
  if(v) cat("Columns, rows: ", length(cnames), ", ", length(rnames), "\n", sep = "")
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

## Find run tools' versions
findtools <- function(runpath){
  # runpath <- "/mnt/NGSAnalyses/RNA-Seq/Mapping/004013_AS3_Expt4_RSS60_Bulk/main.log"
  # tools <- c("BOWTIE_PATH", "TOPHAT", "FASTQC_EXE", "SAMTOOLS", "HTSEQ_COUNT", "FASTX_TRIMMER", "CUTADAPT")
  tools <- c("BOWTIE", "TOPHAT", "FASTQC", "SAMTOOLS", "HTSEQ_COUNT", "FASTX_TRIMMER", "CUTADAPT")
  tools <- paste0("^", casefold(tools, upper = TRUE), ".*=")
  command <- paste(paste(c("grep", tools), collapse = ' -e '), runpath)
  system(command)
}

splitString <- function (text, width) {
  strings <- strsplit(text, " ")[[1]]
  newstring <- strings[1]
  linewidth <- stringWidth(newstring)
  gapwidth <- stringWidth(" ")
  availwidth <- convertWidth(width, "in", valueOnly = TRUE)
  for (i in 2:length(strings)) {
    width <- stringWidth(strings[i])
    if (convertWidth(linewidth + gapwidth + width, "in",
                     valueOnly = TRUE) < availwidth) {
      sep <- " "
      linewidth <- linewidth + gapwidth + width
    }
    else {
      sep <- "\n"
      linewidth <- width
    }
    newstring <- paste(newstring, strings[i], sep = sep)
  }
  newstring
}

run_intervals <- function(
  command, byinter = "5 sec", n = 2,
  intervals = seq(Sys.time(), by = byinter, length.out = n)
){
  #print(data.frame(Times = intervals))
  cat('Run', n, 'each', byinter, 'times\n')
  p1 <- proc.time()
  for(i in 2:length(intervals)) {
    cat("---------------- Executing ---------------\n")
    timestamp()
    # if(i != 2)
    system(command) # remember to remove the if
    sleepTime <- intervals[i] - intervals[i - 1]
    tvar <- sleepTime > 0
    timestamp()
    if (is.na(tvar)){ cat("Broken interval\n"); break }
      cat("Sleeping for a "); print(sleepTime)
    if (tvar) Sys.sleep(sleepTime)
    tvar <- Sys.time() < intervals[i]
    if(isTRUE(tvar)){
      # tvar <- Sys.time() + sleepTime
      cat('Using while...\n')
      while(Sys.time() < intervals[i]){}
    }
    gc()
  }
  cat("Finished at:", date(), "\n")
  print(proc.time() - p1)
}

# Get a row's corresponding matches and order them
row2tab <- function(x, selection = NULL){
  if(is.null(selection)){
    selection <- which(x == max(x), arr.ind = TRUE)
    selection <- rownames(selection)[1]
  }
  y <- reshape2::melt(x[selection, ])
  y <- y[order(y[, ncol(y)], decreasing = TRUE), , drop = FALSE]
  # if(!all(rownames(y) %in% 1:nrow(y))) y$x <- rownames(y)
  # rownames(y) <- NULL
  colnames(y) <- c(selection)#, colnames(y)[-1])
  y
}

# get columns with a max N of groups
get_grouping_cols <- function(
  metadat,
  onames = NULL,
  maxn = 27,
  ckeep = NULL, cignore = NULL,
  types = "character",
  na_rm = FALSE,
  v = FALSE
){
  metadat <- remove.factors(metadat)
  if(is.null(onames)) onames <- colnames(metadat)
  if('numeric' %in% types) maxn <- Inf
  onames <- rev(onames)
  if(v) cat("Selecting types:", paste0(types, collapse = ", "), '\n')
  tvar <- sapply(metadat[, onames, drop = FALSE], class) %in% types # if it's the type and has no NA's
  if(isTRUE(na_rm)) tvar <- tvar & !sapply(metadat[, onames, drop = FALSE], function(x) any(is.na(x)) )
  onames <- onames[unname(tvar)]
  if(!is.null(ckeep)){
    if(length(get_grouping_cols) > 1){
      onames <- onames[onames %in% ckeep]
    }else{ onames <- onames[grepl(ckeep, onames, ignore.case = TRUE)] }
  }
  if(!is.null(cignore)) onames <- onames[!grepl(cignore, onames, ignore.case = TRUE)]
  if(length(onames) == 0){ if(v) cat("Taking all columns\n"); onames <- colnames(metadat) }
  if(v) cat("Filtering", length(onames), "\n")
  nnames <- sapply(onames, function(x) length(unique(metadat[, x])) )
  if(v) cat(sum(duplicated(nnames)), "same N across columns\n")
  onames <- names(sort(nnames))
  onames <- onames[nnames[onames] != nrow(metadat)] # excluding  N = total rows
  onames <- onames[nnames[onames] > 1] # excluding N = 1
  onames <- onames[nnames[onames] <= maxn] # N <= 27 because ggpairs doesn't allow more than this
  if(v) cat("Match 2-group combination:", length(onames), "\n")
  pairs <- try(gtools::combinations(n = length(onames), r = 2, v = onames))
  if(class(pairs) != 'try-error'){
    kha <- apply(pairs, 1, function(x){
      y <- table(metadat[, c(x)]) # as.data.frame.matrix
      diag(y) <- 0; any(rowSums(y) > 0) # take when none but diagonal are 0, aka, same N groups
    })
    if(v && any(!kha)) cat(sum(!kha), "matching 1-to-1 groups:", paste0(pairs[!kha, 2], collapse = ", "), "\n")
    onames <- onames[!onames %in% pairs[!kha, 2]]
  }else{
    if(v) cat("n =", length(onames), "/ unique =", length(unique(onames)), "-", paste0(unique(onames), collapse = ', '), "\n")
  }
  kha <- apply(t(apply(metadat[, onames, drop = FALSE], 1, duplicated)), 2, sum) # get dups across columns
  onames <- onames[(kha != nrow(metadat))]
  if(v) cat("Returning:", length(onames), "\n")
  onames
}

# add gene set percentages
add_pcts <- function(
  mymeta,
  adata,
  gene_pattern = list(
    pctMitochondrial = c('^mt-', '^m-', '^hg19_mt-', '^mm10_mt-'),
    pctRibosomal = paste0(c('^rps', '^rpl'), "[0-9]"),
    pctHSprefix = c("^hsp[0-9]"),
    pctHeatShock = paste0("^hsp", letters[1:5])
  ),
  v = FALSE
){
  if(v) cat("Total samples:", nrow(mymeta), "\n")
  adata <- adata[, rownames(mymeta)]
  gnames <- parse_ens_name(rownames(adata))
  tots <- colSums(adata)
  void <- sapply(names(gene_pattern), function(pname) getpats(x = gnames, accept = gene_pattern[[pname]]) )
  void <- void[sapply(void, length) > 0]
  if(length(void) > 0){
    if(v) cat("Retrieving sets:\n")
    if(v) str(void)
    pcts <- sapply(void, function(x){
      gene_patterned <- getpats(x = rownames(adata), accept = x)
      y <- round(Matrix::colSums(x = adata[gene_patterned, , drop = FALSE])/tots * 100, 2)
      y[is.na(y)] <- max(y, na.rm = TRUE)
      y
    })
    mymeta <- cbind_repcol(pcts, mymeta)
  }
  mymeta
}

# get quantiles
quantileze <- function(
  x,
  probies = seq(0, 1, 0.25),
  names_use = c('q', 'v') # quantiles or values
){
  names_use <- match.arg(names_use)
  y <- unique(x)
  mybreaks <- quantile(y, probs = probies)
  mybreaks <- mybreaks[!duplicated(mybreaks)]
  tvar <- switch(names_use,
    'q' = as.numeric(sub("%", "", names(mybreaks))),
    'v' = round(unname(mybreaks), 2)
  )
  cut(
    x,
    breaks = unname(mybreaks),
    labels = paste0("Q", levels(cut(tvar, tvar))),
    include.lowest = TRUE
  )
}

# https://univ-nantes.io/E114424Z/veneR/blob/master/RNASeq.R
retrieveSexHumanEmbryoKmeans <- function(
  d,
  group = NULL,
  v = FALSE
){
    if(!is.null(group)) group <- droplevels(as.factor(group))
    #return a matrix of count of "male" and "female" predicted cells in each embryo (Kmeans method)
    maleGene <- c("DDX3Y", "EIF1AY", "TTTY15", "RPS4Y1", "RPS4Y2", "SRY")
    tvar <- grepl(paste(paste0("^", maleGene, "$"), collapse = "|"), parse_ens_name(rownames(d)), ignore.case = TRUE)
    maleGene <- rownames(d)[tvar]
    if(v) cat("Using genes:", paste0(maleGene, collapse = ", "), "\n")
    k <- kmeans(t(d[maleGene, ]), 2)
    mORf <- rowSums(k$centers)
    if(mORf[1] < mORf[2]){
      mf <- c("F","M")
    }else{
      mf <- c("M","F")
    }
    if(is.null(group)) return(mf[k$cluster])
    count <- as.factor(mf[k$cluster])
    embryos <- as.factor(levels(group))
    names(count) <- group
    res <- list()
    res$count <- data.frame(matrix(ncol = 2, nrow = length(embryos)))
    colnames(res$count) <- c("Male","Female")
    rownames(res$count) <- embryos

    for(embryo in embryos){
      res$count[embryo,1] <- length(which(count[which(names(count)==embryo)]=="M"))
      res$count[embryo,2] <- length(which(count[which(names(count)==embryo)]=="F"))
    }
    res$freq <- res$count/rowSums(res$count)
    res$pred <- count
    return(res)
}

## Generating colours for a table
colours_from_metadata <- function(
  mytab
){
  tvar <- sapply(mytab, is.character) & sapply(mytab, function(x) length(unique(x)) ) != nrow(mytab)
  x <- unlist(mytab[, tvar])
  mynames <- unique(x); mynames <- mynames[!is.na(mynames)]
  mycouls <- WGCNA::labels2colors(mynames)
  ddf <- data.frame(group = mynames, colour = mycouls, stringsAsFactors = FALSE)
  head(ddf); tail(ddf)
  ddf
  # precols <- read.csv(fcouls, row.names = 1, stringsAsFactors = F)
  # ddf[ddf$colour %in% precols[, 1], ]
  # write.table(ddf, file = fcouls_out, row.names = FALSE, quote = TRUE, sep = ",")
  # system(paste("head", fcouls_out))
  # system(paste("tail", fcouls_out))
  # postcols <- read.csv(fcouls_out, row.names = 1, stringsAsFactors = F)
  # head(postcols); tail(postcols)
}

# order vector following a reference order
this_order <- function(x, ref, which = FALSE){
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
orderthis <- function(
  annot,
  order_by = 'corder',
  cname = NULL,
  grupos = NULL,
  mat = NULL,
  lgenes = NULL,
  samplit = FALSE,
  v = FALSE
){
  annot$corder <- 'ALL'
  if(!order_by %in% c('pca', 'mean', 'hc', colnames(annot))) order_by <- "corder"
  if(v) cat('Setting columns order by:', commas(order_by), '\n')
  if(order_by == 'pca' && is.null(mat)){
    stop("You need an expression matrix, 'mat'")
    return(rownames(annot))
  }
  if(is.null(cname)) cname <- order_by[1]
  if(v) cat('Covariates:', commas(cname, Inf), '\n')
  minsiz <- min(table(annot[, cname]))
  if(isTRUE(samplit)){
    if(v) cat("Sampling\n")
    if(!is.numeric(samplit)) samplit <- NULL else samplit <- -samplit
    annot <- annot[sample_grp(annot, cname = cname, maxln = samplit), ]
  }
  grupos <- grupos[grupos %in% as.character(annot[, cname])]
  if(length(grupos) == 0){ warning("Given groups not in data"); grupos <- NULL }
  if(is.null(grupos)) grupos <- gtools::mixedsort(unique(as.character(annot[, cname])))
  if(v) cat('Groups:', commas(grupos, Inf), '\n')
  if(is.null(lgenes) && !is.null(mat)){
    lgenes <- rep(list(rownames(mat)), length(grupos))
    names(lgenes) <- grupos
  }
  list_order <- lapply(grupos, function(x){
    mynames <- getsubset(c(cname, x), annot, v = v)
    if(v) cat(x, "")
    myorder <- ifelse(isTRUE(order_by %in% colnames(annot) && order_by != "corder"), 'acol', order_by)
    if(v && (!myorder %in% c("acol", "corder"))) cat(myorder); if(v) cat("\n")
    if(cname == order_by) return(mynames)
    if(!is.null(mat)) mygenes <- getfound(lgenes[[x]], rownames(mat), v = FALSE)
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

# Filter a table — adding columns... maybe complex designs
meta_filtering <- function(
  mdata,
  filters = "none",
  cname = "none",
  sepchar = "~",
  v = FALSE
) {
  if(filters == "none" && cname == "none") return(list(annotation = mdata))
  if(v){ cat('Using filter:\n'); str(filters) }
  if(sum(grepl(sepchar, filters)) || length(filters) > 1 || is.list(filters)){
    filters <- translist(filters)
  }

  # subsetting annotation first checking for a file
  if(file.exists(filters[[1]][1])){
    if(v) cat("Filters from file", filters[[1]][1], "\n")
    filters <- list(filters[[1]][1], filters[[1]][-1]) # separating the file from the rest of the filters
    filecon <- readfile(filters[[1]][1], stringsAsFactors = FALSE, v = v); if(v) str(filecon)
    tvar <- head(which(sapply(filecon, function(x) any(x %in% rownames(mdata)) )), 1)
    if(length(tvar) == 0) warning("No column in ", basename(filters[[1]][1]), " was compatible with 'mdata'")
    rownames(filecon) <- filecon[, tvar]
    tvar <- getfound(rownames(filecon), rownames(mdata), element = "cell/sample", v = v)
    mdata <- mdata[tvar, ]
    mdata <- joindf(mdata, filecon)
    if(length(filters) > 1) filters <- filters[-1]
    # if(length(filters[[1]]) > 1) filters <- lapply(filters, function(x) x[-1] )
    filterssuffix <- paste0("_", basename(filters[[1]][1]))
  }

  tvar <- sapply(filters, function(x) grepl(sepchar, x[1]) )
  if(sum(tvar) > 1) stop("currently supporting only one combination")
  if(sum(tvar) == 1){
    filterst <- list()
    for(x in filters[tvar]){
      if(v) cat(x[1], "\n")
      thislist <- as.list(data.frame(unname(t(sapply(x, function(y) unlist(strsplit(y, sepchar)) ))), stringsAsFactors = FALSE))
      filterst <- unname(c(filterst, thislist))
    }
    filters <- c(filterst, filters[!tvar], list(c("combn", filters[[tvar]][-1])))
    addeds <- sapply(filterst, head, 1) #sapply(filterst, function(x) ifelse(!grepl("^tag", x[1]), x, NA) )
    addeds <- addeds[!is.na(addeds)]
  }else{
    addeds <- NULL
  }

  # Add gene tag
  # cname <- "tag_CXCL10_10~tag_FOXP3_10~tag_CD4_10"
  # cname <- "tag_IL9"
  # addgene <- "tag_IL9_0"
  for(addgene in unique(c(cname, sapply(filters, head, 1)))){
    if(addgene %in% colnames(mdata)) next # if it's not in the annotation already
    gg <- unlist(strsplit(gsub("tag_|_[0-9]+", "", addgene), sepchar))
    if(grepl("tag", addgene) && all(gg %in% rownames(cts))){
      if(v) cat("Adding gene tags\n")
      tvar <- as.numeric(gsub("tag_.*_", "", unlist(strsplit(addgene, sepchar))))
      tvar <- ifelse(is.na(tvar), 0, tvar)
      void <- add_gene_tag(lgenes = gg, annot = mdata, mat = cts, thresh = tvar, tag = c('tag', 'p', 'n'), v = v)
      colname <- addgene#paste0("tag", gsub(paste0("tag|_[0-9]+|", sepchar), "", addgene))
      void[, colname] <- apply(void, 1, paste, collapse = "")
      print(reshape2::melt(table(void[, colname])))
      mdata <- cbind(void, mdata); headmat(mdata)
      addeds <- c(addeds, colname)
    }
  }
  tvar <- unique(c(cname, sapply(filters, head, 1)))
  if((!is.null(addeds) && ('combn' %in% tvar)) || grepl(sepchar, cname)){
    addeds <- if(is.null(addeds)) unlist(strsplit(cname, sepchar)) else unique(addeds)
    if(v) cat("Combining columns", addeds, "\n", sep = " ")
    # filters <- filters[!sapply(filters, head, 1) %in% addeds] # remove added from filters
    mdata$combn <- apply(mdata[, addeds, drop = FALSE], 1, paste, collapse = "_")
    if(grepl(sepchar, cname)) colnames(mdata) <- sub("^combn$", cname, colnames(mdata))
  }

  # now if there's further filtering
  # watch out for 'cname' column samples/cells by table(cname, column_filtering_by)
  if(any(sapply(filters, head, 1) %in% colnames(mdata))){
    if(v){ cat("Filtering instructions\n"); str(filters) }
    if(!is.null(cname)){ # watch out before
      void <- lapply(filters, function(x) table(mdata[, c(cname, x[1])], useNA = 'always') ); print(void)
    }
    mdata <- mdata[getsubset(filters, mdata, v = T), ]
    if(!is.null(cname)){ # watch out after
      void <- lapply(filters, function(x) table(mdata[, c(cname, x[1])], useNA = 'always') ); print(void)
    }
    tvar <- sapply(filters, function(x) paste0(paste0(x[-1], collapse = "AND"), "_from_", x[1]) ) # creates "grp1ANDgrp2"
    filters <- paste0("_", paste0(tvar, collapse = "_and_")) # creates "grp1ANDgrp2_from_COLUMN1_and_grp1_from_COLUMN1"
    if(exists('filterssuffix')){
      filters <- paste0(filterssuffix, filters); rm(filterssuffix)
    }
  }
  if(grepl("^express", filters[[1]][1])){
    filters <- sub("expr[A-z]+ ", "", filters)
    sset <- paste0("cellsf <- rownames(subset(annottab, subset = ", filters, "))")
    if(v) cat("Expression:", sset, "\n")
    eval(expr = parse(text = sset))
    mdata <- mdata[cellsf, ]
  }
  return(list(annotation = mdata, filter = filters))
}

finished_file <- function(x, size = 3620){ file.exists(x) && isTRUE(file.size(x) > size) }

# Fisher exact text overlap
test_overlap <- function(
  x,
  y = NULL,
  total = NULL,
  names = c("List1", "List2"),
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
  # colnames(ddf) <- gsub("A$", paste0("_", mynames[1]), colnames(ddf))
  # rownames(ddf) <- gsub("B$", paste0("_", mynames[2]), rownames(ddf))
  colnames(ddf) <- paste0(c('in_', 'out_'), mynames[1])
  rownames(ddf) <- paste0(c('in_', 'out_'), mynames[2])
  # chisq.test(ddf)
  if(is.null(main)) main = paste0(mynames, collapse = " vs ")
  ft <- fisher.test(ddf); if(verbose) ft
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

get_correct_root_state <- function(cds, cell_phenotype, root_type = NULL){
  catgcells <- as.character(pData(cds)[, cell_phenotype])
  if(is.null(root_type)){
    root_type <- names(head(sort(table(catgcells)), 1))
  }
  cell_ids <- which(as.character(catgcells) == root_type)

  closest_vertex <-
    cds@auxOrderingData[[cds@rge_method]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    V(cds@minSpanningTree)$name[as.numeric(names
      (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/R/stats_summary_table.R')
source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/plotting.R')
