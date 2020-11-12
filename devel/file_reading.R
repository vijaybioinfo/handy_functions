#!/usr/bin/R

# detect file format and read it
# Maybe use matchArgs in the future when passing more arguments
# futil = c(list(...), futil); futil <- futil[!duplicated(names(futil))]
readfile <- function(
  myfile,
  usedt = FALSE, # use data.table
  verbose = FALSE,
  ...
){
  if(isTRUE(usedt) && grepl("tsv$|csv$|txt$", myfile)){
    file_content <- data.frame(data.table::fread(myfile), stringsAsFactors = FALSE, check.names = FALSE)
    rownames(file_content) <- file_content[, 1]; file_content <- file_content[, -1, drop = FALSE]
  }else if(grepl("rda$|rta$|rdata$|robj$|rds$", myfile, ignore.case = TRUE)){
    if(verbose) cat(" - From RData file\n");
    file_content <- theObjectSavedIn(myfile)
  }else if(grepl("csv$", myfile, ignore.case = TRUE)){
    if(verbose) cat(" - From CSV file\n")
    file_content <- read.csv(myfile, ...)
  }else if(grepl("h5$|h5$", myfile, ignore.case = TRUE)){
    if(verbose) cat(" - 10X H5\n")
    file_content <- Seurat::Read10X_h5(myfile)
  }else if(grepl("txt$|tsv$", myfile, ignore.case = TRUE)){
    if(verbose) cat(" - From TXT file\n")
    file_content <- read.table(myfile, ...)
  }else if(grepl("loom$", myfile, ignore.case = TRUE)){
    file_content <- connect(filename = myfile, mode = "r")
  }else{ cat("Unknown format:", basename(myfile), "\n") }
  file_content
}

theObjectSavedIn <- function(
  saveFile,
  ob_name = NULL,
  verbose = FALSE
) {
  env <- new.env()
  if(verbose) cat('Reading object ')
  if(grepl('rds$', saveFile, ignore.case = TRUE)){
    if(verbose) cat("RDS\n")
    env$rds <- readRDS(saveFile)
  }else{ if(verbose) cat("RDATA\n"); load(saveFile, envir = env) }
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
  if(verbose) cat('Taking:', ob_name, '\n')
  if(!sum(ob_name %in% loadedObjects)) stop('Object not found')
  env[[ob_name]]
}

# get arguments
matchArgs <- function(x, func){
  x <- x[names(x) %in% methods::formalArgs(func)]
  if(!length(x)) return(NULL)
  x[!duplicated(names(x))]
}

# find a file upstream in a path
findfile <- function(name, path = getwd(), read = FALSE, nup = 4, verbose = FALSE, ...){
  path <- dircheck(path)
  tmp <- ''
  for(i in 1:nup){
    myfile <- paste0(path, tmp, name)
    tmyfile <- system(paste('ls', myfile), intern = T)
    if(length(tmyfile) > 0) myfile <- tmyfile
    if(file.exists(myfile)){
      if(verbose) cat(basename(myfile), 'found!\n')
      if(read) return(readfile(myfile, verbose = verbose, ...)) else return(myfile)
    }
    tmp <- paste0(tmp, '../'); # print(myfile)
  }
  if(verbose) cat(name, 'not found\n')
  return(data.frame(lol = 'void', stringsAsFactors = F))
}
