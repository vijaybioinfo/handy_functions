#!/usr/bin/R

###############################
# FASTQ Sample pre-processing #
###############################

# This functions are designed to manipulate sample names and even merge them

# Get the fastq file prefixes
fastq_unique_names <- function(
  mystrings,
  patties = c(
    "_S[0-9]+", # indexes
    "_L[0-9]{3}.*", # lanes
    "_R[0-9]{1,1}.*" # read
  )
){
  unique(sub("_$", "", gsub(paste0(patties, collapse = "|"), "", unique(mystrings))))
}

fastq_create_map_table <- function (
  mypath,
  verbose = TRUE
) {
  mysamples <- if(any(dir.exists(mypath))){
    if(verbose) cat(mypath, sep = '\n')
    list.files(path = mypath, full.names = TRUE)
  }else{
    mypath
  }
  unames = fastq_unique_names(mysamples)
  if(verbose){
    cat(length(unames), 'samples\n');
    cat(head(unames, 4), "...\t...", tail(unames, 4), sep = '\n')
  }
  data.frame(
    "sample ID" = basename(unames),
    fastq_f = mysamples[seq(from = 1, to = length(mysamples), by = 2)],
    fastq_r = mysamples[seq(from = 2, to = length(mysamples), by = 2)],
    check.names = FALSE
  )
}

# Create bash commands to merge lanes for prefixed sample names
fastq_merge_commands <- function(
  x,
  cmd_only = TRUE,
  rename = NULL # number of parents to rename to
){
  summaries <- lapply(x, function(mysample){
    mysample <- unlist(mysample)
    f2merge <- unlist(lapply(
      X = mysample,
      FUN = function(i) system(paste0('ls ', i, "*_R*.fastq*"), intern = TRUE)
    ))
    suffix <- sapply(f2merge, function(samp) sub(".*R[0-9]{1,1}(.*fastq.*)", "\\1", samp), USE.NAMES = FALSE)
    suffix_final <- unique(suffix)
    suffix_final <- if(length(suffix_final) > 1){
      warning("Check suffixes ", paste0(suffix_final, collapse = ", "))
      tail(suffix_final, 1)
    }else{ suffix_final }
    reads <- sapply(f2merge, function(samp){
      rs <- paste0("R", 1:2)
      names(which(sapply(rs, function(ri) grepl(paste0("_", ri, "_|_", ri, "\\."), samp) )))
    }, USE.NAMES = FALSE)
    newname <- if(is.null(rename)){
      tvar <- basename(mysample)
      sub("FQs|_{2,}", "", tvar[which.min(nchar(tvar))])
    }else if(is.numeric(rename)){
      tvar <- paste0(unique(tail(unlist(strsplit(mysample, "/")), rename)), collapse = "_")
      sub("FQs|_{2,}", "", tvar)
    }else{
      newname <- paste0(rename, newname)
    }
    newname <- paste0(newname, "_", reads, suffix_final)
    if(length(x) < 10) cat(newname, sep ='\n')
    if(length(x) < 10) cat(f2merge, sep ='\n')
    if(length(f2merge) < 3){
      cmds <- paste0('cp ', f2merge, " ", newname)
    }else{
      addition <- ifelse(duplicated(newname), " >> ", " > ")
      cmds <- paste0('cat ', mysample, "*", reads, suffix, addition, newname)
    }
    if(isTRUE(cmd_only)) return(cmds)
    list(
      found = f2merge,
      command = cmds
    )
  })
  # names(summaries) <- unname(unlist(x))
  if(isTRUE(cmd_only)) unlist(summaries) else return(summaries)
}
