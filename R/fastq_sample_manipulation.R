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
    addition <- ifelse(duplicated(newname), " >> ", " > ")
    # if(length(x) < 10) cat(newname, sep ='\n')
    # if(length(x) < 10) cat(f2merge, sep ='\n')
    if(length(f2merge) < 3){
      cmds <- paste0('cp ', f2merge, " ", newname)
    }else{
      cmds <- paste0('cat ', f2merge, addition, newname)
      # cmds <- paste0('cat ', mysample, "*", unique(reads), suffix_final, addition, unique(newname))
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

fastq_group_names <- function(snames){
  y <- if(is.list(snames)) snames else strsplit(x = snames, split = "_")
  for(i in 1:max(sapply(y, length))) print(table(sapply(y, '[', i), useNA = 'always'))
}

fastq_fix_names <- function(x, str_replacements = NULL, sreplace = NULL, ...){
  newnames <- gsub("A", "A", x)
  # Modifications
  for(i in str_replacements) newnames <- gsub(i[1], i[2], newnames, ...)
  # Replacing
  if(!is.null(sreplace)){
    tvar <- newnames %in% names(sreplace)
    if(sum(tvar) != length(sreplace)) warning("Check the very special cases")
    cat(c(newnames[tvar], sum(tvar)), sep = "\n")
    newnames[tvar] <- unname(sreplace[newnames[tvar]])
  }

  newnames
}

fastq_make_metadata <- function(
  snames,
  original_names = NULL,
  name2column
){
  snames <- if(is.list(snames)) snames else strsplit(x = snames, split = "_")
  complete_names <- sapply(snames, paste0, collapse = "_")
  mytab <- data.frame(
    Sample_ID = if(is.null(original_names)) complete_names else original_names
  )
  if(!is.null(original_names)) mytab$Corrected_Sample_Name <- complete_names
  addtab <- sapply(X = name2column, FUN = function(x){
    y <- if(length(x) > 1){
      sapply(lapply(snames, '[', x), paste0, collapse = "-")
    }else{
      sapply(snames, '[', x)
    }
  })
  mytab <- cbind(mytab, addtab)
  return(mytab)
}

fastq_path_swap <- function(stringy, subs = NULL){
  if(is.null(subs)) return(stringy)
  if(!is.list(subs)) subs <- list(subs)
  for(i in subs){
    stringy <- gsub(i[1], i[2], stringy)
  }
  return(stringy)
}

fastq_reorder_checksum <- function(x, sortit = TRUE){
  w <- unique(sapply(x, function(y){
    z <- unlist(strsplit(y, " ")); paste0(rev(z[z != ""]), collapse = " ")
  }))
  if(sortit) return(sort(w))
  w
}
