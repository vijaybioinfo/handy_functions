#!/usr/bin/R

###################
# Merge meta data #
###################

# Mergin metadata tables and prepearing expession data

source("/home/ciro/scripts/handy_functions/devel/overlap.R")

rnaseq_merge_metadata <- function(
  fnames,
  rbinds = NULL, # 0 and 1, same length as fnames
  outdir = NULL,
  mytag = "",
  subset_if_present = NULL, # get these colummns only if present in a table
  rename_columns = list(), # list(c('name', 'newname'))
  replace_in_col = list(), # list(c('column', 'string', 'replacement'))
  exclude = NULL,
  check_names = c('Major_Cell_type', 'Cell_type', 'Disease', 'Treatment', 'Stim'),
  verbose = TRUE
){
  outdir <- if(!grepl("\\/$", outdir) && !is.null(outdir)) paste0(outdir, "/") else outdir

  if(is.null(names(fnames))) names(fnames) <- paste0("T", 1:length(fnames))
  mytables <- if(is.list(fnames)){
    fnames
  }else{
    if(verbose) cat("Reading metadata tables:", fnames, sep = "\n")
    lapply(fnames, read.csv, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
  }
  if(verbose){ cat("Sizes:\n"); print(sapply(mytables, dim)) }

  ## Processing
  if(verbose) cat("Processing tables\n")
  mytables <- lapply(names(mytables), function(x){
    mytab <- mytables[[x]]
    if(verbose) cat(" *", x, "\n")
    if(any(subset_if_present %in% colnames(mytab))){
      if(verbose) cat("  Subsetting\n")
      mytab <- mytab[, colnames(mytab) %in% subset_if_present, drop = FALSE]
    }
    if(any(exclude %in% colnames(mytab))){
      if(verbose) cat("  Excluding\n")
      mytab <- mytab[, !colnames(mytab) %in% exclude, drop = FALSE]
    }; tvar <- sapply(rename_columns, head, 1) %in% colnames(mytab)
    if(any(tvar)){
      if(verbose) cat("  Renaming\n")
      for(i in rename_columns[tvar]){
        colnames(mytab)[colnames(mytab) %in% i[1]] <- i[2]
      }
    }; tvar <- sapply(replace_in_col, head, 1) %in% colnames(mytab)
    if(any(tvar)){
      if(verbose) cat("  Replace in columns\n")
      for(i in replace_in_col[tvar]){
        mytab[, i[1]] <- gsub(i[2], i[3], mytab[, i[1]])
      }
    }
    # rownames(mytab) <- gsub(" {1,}", "", rownames(mytab))
    mytab
  }); names(mytables) <- names(fnames)
  if(!is.null(rbinds)){
    tvar <- data.frame(data.table::rbindlist(mytables[rbinds == 1], fill = TRUE), stringsAsFactors = FALSE, check.names = TRUE)
    rownames(tvar) <- unlist(lapply(mytables[rbinds == 1], rownames))
  }
  is_mergeable <- sapply(1:length(mytables), function(x){
    any(rownames(mytables[[x]]) %in% unlist(sapply(mytables[-x], rownames)))
  }) # getting mergeable table; matching rownames
  if(length(mytables) == 1){
    is_mergeable <- TRUE
  }
  if(length(mytables) == 2 && sum(is_mergeable) != 2){
    is_mergeable <- sapply(mytables, nrow) == max(sapply(mytables, nrow))
    mytables[['mock']] <- mytables[is_mergeable][[1]]; is_mergeable[3] <- TRUE
  }
  if(verbose){
    cat(paste0(c("Mergeable:", basename(fnames[is_mergeable])), "\n"))
    cat(paste0(c("For expansion:", basename(fnames[!is_mergeable])), "\n"))
  }
  mergeable <- mytables[is_mergeable]
  expanding <- mytables[!is_mergeable] # getting possible expanding ones

  ## Planned: Generate further meta data information from sample names

  mergeable <- mergeable[names(mergeable) != ""]
  if(verbose) print(unique(sapply(mergeable, colnames)))

  if(length(expanding)){
    if(verbose) cat("Finding expanding columns\n")
    expanding_column <- lapply(expanding, function(x){
      lapply(mergeable, function(y){ # check in mergeables
        column_eval <- lapply(y, function(z){ # per column
          any(z %in% rownames(x)) # any matching? Ideally all
        })
        column_eval
      })
    })
    column_match <- try(reshape2::melt(expanding_column), silent = TRUE)
    if(class(column_match) != 'try-error'){
      column_match <- column_match[column_match[, 1], ]
      if(nrow(column_match) > 0){
        if(verbose) cat("Expanding\n")
        for(i in 1:nrow(column_match)){
          matching <- unname(unlist(column_match[i, ]))[-1]
          if(verbose) cat(paste0(head(matching), collapse = ", "), "\n") # Present
          tvar <- mergeable[[matching[2]]][, matching[1]]
          table_expand <- expanding[[matching[3]]] # select table
          tmp <- tvar[tvar %in% rownames(table_expand)]
          table_expand <- table_expand[tvar, ] # expand with matching column
          rownames(table_expand) <- rownames(mergeable[[matching[2]]]) # rename expansion
          mergeable[[matching[2]]] <- dplyr::left_join(mergeable[[matching[2]]], table_expand) # expand on matching table
        }
      }
    }
  }

  # Creating super meta data table
  tvar <- list(
    unique(unlist(lapply(mergeable, rownames))),
    unique(unlist(lapply(mergeable, colnames)))
  )
  final_annot <- data.frame(
    matrix(nrow = length(tvar[[1]]), ncol = length(tvar[[2]]), dimnames = tvar),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  # Filling meta data table
  for(i in 1:length(mergeable)){
    final_annot[rownames(mergeable[[i]]), colnames(mergeable[[i]])] <- mergeable[[i]]
  }

  final_metadata <- data.frame(
    cbind(Sample_ID = rownames(final_annot), final_annot),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if(verbose) str(final_metadata)
  tvar <- check_names %in% colnames(final_metadata)
  if(verbose && sum(tvar) > 2){
    tvar <- gtools::combinations(n = sum(tvar), r = 2, v = check_names[tvar])
    cat("Showing groups:\n"); print(head(tvar))
    for(i in 1:nrow(tvar)){
      print(head(table(final_metadata[, tvar[i, ]]), 20))
    }
  }else if(verbose && any(tvar)){
    cat("Showing groups:", paste0(check_names[tvar], collapse = ", "), "\n");
    tvar <- sapply(final_metadata[, check_names[tvar], drop = FALSE], table)
    print(sapply(tvar, head, 20))
  }
  if(!is.null(outdir)){
    fname <- paste0(outdir, "metadata_", mytag, ".csv")
    if(verbose) cat("Writing to:", fname, "\n")
    write.csv(final_metadata, file = fname, row.names = FALSE)
  }

  ## Epilogue: treating missing ones... all table should ideally have the same number of samples
  head(final_metadata[, head(colnames(final_metadata))])
  ov_names_list <- lapply(mergeable, rownames); names(ov_names_list) <- NULL
  ov_names <- overlap_list(ov_names_list)
  ov_names <- ov_names[sapply(ov_names, length) > 0]
  if(length(ov_names) > 1){
    if(verbose) cat("Merged table is a subset\n")
    if(verbose) str(ov_names)
    filt_final_metadata <- final_metadata[ov_names[[which.max(sapply(ov_names, length))]], ]
    if(verbose) cat("Taking", nrow(filt_final_metadata), "samples\n")
    if(!is.null(outdir)){
      fname <- paste0(outdir, "metadata_filtered_", mytag, ".csv")
      if(verbose) cat("Writing to:", fname, "\n")
      write.csv(filt_final_metadata, file = fname, row.names = FALSE)
    }
  }

  return(final_metadata)
}


# Getting tables of expression
rnaseq_merge_edata <- function(
  run_dirs,
  mdata = NULL, # If you want to cross check how many samples overlap
  outdir = NULL,
  keep1st_seen_sample = NA,
  # TRUE, keeps the first seen sample names while column binding
  # FALSE, keeps the second
  # NA, append the folder names as prefix
  verbose = TRUE
){
  outdir <- if(!grepl("\\/$", outdir) && !is.null(outdir)) paste0(outdir, "/") else outdir
  edatas <- lapply(X = c('TPM', 'raw'), FUN = function(ctype){
    mats <- lapply(run_dirs, function(dname){
      dname <- paste0(dname, ifelse(grepl("\\/$", dname), "", "/"))
      matname <- paste0(dname, '4.Output/counts/', ctype, '_counts.csv')
      if(verbose) cat(matname, "\n")
      mat <- data.frame(data.table::fread(matname), check.names = FALSE, stringsAsFactors = FALSE)
      rownames(mat) <- mat[, 1]; mat <- as.matrix(mat[, -1, drop = FALSE])

      config <- rjson::fromJSON(paste(readLines(paste0(dname, "conf_RNA_Seq.json")), collapse=""))
      config <- sapply(config, function(x) if(is.list(x)) x else as.list(x) )
      gannotname <- config$config$annotation_file
      gannot <- read.csv(gannotname, row.names = 1, stringsAsFactors = FALSE)

      newnames <- paste0(gsub("'", '', gannot[rownames(mat), 'gene_name']), "_", rownames(mat))
      rownames(mat) <- newnames

      return(mat)
    })
    if(is.na(keep1st_seen_sample)){
      mynames <- unique(unlist(sapply(mats, colnames), use.names = FALSE))
      if(length(mynames) < sum(sapply(mats, ncol))){
        warning("Reapeated across runs\n")
        mats <- lapply(1:length(mats), function(x){
          mat <- mats[[x]]
          colnames(mat) <- paste0(names(mats)[x], ".", colnames(mat))
          return(mat)
        })
      }
    }
    names(mats) <- names(run_dirs)
    edata <- mats[[1]]
    if(length(mats) >= 2){
      for(i in 2:length(mats)){
        if(isTRUE(keep1st_seen_sample)){
          if(verbose) cat("Preserving first found samples\n")
          edata <- cbind(edata, mats[[i]][, !colnames(mats[[i]]) %in% colnames(edata)])
        }else if(is.na(keep1st_seen_sample)){
          if(verbose) cat("All samples kept\n")
          edata <- cbind(edata, mats[[i]])
        }else{
          if(verbose) cat("Updating samples\n")
          edata <- cbind(edata[, !colnames(edata) %in% colnames(mats[[i]])], mats[[i]])
        }
      }
    }else{ if(verbose) cat("Only one matrix\n") }
    if(verbose) cat("Dimensions", dim(edata), "\n")
    if(verbose && !is.null(mdata)){
      cat("Intersection with metadata:", length(intersect(colnames(edata), rownames(mdata))), "\n")
      cat("Intersection with filtered metadata:", length(intersect(colnames(edata), rownames(mdata))), "\n")
    }
    if(!is.null(outdir)){
      fname <- paste0(outdir, ctype, ".Rdata")
      if(verbose) cat("Writing to:", fname, "\n\n")
      save(edata, file = fname)
    }
    return(edata)
  })
}
