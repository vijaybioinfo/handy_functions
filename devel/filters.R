#!/usr/bin/R

# Filter a table — adding columns... maybe complex designs
# If it's a file, you can specify it be added with 'append_file'
# You can pass 'expression instructions' where'expression ' is necessary
filters_complex <- meta_filtering <- function(
  mdata,
  filters = "none",
  cname = "none",
  sepchar = "~",
  cts = NULL,
  verbose = FALSE
) {
  if(is.null(filters[[1]][1])) filters <- "none"
  if(is.null(cname)) cname <- "none"
  if(filters[[1]][1] == "none" && cname[1] == "none") return(list(annotation = mdata))
  if(is.vector(mdata)) mdata <- data.frame(feature = mdata, row.names = mdata, stringsAsFactors = FALSE)
  if(is.matrix(mdata)) mdata <- data.frame(feature = rownames(mdata), row.names = rownames(mdata), stringsAsFactors = FALSE)
  if(any(filters %in% rownames(mdata))){
    filters <- show_found(x = filters, y = rownames(mdata), verbose = verbose)
    return(list(annotation = mdata[filters, ]))
  }
  tvar <- cname %in% colnames(mdata);
  if(sum(!tvar)) warning("Column name(s) not found: ", show_commas(cname[!tvar]))
  cname <- cname[tvar]
  if(verbose) cat("Watch-out columns:", show_commas(cname), "\n")
  if(verbose){ cat('Using filter:\n'); str(filters) }
  if(sum(grepl(sepchar, filters)) || length(filters) > 1 || is.list(filters)){
    filters <- filters_pattern(filters)
  }

  # subsetting annotation first checking for a file
  if(file.exists(filters[[1]][1])){
    if(verbose) cat("Filters from file", filters[[1]][1], "\n")
    filters <- c(filters, filters[[1]][1]) # separating the file from the rest of the filters
    filters[[1]] <- filters[[1]][-1]; filters <- rev(filters)
    filecon <- readfile(filters[[1]][1], stringsAsFactors = FALSE, verbose = verbose); if(verbose) str(filecon)
    tvar <- head(which(sapply(filecon, function(x) any(x %in% rownames(mdata)) )), 1)
    if(length(tvar) == 0) warning("No column in ", basename(filters[[1]][1]), " was compatible with 'mdata'")
    rownames(filecon) <- filecon[, tvar]
    tvar <- intersect(x = rownames(filecon), y = rownames(mdata))
    if(verbose) cat(length(tvar), "intersecting\n")
    mdata <- mdata[tvar, ]
    if("append_file" %in% unlist(filters)) mdata <- joindf(mdata, filecon)
    filterssuffix <- paste0("_", basename(filters[[1]][1]))
    if(length(filters) > 1) filters <- filters[-1]
  }

  # Add gene tag
  for(addgene in unique(c(cname, sapply(filters, head, 1)))){
    if(addgene %in% colnames(mdata)) next # if it's not in the annotation already
    gg <- unlist(strsplit(gsub("tag_|_[0-9]+", "", addgene), sepchar))
    if(grepl("tag", addgene) && all(gg %in% rownames(cts))){
      if(verbose) cat("Adding gene tags\n")
      tvar <- as.numeric(gsub("tag_.*_", "", unlist(strsplit(addgene, sepchar))))
      tvar <- ifelse(is.na(tvar), 0, tvar)
      tags_df <- features_add_tag(
        lgenes = gg,
        annot = mdata,
        mat = cts,
        thresh = tvar, tag = c('tag', 'p', 'n'),
        verbose = verbose
      )
      colname <- addgene
      tags_df[, colname] <- apply(tags_df, 1, paste, collapse = "")
      print(reshape2::melt(table(tags_df[, colname])))
      mdata <- cbind(tags_df, mdata); headmat(mdata)
    }
  }
  allcnames <- unique(c(cname, sapply(filters, head, 1)))
  allcnames <- allcnames[grepl(sepchar, allcnames)]
  allcnames <- allcnames[!allcnames %in% colnames(mdata)]
  if(length(allcnames) > 0){
    for(cname_i in allcnames){
      addeds <- unlist(strsplit(cname_i, sepchar))
      if(verbose) cat("Combining columns:", addeds, sep = "\n")
      mdata$combn <- apply(mdata[, addeds, drop = FALSE], 1, paste, collapse = "_")
      colnames(mdata) <- sub("^combn$", cname_i, colnames(mdata))
    }
  }

  # now if there's further filtering
  # watch out for 'cname' column samples/cells by table(cname, column_filtering_by)
  if(verbose) cat("Filtering now\n")
  if(any(sapply(filters, head, 1) %in% colnames(mdata))){
    if(verbose) cat("Based on a list\n");
    if(!is.null(cname) && verbose){ # watch out before
      tvar <- filters[!sapply(filters, head, 1) %in% cname[1]]
      void <- lapply(tvar, function(x) table(mdata[, c(cname[1], x[1])], useNA = 'always') )
      void <- reshape2::melt(void)
      if('value' %in% colnames(void)) print(head(void[void$value > 0, ], 20)) else print(head(void, 20))
    }
    mdata <- mdata[filters_subset_df(filters, mdata, verbose = verbose), ]
    if(!is.null(cname) && verbose){ # watch out after
      tvar <- filters[!sapply(filters, head, 1) %in% cname[1]]
      void <- lapply(tvar, function(x) table(mdata[, c(cname[1], x[1])], useNA = 'always') )
      void <- reshape2::melt(void)
      if('value' %in% colnames(void)) print(head(void[void$value > 0, ], 20)) else print(head(void, 20))
    }
    tvar <- sapply(filters, function(x) paste0(paste0(x[-1], collapse = "AND"), "_from_", x[1]) ) # creates "grp1ANDgrp2"
    filters <- paste0("_", paste0(tvar, collapse = "_and_")) # creates "grp1ANDgrp2_from_COLUMN1_and_grp1_from_COLUMN1"
    if(exists('filterssuffix')){
      filters <- paste0(filterssuffix, filters); rm(filterssuffix)
    }
  }

  tvar <- sapply(filters, function(x) any(grepl("^express", x)) )
  if(any(tvar)){
    if(verbose) cat("Based on an expression\n");
    filters <- sub("expr[A-z]+ ", "", filters[tvar][[1]])
    sset <- paste0("cellsf <- rownames(subset(annottab, subset = ", filters, "))")
    if(verbose) cat("Expression:", sset, "\n")
    eval(expr = parse(text = sset))
    mdata <- mdata[cellsf, ]
  }
  return(list(annotation = mdata, filter = filters))
}

# transform string/vector/list to list for subset
filters_pattern <- translist <- function(pat){
  if(is.null(pat)) return(character(0))
  if(grepl(":|;", pat[[1]][1])){
    pat <- strsplit(unlist(strsplit(pat, ";")), ":")
    pat <- lapply(pat, function(x) unlist(strsplit(gsub(x, pattern = " ", replacement = ""), ",")) )
    # return(pat)
  }
  if(sum(grepl("^list", pat))) pat <- eval(expr = parse(text = pat))
  if(!is.list(pat)) pat <- list(pat)
  # tvar <- any(sapply(pat, function(x) any(grepl('~', x)) && length(x) > 1 ))
  # if(tvar) return(pat) # if it has the separator and is a vector > 1, return it
  pat <- lapply(pat, function(x){
    if(grepl('c\\(.*)', x[1])){
      eval(expr = parse(text = x))
    }else if(grepl('~', x[1])){
      unlist(strsplit(x, '~'))
    }else{ return(x) }
  })
  pat
}
filters_pattern2list <- function(x){
  sapply(x, function(y){
    if(!grepl(":|;", y)) return(y)
    z <- gsub(";", "'), c\\('", gsub(":", "', '", y))
    paste0("\"list(c('", z, "'))\"")
  }, USE.NAMES = FALSE)
}

# subset data.frame with c('col','values') or a list(c('col1','values1'), c('col2','values2'))
filters_subset_df <- getsubset <- function(
  x = NULL,
  df,
  op = 'and',
  verbose = FALSE
){
  df <- remove.factors(df)
  if(is.null(x)) return(rownames(df))
  x <- filters_pattern(x)
  if(!is.null(x$op)){
    op = x[['op']]; x <- x[names(x) != "op"]
  }
  # if(is.null(x)) x <- c(colnames(df)[1], unique(df[, 1]))
  if(length(x) == 1 && class(x) != 'list') x <- c(x, unique(df[, 1]))
  x <- x[!is.na(x)]
  if(class(x) != 'list') x <- list(x)
  if(verbose){
    maxchar <- max(c(14, max(nchar(sapply(x, head, 1)))))
    tvar <- paste(addspaces('Category', maxchar), ':=>\t Size \t Classes\n')
    cat(tvar)
    cat(paste0(rep("-", nchar(tvar) + 2), collapse = ''), '\n')
  }
  x <- x[sapply(x, function(s){
    if(!s[1] %in% colnames(df)){ warning("No '", s[1], "' in ", show_commas(colnames(df))); FALSE }
    TRUE
  })]
  x <- lapply(x, function(s){
    stmp <- sub("^\\-", "", s[-1])
    tvar <- stmp %in% unique(df[, s[1]])
    if(!all(tvar)) warning("No '", show_commas(stmp[!tvar]), "' in ", s[1])
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
    if(verbose){
      cat(addspaces(s[1], maxchar), ':=>\t', sum(thisrows), '\t', show_commas(s[-1]), '\n')
    }
    return(thisrows)
  })
  if(op == 'and') tmp <- rowSums(tmp) == length(x) # all x elements are in the row
  if(op == 'or') tmp <- apply(tmp, 1, any) # only some elements are in the row
  tvar <- rownames(df[tmp, , drop = FALSE])
  if(length(tvar) == 0) warning('In total no selection')
  if(verbose && length(x) > 1) cat(addspaces("Total elements", maxchar), ':=>\t', length(tvar), '\t', length(x), '\n')
  return(tvar)
}

filters_summary <- summary_subset <- function(x){
  tvar <- sapply(x, function(x) paste0(paste0(x[-1], collapse = "AND"), "_from_", x[1]) )
  if(length(tvar) > 0) paste0(tvar, collapse = "_and_") else ""
}

# Create a column for cells expressing a gene
features_get_tag <- function(gg, sig = '+', prefix = FALSE){
  if(is.list(gg)){ cat("'gg' already a list\n"); return(gg) }
  lapply(gg, function(x){
      if(isTRUE(prefix)) return(c(paste0('tag_', x), paste0(sig, x)))
      c(paste0('tag_', x), paste0(x, sig))
  })
}
features_add_tag <- add_gene_tag <- function(
  lgenes,
  annot,
  mat,
  thresh = 0,
  tag = c('tag', '+', '-'),
  prefix = FALSE,
  verbose = FALSE
){
  if(length(thresh) != length(lgenes)){ # check if there is a threshold for each gene
    if(length(thresh) == 1) thresh <- rep(thresh, length(lgenes))
    if(!is.null(names(thresh))) thresh <- thresh[lgenes]
    thresh <- thresh[1:length(lgenes)]
    if(sum(is.na(thresh)))
      warning('No threshold specified for ', show_commas(lgenes[is.na(thresh)]), '; using 0\'s')
    thresh[is.na(thresh)] <- 0
  }
  names(thresh) <- lgenes
  newcolms <- cbindList(lapply(lgenes, function(thisgene){
    if(verbose) cat('----\nGene:', thisgene, '- threshold:', thresh[thisgene], '\n')
    if(prefix){
      tags <- paste0(tag, "_", casefold(thisgene, upper = TRUE))
    }else{
      tags <- paste0(tag[1], "_", casefold(thisgene, upper = TRUE))
      tags <- c(tags, paste0(casefold(thisgene, upper = TRUE), tag[-1]))
    }
    annot[, tags[1]] <- ifelse(as.vector(t(mat[thisgene, rownames(annot)] > thresh[thisgene])), tags[2], tags[3])
    tvar <- rownames(annot[annot[, tags[1]] == tags[3], ]) # check values
    tvar <- unlist(mat[thisgene, tvar])
    # if(verbose){ cat("Negative summary:\n"); print(summary(tvar)) }
    tvar <- rownames(annot[annot[, tags[1]] == tags[2], ]); tvar <- unlist(mat[thisgene, tvar])
    # if(verbose){ cat("Positive summary:\n"); print(summary(tvar)) }
    if(verbose){ cat('Proportions'); print(table(annot[, tags[1]])) }
    return(annot[, tags[1], drop = FALSE])
  }))
  newcolms
}

# when you have a vector of filters you want to place separated by
filters_further <- function(
  x,
  cnames = NULL,
  th = "",
  op = "&",
  verbose = FALSE
){
  if(is.list(x)){
    if(verbose) cat("It is a list of", length(x), "elements\n")
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
    if(verbose) cat("Digesting:", x, "\n")
    xcolumns <- paste0("(", paste0(cnames, collapse = "|"), ")")
    x <- gsub(xcolumns, paste0("x[, '", "\\1", "']"), x)
  }else{
    warning("Nothing to filter\n")
  }
  x
}

filters_thresholds <- getDEGenes <- function(
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
  verbose = FALSE
){
  x <- as.data.frame(x, stringsAsFactors = FALSE, check.names = FALSE)
  if(sum(gene_name %in% colnames(x))){
    rownames(x) <- x[, gene_name]
  }
  if(!is.null(catch)){
    pv <- max(x[catch, pvtype], na.rm = TRUE)
    # fc <- max(abs(x[catch, lfc.type]))
  }
  if(verbose) cat(pvtype,': ', pv, '\n', sep = "")
  if(verbose) cat(lfc.type,': ', fc, '\n', sep = "")

  assiggy <- "g <- rownames(x[which(x[, pvtype] <"
  mid <- "x[, lfc.type]"
  suf <- "fc"
  if(is.na(upreg)){
    if(verbose) cat('Both\n')
    express <- paste0(assiggy, th, "pv & abs(", mid, ") >", th, suf)
  }else if(isTRUE(upreg)){
    if(verbose) cat('Upregulated\n')
    express <- paste0(assiggy, th, "pv & ", mid, ">", th, suf)
  }else{
    if(verbose) cat('Downregulated\n')
    express <- paste0(assiggy, th, "pv & ", mid, "<", th, " -", suf)
  }
  if(!is.null(further)){ # Extra condition for which filter
    tvar <- filters_further(x = further, th = th, cnames = colnames(x), verbose = verbose)
    if(tvar != further) express <- paste(express, "&", tvar)
  }
  express <- paste0(express, "), ])")
  if(verbose) cat("------------\n", express, "\n------------\n")
  eval(expr = parse(text = express))
  n.gs <- length(g) # Get number retrieved genes
  if(verbose) cat('Retrieved genes:', show_commas(g), '\n')
  if(n.gs == 0 || is.null(g) || is.na(g)){
    return(character(0))
  }else{
    return(g)
  }
}

sample_even <- sample_grp <- function(
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
    cat("Given group:", show_commas(names(grsize)), "\n")
    cat("Init. sizes:", show_commas(grsize), "\n")
    cat("Final sizes:", show_commas(finalgrsize), "\n")
    cat("% from total:", (sum(finalgrsize) / nrow(annot)) * 100, "\n")
    cat("Returning:", length(scells), "\n")
  }
  scells
}

### Figure configuration ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig_config_template = list(
  name = "figure_1A",
  metadataf = "/path/to/file{.csv,.rdata,.tsv}",
  edataf = "none", # optional in case metadataf is an object containing this data
  sample_filter = NULL,
  features = NULL,
  type = "violin|dotplot|etc",
  axis_x = list(name = "column", order = c("group1", "group2")),
  axis_y = 'Label',
  facet = NULL,
  transf = NULL, # log2, log2cpm
  individually = TRUE, # if you want to plot each feature
  ncol = 2, size = c(5, 5)
)
fig_config_check <- function(config, template = fig_config_template){
  config_names <- names(fig_config_template)
  config_names <- config_names[!config_names %in% names(config)]
  for(config_name in config_names){
    config[[config_name]] <- NULL
  }
  if(!is.list(config$axis_x)) config$axis_x <- list(name = config$axis_x)
  return(config)
}
fig_replace <- function(config, replacement){
  for(i in 1:length(replacement)){
    config[[names(replacement)[i]]] <- replacement[[i]]
  }
  return(config)
}
fig_features <- function(config, features, verbose = TRUE){
  if(!is.null(config$features)) config <- list(config)
  one_config <- length(config) == 1
  if(verbose) cat(length(config), "configuration(s)\n")
  config <- lapply(X = config, FUN = function(x){
    cat(x$name, "\n")
    x$features <- show_found(x = x$features, features, verbose = verbose)
    cat("\n\n")
    x
  })
  if(one_config) config[[1]] else config
}

fig_set_data <- function(
  config,
  keep_object = TRUE, # takes too much space
  verbose = TRUE
){
  # Digesting
  config <- fig_config_check(config = config)
  config$edataf <- if(is.null(config$edataf)) "none" else config$edataf
  config$edatafbk <- if(is.null(config$edatafbk)) "none" else config$edatafbk
  config$metadatafbk <- if(is.null(config$metadatafbk)) "file_name" else config$metadatafbk

  cat("# Fetching data # ------------------------------------------------------------\n")
  config$metadata <- if(isTRUE(config$metadataf != config$metadatafbk)){
    if(verbose) cat("Loading [meta]data\n", config$metadataf, "\n")
    readfile(config$metadataf, row.names = 1, stringsAsFactors = FALSE, verbose = verbose)
  }else{ config$metadata }
  config$edata <- if(isTRUE(config$edataf != config$edatafbk)){
    if(verbose) cat("Loading features\n", config$edataf, "\n")
    readfile(config$edataf)
  }else{ config$edata }

  myobject <- casefold(sapply(config, class))
  config$odata <- if(any(myobject %in% "seurat")){
    config[[which(myobject %in% "seurat")[1]]]
  }

  if(!is.null(config$odata)){ # Seurat
    if(verbose) cat("Using", class(config$odata), "object\n")
    tvar <- as.numeric(gsub("(..).*", "\\1", config$odata@version))
    if(tvar < 3) config$odata <- suppressMessages(UpdateSeuratObject(config$odata))
    if(is.null(config$metadata) || !is.data.frame(config$metadata)) config$metadata <- config$odata@meta.data
    tvar <- is.null(config$edata) || length(dim(config$edata)) != 2 || config$edataf == "none"
    if(tvar) config$edata <- config$odata@assays$RNA@data
  }
  config$metadatafbk <- config$metadataf
  config$edatafbk <- config$edataf

  ## Operations
  if(verbose){
    tvar <- unique(c(head(1:ncol(config$metadata)), tail(1:ncol(config$metadata))))
    cat("Metadata\n"); str(config$metadata[, tvar])
    cat("Features data\n"); str(config$edata)
    cat("Config\n"); str(config, max.level = 1)
  }

  # Check it it'll be a folder or just a prefix # ------------------------------
  config$name <- dircheck(config$name)
  if(verbose) cat(config$name, "\n")
  rownames(config$edata) <- features_parse_ensembl(rownames(config$edata))
  myfeatures <- if(!is.null(config$features)) features_parse_ensembl(config$features)
  cat("# Filtering features and agents (samples/cells) # ----------------------------\n")
  myfeatures <- rownames(filters_complex(
    mdata = config$edata,
    filters = myfeatures,
    verbose = verbose
  )$annotation)
  # It can also add a new category based on feature levels
  tvar <- if(is.null(config$axis_x$name)){
    unique(unlist(sapply(config$axis_x, function(x) x$name )))
  }else{ config$axis_x$name }
  tvar <- filters_complex(
    mdata = config$metadata,
    filters = config$sample_filter,
    cname = tvar,
    cts = config$edata,
    verbose = verbose
  )
  ssamples <- rownames(tvar$annotation)
  edata_ss <- if(!is.null(config$transf)){
    cat("# Applying transformations # -------------------------------------------------\n")
    tvar <- c(gsub("log[0-9]{,1}", "", config$transf), gsub("(log[0-9]{,1}).*", "\\1", config$transf))
    edata_ss <- count_transformation(
      cts = config$edata[rownames(config$edata) %in% myfeatures, colnames(config$edata) %in% ssamples],
      transf = tvar[1], verbose = verbose
    )
    count_transformation(cts = edata_ss[[1]][myfeatures, ], transf = tvar[2], verbose = verbose)[[1]]
  }else{ config$edata[rownames(config$edata) %in% myfeatures, colnames(config$edata) %in% ssamples] }
  cat("# Sorting the identities # ---------------------------------------------------\n")
  ddf <- config$metadata[ssamples, ]
  if(length(myfeatures) != nrow(config$edata)) ddf <- cbind(ddf, t(as.matrix(edata_ss))[ssamples, ])
  ddf <- fig_set_identities(mdata = ddf, idents = config$axis_x, verbose = verbose)

  config$pdata = ddf
  config$samples = ssamples
  config$features = myfeatures
  if(!keep_object) config <- config[!names(config) %in% "odata"]

  return(config)
}

fig_set_identities <- function(mdata, idents = NULL, verbose = FALSE){
  if(is.null(idents)) return(mdata)
  idents <- if(is.null(idents$name)) idents else list(idents)
  if(is.null(names(idents))){
    new_cnames <- paste0("Identity", c("", 1:20));
    tvar <- new_cnames %in% colnames(mdata)
    if(sum(tvar)) warning(show_commas(new_cnames[tvar]), " exists in data")
    names(idents) <- new_cnames[1:length(idents)]
  }
  mdata$Identity <- NULL
  for(i in names(idents)){
    if(verbose) cat("-- Column(s):", show_commas(idents[[i]]$name), "\n")
    tvar <- idents[[i]]$name %in% colnames(mdata)
    if(any(!tvar)){
      warning("Column(s) not found: ", show_commas(idents[[i]]$name[!tvar]))
      warning("Possible values: ", show_commas(colnames(mdata), Inf))
    }
    mdata$tmp <- do.call(paste, c(remove.factors(mdata[, idents[[i]]$name, drop = FALSE]), sep = "_"))
    mdata$tmp <- if(!is.null(idents[[i]]$identnames)){
      if(verbose) cat("Replacing identities\n")
      factor(idents[[i]]$identnames[as.character(mdata$tmp)], levels = unique(unname(idents[[i]]$identnames)))
    }else{ mdata$tmp }
    mdata$tmp <- if("REST" %in% idents[[i]]$order){
      tvar <- as.character(mdata$tmp)
      ifelse(!tvar %in% idents[[i]]$order, "REST", tvar)
    }else{ mdata$tmp }
    mdata$tmp <- if(!is.null(idents[[i]]$order) && any(mdata$tmp %in% idents[[i]]$order)){
      factor(mdata$tmp, idents[[i]]$order)
    }else{ factormix(mdata$tmp) }
    if(verbose) cat("Order:", show_commas(levels(mdata$tmp)), "\n")
    colnames(mdata) <- gsub("^tmp$", i, colnames(mdata))
  }
  return(mdata)
}
