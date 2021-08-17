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
  if(length(filters) == 0) filters <- "none"
  if(is.null(filters[[1]][1])) filters <- "none"
  if(is.null(cname)) cname <- "none"
  if(filters[[1]][1] == "none" && cname[1] == "none") return(list(annotation = mdata))
  if(is.vector(mdata))
    mdata <- data.frame(feature = mdata, row.names = mdata, stringsAsFactors = FALSE)
  if(is.matrix(mdata))
    mdata <- data.frame(mdata, row.names = rownames(mdata), stringsAsFactors = FALSE)
  if(any(filters %in% rownames(mdata))){
    filters <- show_found(x = filters, y = rownames(mdata), verbose = verbose)
    return(list(annotation = mdata[filters, ]))
  }
  if(verbose) cat("Watch-out columns:", show_commas(cname), "\n")
  if(verbose){ cat('Using filter:\n'); str(filters) }
  if(sum(grepl(sepchar, filters)) || length(filters) > 1 || is.list(filters)){
    filters <- filters_pattern(filters)
  }

  ## Checking for a file ## ----------------------------------------------------
  if(file.exists(filters[[1]][1])){
    if(verbose) cat("Filters from file", filters[[1]][1], "\n")
    filters <- c(filters, filters[[1]][1]) # separating the file from the rest of the filters
    filters[[1]] <- filters[[1]][-1]; filters <- rev(filters)
    filecon <- readfile(filters[[1]][1], stringsAsFactors = FALSE, verbose = verbose)
    if(verbose) str(filecon)
    tvar <- head(which(sapply(filecon, function(x) any(x %in% rownames(mdata)) )), 1)
    if(length(tvar) == 0)
      warning("No column in ", basename(filters[[1]][1]), " was compatible with 'mdata'")
    rownames(filecon) <- filecon[, tvar]
    tvar <- intersect(x = rownames(filecon), y = rownames(mdata))
    if(verbose) cat(length(tvar), "intersecting\n")
    mdata <- mdata[tvar, ]
    if("append_file" %in% unlist(filters)) mdata <- joindf(mdata, filecon)
    filterssuffix <- paste0("_", basename(filters[[1]][1]))
    if(length(filters) > 1) filters <- filters[-1]
  }

  ## Add gene tag ## -----------------------------------------------------------
  # Example: tag_GENE or tag_GENE_10 > GENEp and GENEn
  # tag_CD4~tag_CD8B > CD4p_CD8Bp, CD4p_CD8Bn, CD4n_CD8Bp, CD4n_CD8Bn
  for(addgene in unique(c(cname, sapply(filters, head, 1)))){
    if(addgene %in% colnames(mdata)) next # if it's not in the annotation already
    gg <- unlist(strsplit(gsub("tag_|_[0-9]+", "", addgene), sepchar))
    if(grepl("tag", addgene) && all(gg %in% rownames(cts))){
      if(verbose) cat("Adding gene tags\n")
      explevel <- as.numeric(gsub("tag_.*_", "", unlist(strsplit(addgene, sepchar))))
      explevel <- ifelse(is.na(explevel), 0, explevel)
      tags_df <- features_add_tag(
        gg, mdata, cts,
        thresh = tvar, tag = c('tag', 'p', 'n'),
        verbose = verbose
      )
      tags_df[, addgene] <- apply(tags_df, 1, paste, collapse = "_")
      print(reshape2::melt(table(tags_df[, addgene])))
      mdata <- cbind(tags_df, mdata); headmat(mdata)
    }
  }

  ## Combined columns ## -------------------------------------------------------
  # Example: columnA: elementA1, elementA2
  # columnB: elementB5, elementB8
  # columnA~columnB, elementA1_elementB5
  # disease~stim, Mild_Unstim (if you are comparing Mild_Unstim vs Sever_Unstim)
  allcolumns <- unique(c(cname, sapply(filters, head, 1)))
  allcolumns <- allcolumns[grepl(sepchar, allcolumns)]
  allcolumns <- allcolumns[!allcolumns %in% colnames(mdata)]
  if(length(allcolumns) > 0){
    if(verbose) cat("Combining columns\n")
    for(cname_i in allcolumns){
      addeds <- unlist(strsplit(cname_i, sepchar))
      if(verbose) cat(" -", addeds, "\n")
      mdata$combn <- apply(mdata[, addeds, drop = FALSE], 1, paste, collapse = "_")
      colnames(mdata) <- sub("^combn$", cname_i, colnames(mdata))
    }
  }

  ## Based on a list ## --------------------------------------------------------
  # watch out for 'cname' column samples/cells by table(cname, column_filtering_by)
  if(verbose) cat("Filtering:")
  did_filter <- nrow(mdata); filters0 <- filters
  filters_found <- sapply(filters, head, 1) %in% colnames(mdata)
  if(any(filters_found)){
    if(any(!grepl("^expr", sapply(filters[!filters_found], head, 1)))){
      stop("Missing column ", sapply(filters[!filters_found], head, 1))
    }
    filters <- filters[filters_found]
    if(verbose) cat("- based on a list\n")
    if(isTRUE(cname %in% colnames(mdata)) && verbose && (length(filters) > 1)){ # watch out before
      tvar <- filters[!sapply(filters, head, 1) %in% cname[1]]
      void <- lapply(tvar, function(x) table(mdata[, c(cname[1], x[1])], useNA = 'always') )
      void <- reshape2::melt(void)
      if('value' %in% colnames(void)) print(head(void[void$value > 0, ], 20)) else print(head(void, 20))
    }
    mdata <- mdata[filters_subset_df(filters, mdata, verbose = verbose), ]
    if(isTRUE(cname %in% colnames(mdata)) && verbose && (length(filters) > 1)){ # watch out after
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
  }; filters <- filters0

  tvar <- sapply(filters, function(x) any(grepl("^expr", x)) )
  if(any(tvar)){
    if(verbose) cat("- based on an expression\n")
    filters <- unlist(filters[tvar])
    myfilters <- filters[!grepl("^expr", filters)];
    for(myfilter in myfilters){
      sset <- paste0("cellsf <- rownames(subset(mdata, subset = ", myfilter, "))")
      if(verbose) cat("Expression:", sset, "\n")
      eval(expr = parse(text = sset)); mdata <- mdata[cellsf, ]
    }
  }
  tvar <- cname %in% colnames(mdata);
  if(sum(!tvar) > 0) cat("Column name(s) not found: ", show_commas(cname[!tvar]), "\n")
  if(did_filter == nrow(mdata) && verbose){
    cat("- no filter applied\n"); filters = ""
  }
  return(list(annotation = mdata, filter = filters))
}

# transform string/vector/list to list for subset
# column~element
# column:element;column2:element3
# Column == 'Element1'
# "list(c('column', 'element'), c('column2', 'element3'))"
# "list(c('column', '-element'), c('column2', 'element3'))"
# "list(c('RNA_snn_res.0.6', '4', '6'), c('orig.celltype', 'CD4'))"
# "expr: RNA_snn_res.0.6 %in% c('4', '6') & orig.celltype == 'CD4'"
filters_pattern <- translist <- function(pat){
  if(is.null(pat)) return(character(0))
  if(grepl(":|;", pat[[1]][1])){
    pat <- strsplit(unlist(strsplit(pat, ";")), ":")
    pat <- lapply(pat, function(x) unlist(strsplit(gsub(x, pattern = " ", replacement = ""), ",")) )
  }
  if(sum(grepl("^list", pat))) pat <- eval(expr = parse(text = pat))
  if(!is.list(pat)) pat <- list(pat)
  pat <- lapply(pat, function(x){
    if(any(grepl('~', x)) && length(x) > 1) return(x)
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
  if(isTRUE(x == "")) return("")
  tvar <- sapply(x, function(x) paste0(paste0(x[-1], collapse = "AND"), "_from_", x[1]) )
  if(length(tvar) > 0) paste0(tvar, collapse = "_and_") else ""
}

filters_vector = function(
  x,
  include = "", # "" means everything
  exclude = "noneofthem", # "noneofthem" so nothing matches
  rename = NULL,
  verbose = TRUE
){
  obs_in = if(!is.null(include)){
    if(length(include) > 1) x[x %in% include] else grep(pattern = include, x = x, value = TRUE)
  }else{ x }
  obs_ex = if(!is.null(exclude)){
    if(length(exclude) > 1) x[x %in% exclude] else grep(pattern = exclude, x = x, value = TRUE)
  }else{ exclude }
  if(verbose){ cat("Excluding:\n"); print(format(unname(obs_ex), justify = "centre")) }
  obs_in = obs_in[!obs_in %in% obs_ex]; names(obs_in) <- obs_in
  for(i in names(rename)) names(obs_in) <- gsub(i, rename[[i]], names(obs_in))
  if(verbose){ cat("Final names:\n"); print(format(obs_in, justify = "centre")) }
  return(obs_in)
}

# get columns with a max N of groups
filters_columns <- function(
  mdata,
  onames = NULL,
  maxn = 50,
  types = "all",
  na_rm = FALSE,
  duplicate_rm = TRUE,
  verbose = TRUE,
  ... # arguments for filters_vector
){
  if(is.null(onames)) onames <- colnames(mdata); onames <- rev(onames)
  if(any(c("numeric", "integer") %in% types)) maxn <- Inf
  tvar <- sapply(mdata[, onames, drop = FALSE], class)
  if("all" %in% types) types <- unique(tvar)
  if(verbose) cat("Types:", show_commas(types), '\n')
  tvar <- setNames(tvar %in% types, names(tvar)) # if it's the type and has no NA's
  if(isTRUE(na_rm)) tvar <- tvar & !sapply(mdata[, names(tvar), drop = FALSE], function(x) any(is.na(x)) )
  onames <- onames[onames %in% names(which(tvar))]
  onames <- filters_vector(#include = include, exclude = exclude, rename = rename,
    x = onames, verbose = verbose > 1, ...
  )
  if(length(onames) == 0){ if(verbose) cat("Taking all columns\n"); onames <- colnames(mdata) }
  if(verbose) cat("Filtering", length(onames), "\n")
  nnames <- sapply(setNames(onames, unname(onames)), function(x){
    if(is.numeric(mdata[, x])) -1 else length(unique(mdata[, x]))
  })
  onames <- onames[order(nnames)]
  onames <- onames[nnames[onames] != nrow(mdata)] # excluding  N = total rows
  onames <- onames[nnames[onames] > 1] # excluding N = 1
  onames <- onames[nnames[onames] <= maxn] # N <= 27 because ggpairs doesn't allow more than this
  onames <- onames[!is.na(onames)]
  if(length(onames) == 1) return(onames)
  if(duplicate_rm & length(onames) > 2){
    tvar <- is.finite(nnames) & nnames > 1;
    dnames <- onames[tvar]
    if(verbose) cat("Columns wiht same N:", length(dnames), "\n")
    if(verbose) cat("Checking if groups in columns are the same\n")
    if(verbose > 1) print(reshape2::melt(sort(nnames[unname(dnames)])))
    pairs <- try(gtools::combinations(n = length(dnames), r = 2, v = dnames), silent = TRUE)
    if(class(pairs) != 'try-error'){
      kha <- apply(pairs, 1, function(x){
        y <- table(mdata[, c(x)])
        diag(y) <- 0; any(rowSums(y) > 0) # take when none but diagonal are 0, aka, same N groups
      })
      if(verbose && any(!kha)) cat(sum(!kha), "matching 1-to-1 groups:", show_commas(pairs[!kha, 2]), "\n")
      if(any(!kha)){
        dnames <- dnames[!dnames %in% pairs[!kha, 2]]
        onames <- c(onames[!onames %in% dnames], dnames)
      }else{ cat("Keeping all explored\n") }
    }else if(verbose){
      cat("n =", length(onames), "/ unique =", length(unique(onames)), "-", show_commas(unique(onames)), "\n")
    }# by row, is the element duplicated and then are all elements in a given row duplicated?
  }
  kha <- apply(t(apply(mdata[, onames, drop = FALSE], 1, duplicated)), 2, sum)
  onames <- onames[(kha != nrow(mdata))]
  if(verbose) cat("Returning:", length(onames), "\n")
  onames
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
  # Used cbindList > list2evendf > mat_names
  newcolms <- lapply(lgenes, function(thisgene){
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
  })
  as.data.frame(newcolms)
}

features_find = function(
  features,
  universe,
  thesarus = TRUE,
  verbose = FALSE
) {
  features_f <- features[file.exists(features)]
  features_df = if(length(features_f)){
    if(verbose) cat("From file(s):", features_f, sep = "\n")
    tmp = lapply(X = features_f, FUN = read.csv, stringsAsFactors = FALSE)
    if(length(features) > 1){
      y <- features[!features %in% features_f]
      if(verbose) cat("Extra", length(y), "features\n")
      tmp <- c(tmp, data.frame(gene_name = y, stringsAsFactors = FALSE))
    }; data.table::rbindlist(tmp, fill = TRUE)
  }else{
    data.frame(gene_name = features, stringsAsFactors = FALSE)
  }
  features_df$found <- apply(X = features_df, MARGIN = 1, FUN = function(x){
    y <- unlist(strsplit(x, "~")) # if() y <- casefold(y, upper = TRUE)
    z <- y[y %in% universe]
    if(length(z) == 0){
      z <- unique(c(y, unlist(Seurat::GeneSymbolThesarus(
        symbols = y, timeout = 15, several.ok = TRUE))))
      z <- y[y %in% universe]
    }; ifelse(length(z) == 0, "", z)
  })
  if(any(features_df$found == "") && verbose){
    cat("Not found\n"); str(features_df[features_df$found == "", ])
  }
  return(features_df)
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
  verbose = FALSE
){
  annot$tmp123 = make.names(as.character(annot[, cname]))
  grsize <- table(annot$tmp123[!is.na(annot[, cname])])
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
  if(verbose){ cat("Factors:"); str(factored)}
  scells <- lapply(names(grsize), function(x){
    sample(rownames(annot[which(annot$tmp123 == x), ]), factored[x])
  }); names(scells) <- names(grsize)
  finalgrsize <- sapply(scells, length)
  scells <- unlist(scells)
  if(verbose){
    cat("Given group:", show_commas(names(grsize)), "\n")
    cat("Init. sizes:", show_commas(grsize), "\n")
    cat("Final sizes:", show_commas(finalgrsize), "\n")
    cat("% from total:", (sum(finalgrsize) / nrow(annot)) * 100, "\n")
    cat("Returning:", length(scells), "\n")
  }
  scells
}

## Utilities
remove.factors <- function (df) {
  for (varnum in 1:length(df)) {
    if ("factor" %in% class(df[, varnum])) {
      df[varnum] = as.character(df[, varnum])
    }
  }
  return(df)
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

addspaces <- function(x, m) paste0(x, paste0(rep(' ', m-nchar(x)), collapse = ''))
