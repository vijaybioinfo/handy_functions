#!/usr/bin/R

###############################
# Handy functions for n-tasks #
###############################

# This functions are designed to be use across all scripts from Vijay's Lab
# NOTE: highly sensitive, dont change input parameters names so confidently

ll <- function(x = getwd()) system(paste('ls -hlt', x))

# get a name or a "blank" character ('')
newname <- function(x, default = '') ifelse(x != default, x, '')

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

# Number of character occurences in a string
countChars <- function(char, string) {
    tvar <- gsub(char, "", string)
    return (nchar(string) - nchar(tvar))
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

# found in column
colfound <- function(x, tab) sum(sapply(x, function(y) grepl(paste0('^', y, '$'), colnames(tab))))

found_partial <- function(x, y){
  # tvar <- paste0(x, collapse = "|")
  z <- unlist(sapply(x, function(tvar) y[grepl(tvar, y)] ), use.names = FALSE)
  return(z)
}

# set rownames from column
setrows <- function(tab, x = 'gene_name', v = FALSE){
  x <- show_found(x, colnames(tab))
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
    gg <- show_found(genes_use[[name]], rownames(object@assays$RNA@data), element = 'Genes', v = v)
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
    gg <- show_found(thesegenes, rownames(object@assays$RNA@data), element = 'Genes', v = v)
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

# midpoints
mids <- function(x) c(x[1]/2, x[-length(x)] + diff(x)/2)

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
  gg <- show_found(x, sour, v = v)
  if(!length(gg)){
    if(v) cat('Sampling... ')
    set.seed(myseed); gg <- sample(sour, ln)
    if(v) cat('Got:', commas(gg), '\n')
  }
  gg
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


# When, for some reason, Ctrl+L stops working
clr <- function() system("clear")

# Erase all excluding a set of objects
cleanEnv <- cleanenv <- function(keep = NULL){
  if(is.null(keep)) keep <- 'no variables'
  myvars <- ls(pos = 1L)
  cat('Keeping:',paste0(keep, collapse = ', '),'\n')
  rm(list = myvars[!myvars %in% keep], pos = 1L)
  eval.parent(gc(), n = 1); gc()
  source('https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/code.R')
}

# limit length of a vector
vlimit <- function(x, mylim = Inf){
  if(length(x) > mylim) return(x[1:mylim]) else return(x)
}


# repeat char/number in a data.frame with x rows
dfrep <- function(x, n) data.frame(Identity=rep(as.character(n), length(x)), row.names = x, stringsAsFactors = F)

# Call when a process is taking too long and you wan to be notified when it finishes
process_has_ended <- function(x = NULL, email = "${USER}"){
  x <- ifelse(is.null(x), "Process", x)
  command <- paste0('echo "', x,  ' finished" | mail -s "Job status" "', email, '@lji.org"')
  cat(command, '\n')
  tvar <- system(command, intern = TRUE)
}
# process_has_ended()
