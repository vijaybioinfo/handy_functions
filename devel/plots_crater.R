#!/usr/bin/R

#' Crater plot
#'
#' This function creates a crater plot out of two Differential Gene Expression
#' Analysis.
#'
#' @param tests_list Named list of data.frames or named vector of files' paths.
#' Indicate the row names column after a '~'; default is 1. Eg.: /path/tab.csv~3
#' @param comp_names Vector of labels for comparison 1 and 2.
#' @param edataf Path to expression matrix (for stats calculation).
#' @param annotf Path to annotation data, sample names in the first column or
#' indicate the column after a '~'; eg.: /path/annot.csv~4
#' @param sample_filter Subsetting samples; eg.: c('column', 'group1'), or
#' list(c()).
#' @param feature_subset Vector of features to take into acount.
#' @param lfc Log Fold Change column.
#' @param pv [adjusted] P-values column.
#' @param lfcthresh Fold change threshold. To keep features over this threshold
#' coloured and sized, you need to pass the number as character.
#' @param keep_nas Keep features with NA in 'pv' or 'lfc' or give the name(s,
#' separated by a pipe '|', example "IL31|IL1R2").
#' @param topgenes Number, or a vector of genes [you can add a string
#' 'top(number)' if you need more genes added to the vector].
#' @param gene_pattern Pattern for genes you want to show on the plot.
#' @param column4stats Columns you want stats to be grouped by.
#' @keywords DGEA
#'
#' @return Returns table of crater data and stats if columns for it are given
#'
#' @importFrom ggplot2 super_stats_table
#'
#' @export
#'
#' @examples
# ' addinfo <- crater_plot(tests_list = c(comp1 = '/path1', comp2 = '/path2'))
#'

source("https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/filters.R")
# filters_subset_df
source("https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/utilities.R")
# show_commas, features_parse_ensembl, bordering, summ_tables
source("https://raw.githubusercontent.com/vijaybioinfo/handy_functions/master/devel/plots.R")
# make_breaks, make_title, plot_blank

crater_plot <- function(
  tests_list,
  edataf,
  comp_names = c("Comparison_1", "Comparison_2"),
  annotf = NULL,
  sample_filter = NULL,
  feature_subset = NULL,
  lfc = "log2FoldChange",
  pv = "padj",
  lfcthresh = 1,
  keep_nas = FALSE,
  topgenes = 5,
  gene_pattern = NULL,
  column4stats = NULL,
  gene_filter = list(mean = c("<=1", NA)),
  outputname = "./",
  return_out = FALSE,
  plot_interactive = FALSE,
  couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000'),
  verbose = FALSE
){
  suppressPackageStartupMessages(library(ggrepel))
  suppressPackageStartupMessages(library(dplyr))

  # Digest parameters
  if(length(tests_list) != 2) stop("You need two comparison to make a crater plot")
  if(is.character(tests_list[[1]])){
    if(length(unique(names(tests_list))) == 1){
      if(verbose) cat("Potentially identical comparison in both axes\n")
      names(tests_list) <- if(tests_list[[1]] != tests_list[[2]]){
        if(verbose) cat("Fixing names from file paths\n")
        if(length(unique(dirname(tests_list))) == 1){
          gsub("\\.[A-z]{,5}$", "", basename(tests_list))
        }else{ str_diff(strings = dirname(tests_list)) }
      }else{
        if(verbose) cat("Fixing names\n")
        paste0(names(tests_list), "_", 1:length(names(tests_list)))
      }
    }
  }
  if(is.null(names(tests_list))) names(tests_list) <- gsub(" {1,}", "_", comp_names)
  comp_names <- names(tests_list) <- gsub("\\-", "_", names(tests_list))

  # Creating file's name
  selection <- if(!is.null(sample_filter)){
    selection <- sapply(sample_filter, function(x) paste0(paste0(x[-1], collapse = ", "), " from ", x[1]) )
    paste0(selection, collapse = "; ")
  }else{ NULL }
  tmp <- paste0(casefold(comp_names), collapse = "_against_")
  if(!grepl("_$", outputname) && !dir.exists(outputname)) outputname <- paste0(outputname, "_")
  if(!grepl("/$", outputname) && dir.exists(outputname)) outputname <- paste0(outputname, "/")
  outputname <- paste0(outputname, 'FC', lfcthresh, "_", tmp)
  if(verbose) cat("Output at:", outputname, "\n")

  # Reading
  if(verbose) cat("Setting\n - Annotation\n"); rnames <- 1
  annot <- if(is.character(annotf)){
    if(grepl("^~", annotf)) annotf <- sub("~", paste0("/home/", system("echo $USER", intern = TRUE)), annotf)
    if(grepl("~", annotf)){ annotf <- unlist(strsplit(annotf, "~")); rnames <- as.numeric(annotf[-1]); annotf <- annotf[1] }
    readfile(annotf, row.names = ifelse(length(rnames), rnames, 1), check.names = FALSE)
  }else{ annotf }
  if(verbose) cat(" - Matrix of expression\n")
  edata <- if(is.character(edataf)){
    if(grepl("^~", edataf)) edataf <- sub("~", paste0("/home/", system("echo $USER", intern = TRUE)), edataf)
    if(grepl("~", edataf)){ edataf <- unlist(strsplit(edataf, "~")); rnames <- as.numeric(edataf[-1]); edataf <- edataf[1] }
    readfile(edataf, row.names = ifelse(length(rnames), rnames, 1), check.names = FALSE)
  }else{ edataf }
  rm(annotf, edataf); if(verbose) cat(" - Result tables\n")
  fstat <- if(is.character(tests_list[[1]])) lapply(tests_list, readfile, row.names = 1, stringsAsFactor = FALSE) else tests_list

  allgenes <- intersect(rownames(fstat[[1]]), rownames(fstat[[2]]))
  tvar <- intersect(feature_subset, allgenes)
  if(length(tvar)){ if(verbose) cat("Taking", length(tvar), "of", length(allgenes), "features\n"); allgenes <- tvar }
  # Get the Fold Changes and P-values
  if(verbose) cat("Getting LFCs and P-values from intersecting genes\n")
  fstat <- lapply(fstat, function(x){
    nas <- list(pv = is.na(x[, pv]), lfc = is.na(x[, lfc]))
    if(isTRUE(keep_nas)){ # Setting NAs to non-significant and 0 FC
      x[nas$pv, pv] <- 1; x[nas$lfc, lfc] <- 0
    }else if(is.character(keep_nas)){ # Keeping SOME genes if NA
      x[grepl(keep_nas, rownames(x)) & nas$pv, pv] <- 1
      x[grep(keep_nas, rownames(x)) & nas$pv, lfc] <- 0
    };
    data.frame(x[allgenes, ], stringsAsFactor = FALSE)
  })
  myfcs <- lapply(fstat, "[", c(lfc, pv))
  myfcs <- cbind(myfcs[[1]], myfcs[[2]]);
  colnames(myfcs) <- paste0(rep(names(fstat), each = 2), ".", names(myfcs))
  mydatafc <- data.frame(myfcs, stringsAsFactors = FALSE, row.names = rownames(fstat[[1]]))
  names(mydatafc) <- gsub(paste0("\\.", lfc), "", names(myfcs))
  mydatafc <- mydatafc[, c(names(fstat), grep(pv, names(mydatafc), value = TRUE))]
  mydatafc$min_padj <- apply(mydatafc[, grep(paste0("\\.", pv, "$"), names(mydatafc))], 1, min, na.rm = TRUE)
  mydatafc$min_fc <- apply(mydatafc[, grep(paste0("\\.", lfc, "$"), names(mydatafc))], 1, min, na.rm = TRUE)
  if(verbose) print(head(mydatafc))
  mydatafc <- mydatafc[complete.cases(mydatafc), ]
  mydatafc <- mydatafc[intersect(rownames(mydatafc), rownames(edata)), ]
  head(mydatafc)

  # Filtering
  thesecells <- if(!is.null(annot)) rownames(annot) else colnames(edata)
  if(!is.null(sample_filter) && !is.null(annot)){
    if(verbose) cat("Filtering samples\n")
    thesecells <- filters_subset_df(sample_filter, annot, v = TRUE)
    annot <- annot[thesecells, ]
  }
  thesecells <- intersect(thesecells, colnames(edata))
  if(verbose) cat("Sample number:", length(thesecells), "\n")
  if(length(thesecells) == 0) stop("No samples selected")
  mydatafc$mean <- rowMeans(edata[rownames(mydatafc), thesecells])

  # All of the filters must be passed
  mydatafc$filters <- rowSums(sapply(names(gene_filter), function(i){
    myfilter <- gene_filter[[i]]
    express <- paste0("tvar <- mydatafc[, '", i, "']", myfilter[[1]])
    eval(expr = parse(text = express))
    if(verbose) cat("Pass", i, myfilter[[1]], "filter:", sum(tvar), "- set to", as.numeric(myfilter[[2]]), "\n")
    return(tvar)
  })) > 0;
  mydatafc$lfcthresh <- rowSums(abs(mydatafc[, 1:2]) <= as.numeric(gsub("[A-z]", "", lfcthresh))) == 2
  if(is.character(lfcthresh)){
    lfcthresh <- as.numeric(gsub("[A-z]", "", lfcthresh))
    if(verbose) cat("Colouring and sizing only above", lfcthresh, lfc)
    if(verbose) cat(" -", sum(mydatafc$lfcthresh), "\n")
    mydatafc$filters <- (mydatafc$filters | mydatafc$lfcthresh)
  }
  if(verbose) cat("Modifying", sum(mydatafc$filters), "/", nrow(mydatafc), "samples\n")
  mycaption <- NULL;
  convert <- c(">" = "<=", "<" = ">=", "<=" = ">", ">=" = "<", "==" = "!=", "!=" = "=")
  for(i in names(gene_filter)){
    if(verbose) cat(i)
    myfilter <- gene_filter[[i]]
    mydatafc[mydatafc$filters, i] <- as.numeric(myfilter[[2]])
    # outputname <- paste0(outputname, "_", i, gsub(">|<|=", "", myfilter[[1]]))
    tvar <- names(convert) == gsub("[0-9]{1,}|\\.", "", myfilter[[1]])
    mycaption <- paste0(mycaption, i, sub(names(convert)[tvar], convert[tvar], myfilter[[1]]), "; ")
    if(verbose) cat("; ")
  }; if(verbose) cat("\n") #outputname <- paste0(outputname, "_", nrow(mydatafc), "features")

  # Renaming some columns
  mydatafc$min_padj <- apply(mydatafc[, grep(paste0("\\.", pv, "$"), names(mydatafc))], 1, min, na.rm = TRUE)
  mydatafc$min_padj[is.na(mydatafc$mean)] <- 1 # is.infinite(mydatafc$min_padj)
  mydatafc$min_padj[mydatafc$min_padj == 0] <- min(mydatafc$min_padj[mydatafc$min_padj > 0])
  mydatafc$significance <- -log10(mydatafc$min_padj)
  mydatafc$gene_name <- features_parse_ensembl(rownames(mydatafc))
  mydatafc$log2_mean <- log2(mydatafc$mean + 1)
  mydatafc$log2_mean[which(is.na(mydatafc$log2_mean) || mydatafc$log2_mean == 0)[1]] <- 0

  # Plotting border genes and with pattern
  borgenes <- topgenes
  topgenes <- if(any(grepl("^top[0-9]+", as.character(topgenes)))){
    as.numeric(gsub("top", "", grep("^top[0-9]+", topgenes, value = TRUE)))
  }else if(is.character(topgenes)){ topgenes }else{ topgenes[1] }
  borgenes <- list(
    if(is.numeric(topgenes)) bordering(mydatafc[!mydatafc$filters, ], cnames = 1:2, n = topgenes) else borgenes,
    if(!is.null(gene_pattern)) grep(gene_pattern, rownames(mydatafc), value = TRUE) else gene_pattern,
    borgenes
  )
  borgenes <- features_matching(unique(unlist(borgenes)), rownames(mydatafc))
  if(verbose) cat("Showing", length(borgenes), "genes\n")
  mydatafc_repel <- mydatafc[borgenes, ]
  ncolor <- make_breaks(c(0, mydatafc$log2_mean), n = length(couls)*2, push = c(max = 0.1))
  set.seed(27); aesy <- aes_string(x = comp_names[1], y = comp_names[2], color = 'log2_mean', size = 'significance')
  mysubtitle <- paste0(selection, ifelse(is.null(selection), "", "\n"))
  mysubtitle <- paste0(mysubtitle, 'Passed filters: ', sum(!mydatafc$filters))
  p <- ggplot(data = mydatafc) +
    geom_point(data = mydatafc[mydatafc$filters, ], mapping = aesy) +
    geom_point(data = mydatafc[!mydatafc$filters, ], mapping = aesy) +
    scale_color_gradientn(colours = couls, values = scales::rescale(ncolor), na.value = "#c4c4c4") +
    geom_hline(yintercept = c(-lfcthresh, lfcthresh), linetype = "dashed", alpha = 0.4) +
    geom_vline(xintercept = c(-lfcthresh, lfcthresh), linetype = "dashed", alpha = 0.4) +
    labs(
      x = make_title(comp_names[1]), y = make_title(comp_names[2]), caption = mycaption,
      subtitle = mysubtitle, color = 'Expression', size = 'Significance'
    ) + theme_bw() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(size = 2)
    )
  p <- p + guides(size = guide_legend(keywidth = 0.5, keyheight = 0.5, default.unit = "inch"))

  aesy <- aes_string(x = comp_names[1], y = comp_names[2], label = 'gene_name', size = 'significance')
  mydatafc_repel$significance[mydatafc_repel$significance == 0] <- 0.1; set.seed(27)
  p <- p + geom_text_repel(data = mydatafc_repel, mapping = aesy, color = 'black')
  return_list <- list(comp_stat = mydatafc, plot = p)

  if(!isTRUE(return_out)){
    if(verbose) cat("Plotting\n")
    # pdf(paste0(outputname, '.pdf'), 10, 10)
    # print(p)
    # dev.off()
    # pdf(paste0(outputname, '_blank.pdf'), 10, 10)
    # print(plot_blank(p))
    # dev.off()
    if(isTRUE(plot_interactive)){
      if(verbose) cat("Generating interactive files\n")
      suppressPackageStartupMessages(library(dplyr))
      suppressPackageStartupMessages(library(plotly))
      pathbk <- getwd(); setwd(dirname(outputname))
      ddfpassed <- mydatafc[!mydatafc$filters, ]
      ddfpassed$x <- ddfpassed[, comp_names[1]]; ddfpassed$y <- ddfpassed[, comp_names[2]]
      gp <- plot_ly(
        data = ddfpassed, x = ~x, y = ~y, color = ~log2_mean,
        text = ~paste("Gene: ", gene_name,
          '<br>Significance:', significance, '<br>Mean:', mean)
      ) %>% layout(
        xaxis = list(title = comp_names[1], showgrid = TRUE),
        yaxis = list(title = comp_names[2], showgrid = TRUE)
      )
      kk <- try(htmlwidgets::saveWidget(
        as_widget(gp), paste0(basename(outputname), '.html')
      ), silent = TRUE)
      # print(xtable::xtable(ddfpassed), type = "html", file = paste0(basename(outputname), '_table.html'))
      setwd(pathbk)
    }
  }

  # Creating table - maybe for supplementary
  if(!is.null(column4stats) && !is.null(annot)){
    if(verbose) cat("Creating stats table for", show_commas(column4stats, 4), "\n")
    fstat <- lapply(names(fstat), function(x){
      y <- fstat[[x]][, c(lfc, pv)]
      colnames(y) <- paste0(x, "_", colnames(y))
      y$gene_name <- paste0("'", rownames(y))
      y
    })
    mystats <- lapply(column4stats, function(x){
      stats_summary_table(
        mat = edata,
        groups = make_list(annot, x, grouping = TRUE),
        rnames = rownames(mydatafc),
        moments = c("mn", "p"),
        expr_cutoff = 10,
        v = TRUE
      )
    })
    mydatafc$gene_name = rownames(mydatafc)
    mytab <- summ_tables(ltabs = c(list(mydatafc), mystats), v = TRUE)

    smytab <- mytab
    smytab <- smytab[order(-smytab[, colnames(mydatafc)[1]]), ]
    dim(smytab); head(smytab)
    return_list$stats <- smytab
    if(!isTRUE(return_out)) write.csv(smytab, file = paste0(outputname, "_stats.csv"))
  }
  if(isTRUE(return_out)) return(return_list)
}

str_diff <- function(strings, str_split = "", sufix = TRUE){
  splitnames <- strsplit(strings, split = str_split)
  sapply(list(1:2, 2:1), function(x){
    y <- splitnames[[x[1]]] == splitnames[[x[2]]]
    y <- if(isTRUE(sufix)){
      max(which(diff(which(y)) == 1)):length(splitnames[[x[1]]])
    }else{ !y }
    z <- paste0(splitnames[[x[1]]][y], collapse = "")
    z <- gsub("\\/", "_", z)
    z
  })
}
