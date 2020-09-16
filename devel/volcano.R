#!/usr/bin/R

# volcano adding a new continuous colour to the points

# source('/mnt/BioHome/ciro/scripts/functions/volcano_variant_color.R')
volplot <- function(
  x,
  pvalth = 0.05,
  lfcth = 0.5,
  pvaltype = 'padj',
  lfctype = 'log2FoldChange',
  col_feature = NULL,
  size_feature = NULL,
  gene_name = NULL,
  group = NULL,
  do_fdr_log = TRUE,
  interact = FALSE,
  ngenes = 5,
  clipp = FALSE,
  legends = TRUE,
  titl = "Volcano Plot",
  subt = NULL,
  check_genes = FALSE, # named list or a vector; if only text: list(text = c())
  gname = "Markers",
  grepel = TRUE,
  only_degs = FALSE,
  return_plot = FALSE,
  mycolours = NULL, # groups colours
  couls = c('#ffdf32', '#ff9a00', '#ff5a00', '#ff5719','red2','#b30000', '#670000'), # continuous colours
  v = FALSE
){
  x <- data.frame(x, stringsAsFactors = F, check.names = FALSE)
  # add a grouping column; default value is "not significant"
  if(is.null(gene_name)){
    x["gene_name"] <- rownames(x)
  }else{
    if("gene_name" %in% colnames(x) && gene_name != "gene_name"){
      colnames(x) <- sub("gene_name", "pgene_name", colnames(x))
    }
    colnames(x) <- sub(paste0("^", gene_name, "$"), "gene_name", colnames(x)); gene_name <- "gene_name"
  }
  x$saved_name <- gsub("'", "", rownames(x)); x$gene_name <- parse_ens_name(x$gene_name)
  kk <- getfound(c(pvaltype, lfctype, col_feature, size_feature, gene_name), colnames(x), v = v)
  rownames(x) <- x$gene_name
  colnames(x) <- sub(paste0("^", pvaltype, "$"), "FDR", colnames(x))
  colnames(x) <- sub(paste0("^", lfctype, "$"), "Fold", colnames(x))
  x <- x[, !duplicated(colnames(x))]
  if(is.null(group)){
    if(v) cat("Setting groups\n")
    x["group"] <- "Not_significant"
    # change the grouping for the entries with significance but not a large enough Fold change
    x[which(abs(x['FDR']) < pvalth & abs(x['Fold']) > lfcth ), "group"] <- "Significant"

    # change the grouping for the entries a large enough Fold change but not a low enough p value
    x[which(abs(x['FDR']) > pvalth & abs(x['Fold']) > lfcth ), "group"] <- "FC"

    # change the grouping for the entries with both significance and large enough fold change
    x[which(abs(x['FDR']) < pvalth & abs(x['Fold']) > lfcth ), "group"] <- "DEG"
  }else{
    if(v) cat("Groups column given\n")
    x["group"] <- x[group]
  }
  # wont work if the name has "log" as in logged already
  if(isTRUE(do_fdr_log) && !grepl("log", casefold(pvaltype))){
    x[, 'FDR'] <- (-log10(x[, 'FDR']))
    pvaltype <- paste0("FDR: -log10(", pvaltype, ")")
    lpvalth <- round(-log10(pvalth), 3)
  }else{ lpvalth <- pvalth }
  # x <- x[order(x['FDR']), ]
  # if(max(abs(diff(head(x[, 'FDR'], 40)))) > 100) # checking if there's a huge difference between the top
  x[is.infinite(x[, 'FDR']), 'FDR'] <- max(x[is.finite(x[, 'FDR']), 'FDR'], na.rm = TRUE)
  x$group1 <- x$group

  if(v) cat("Find given markers and label in groups\n")
  lcheck_genes <- if(is.list(check_genes)) check_genes else list(check_genes)#, c("IL9", "IL17F", "IFNG", "IL5"))
  if(is.null(names(lcheck_genes))) names(lcheck_genes) <- paste0(gname, 1:length(lcheck_genes))
  lcheck_genes <- overlap_calc(lcheck_genes)
  lcheck_genes <- lcheck_genes[sapply(lcheck_genes, length) > 0]
  all_genes <- character(0)
  gname <- names(lcheck_genes)
  if(names(lcheck_genes)[1] == "text" & length(lcheck_genes) == 1){
    if(v) cat("Check names only as text\n")
    lcheck_genes <- FALSE; ngenes = 0
    all_genes <- unname(unlist(check_genes))
  }else{
    for(i in 1:length(lcheck_genes)){
      check_genes <- lcheck_genes[[i]]
      if(check_genes[1] != "" && is.character(check_genes[1])){
        if(v) cat(names(lcheck_genes[i]))
        tvar <- check_genes %in% rownames(x)
        if(sum(!tvar)) warning("Some genes are not in data: ", commas(check_genes[!tvar]))
        tvar <- check_genes[tvar]
        if(length(tvar)) x[tvar, "group"] <- names(lcheck_genes[i]) else tvar <- ''
      }else{ tvar <- '' }; if(v) cat(sum(tvar != ""), "features\n")
      all_genes <- c(all_genes, tvar)
    }
    all_genes <- unique(all_genes)
  }

  if(ngenes > 0){
    if(v) cat("Top", ngenes, "genes\n")
    check_genes <- unique(c(bordering(x, cnames = "FDR", ngenes), bordering(x, cnames = "Fold", ngenes)))
    # check_genes <- head(rownames(x[with(x, order(Fold, FDR)),]), ngenes)
    # check_genes <- unique(c(check_genes, head(rownames(x[with(x, order(-Fold, FDR)),]), ngenes)))
    # check_genes <- unique(c(rbind(check_genes, head(rownames(x[with(x, order(-FDR)), ]), ngenes))))
  }else if(!is.character(check_genes[1])){
    check_genes <- rownames(x[with(x, order(-FDR)), ])[1]
  }
  check_genes <- check_genes[!is.na(check_genes)]
  if(all_genes[1] != "") check_genes <- unique(c(head(all_genes, Inf), check_genes))
  check_genes <- x[check_genes, ]
  # check_genes <- x[x$group == 'DEG', ]

  if(v) cat("Found", commas(unique(x[, "group"]), 10), "\n")
  mycats <- unique(unlist(x["group"]))
  if(v) cat("Setting colours and labels\n")
  tcols <- c("#b2c5ca", "#008080", "#b7dac5", "#b22028")
  names(tcols) <- c("Not_significant", "Significant", "FC", "DEG")
  mycolours <- c(tcols, v2cols(mycats[!mycats %in% names(tcols)], mycolours, fw = TRUE, v = v))
  signf <- paste0(lfctype, ": ", lfcth, " / ", pvaltype, ": ", pvalth)
  if(is.character(x$group1)) signf <- paste(signf, "\nNo. of DEGs:", sum(x$group1 %in% 'DEG'))
  if(!is.null(subt)){
    subt <- paste(subt, " | ", signf)
  }else{ subt <- signf }
  capty <- NULL; if(length(check_genes) > 1){
    tvar <- table(check_genes[!check_genes$group %in% c("Not_significant", "DEG"), "group"])
    capty <- paste0(paste0(names(tvar), ": ", unname(tvar)), collapse = "; ")
  }

  if(!interact){
    if(v) cat("-- Static plot\n")
    p <- ggplot(x, aes(x = Fold, y = FDR, label = gene_name))
    if(isTRUE(clipp)) clipp <- 0.999
    if(is.numeric(clipp) || is.character(clipp)){
      if(is.character(clipp)){
        clippit <- as.numeric(sub("c", "", clipp))
        clippit[2] <- max(abs(x[, 'Fold'])) + min(abs(diff(abs(x[, 'Fold']))))
      }else if(clipp < 1){
        clippit <- quantile(abs(x[, 'FDR']), prob = clipp)
        clippit[2] <- quantile(abs(x[, 'Fold']), prob = clipp)
      }else{
        clippit <- x[tail(order(abs(x[, 'FDR'])), clipp + 1)[2], 'FDR']
        clippit[2] <- x[tail(order(abs(x[, 'Fold'])), clipp + 1)[2], 'Fold']
      }
      p <- p + coord_cartesian(
        ylim = c(ifelse(any(x['FDR'] < 0, na.rm = TRUE), -clippit[1], 0), clippit[1]),
        xlim = c(ifelse(any(x['Fold'] < 0, na.rm = TRUE), -clippit[2], 0), clippit[2])
      )
    }
    if(length(couls) == 2){
      mypalette <- colorRampPalette(couls, space = 'Lab')
      couls <- mypalette(4)
    }
    colbar <- scale_color_manual(values = mycolours)
    aesy <- aes_string(colour = 'group')
    if(isTRUE(col_feature %in% colnames(x))){
      if(v) cat("Colouring based on", col_feature)
      brks <- make_breaks(x[, col_feature], n = length(couls), push = c(max = 0.1))
      colbar <- scale_color_gradientn(name = col_feature, colours = couls, na.value = 'gray',
        breaks = brks, labels = brks, limits = c(brks[1], max(brks)), guide = "colorbar")
      aesy <- aes_string(colour = col_feature)
    }
    if(isTRUE(col_feature %in% colnames(x)) && isTRUE(size_feature %in% colnames(x))){
      if(v) cat(" and size on", size_feature)
      scasize <- scale_size(range = c(0, 3))#, breaks = seq(0, 3, length.out = 4),
        # labels = make_breaks(x = x[, size_feature], n = 4))
      aesy <- aes_string(colour = col_feature, size = size_feature)
    }; if(v) cat("\n")
    # if(is.character(grepel)){
    #   tvar <- paste(names(aesy), "=", gsub("~", "", as.character(aesy)), collapse = ", ")
    #   tvar <- paste("aesy <- aes_string(", tvar, ",", grepel, "= 'gene_name' )"); print(tvar)
    #   eval(expr = parse(text = tvar))
    # }
    p <- p + geom_point(data = x, mapping = aesy) + colbar
    p <- p + geom_hline(yintercept = lpvalth, linetype = "dashed", color = "gray") +
      geom_vline(xintercept = lfcth, linetype = "dashed", color = "gray")
    if(any(x['Fold'] < 0, na.rm = TRUE)) p <- p + geom_vline(xintercept = -lfcth, linetype = "dashed", color = "gray")
    if(any(x['FDR'] < 0, na.rm = TRUE)) p <- p + geom_hline(yintercept = -lpvalth, linetype = "dashed", color = "gray")
    if(any(-x$Fold > 0, na.rm = TRUE)) p <- p + xlim(c(-max(abs(x$Fold), na.rm = TRUE), max(abs(x$Fold), na.rm = TRUE)))
    if(any(-x$FDR > 0, na.rm = TRUE)) p <- p + ylim(c(-max(abs(x$FDR), na.rm = TRUE), max(abs(x$FDR), na.rm = TRUE)))
    if(only_degs){
      if(v) cat("DEGs only colouring\n")
      # subdata <- dplyr::filter(x, abs(Fold) <= lfcth | FDR >= pvalth)
      p <- p + stat_density2d(data = x[!x$group %in% c('DEG', gname), ],
          aes(fill = ..density..^0.1), geom = "raster",  contour = FALSE, show.legend = FALSE) +
        scale_fill_gradientn(colours = colorRampPalette(c("white", blues9))(256), breaks = NULL) +
        geom_point(data = x[x$group %in% c('DEG', gname), ], aes(color = group))
    }else if(gname %in% check_genes$group){
      p <- p + geom_point(data = check_genes[check_genes$group %in% gname, ],
        aes(x = Fold, y = FDR, fill = group), shape = 21, size = 2) +
        scale_fill_manual(values = mycolours[names(mycolours) %in% gname])
    }

    if(length(rownames(check_genes)) && isTRUE(grepel)){
      if(v) cat("Repeling names\n")
      p <- p + ggrepel::geom_text_repel(
          # data = check_genes[bordering(dtab = check_genes, cnames = "Fold", n = 2*ngenes), ],
          data = check_genes,
          aes(label = gene_name),#  size = 7,
          box.padding = unit(0.35, "lines"),
          point.padding = unit(0.3, "lines")
        )
    }
    p <- p + theme_classic() + #+ theme_classic(base_size = 2)
      labs(title = titl, subtitle = subt, x = lfctype, y = pvaltype, caption = capty) +
      theme(legend.position = ifelse(isTRUE(legends), 'right', 'none'),
        plot.title = element_text(face = "bold", hjust = 0.5),
        panel.background = element_blank())
      if(exists('scasize')) p <- p + scasize #scale_size_binned()
      if(!legends) p <- p + theme(plot.margin = unit(c(0.1, 2, 0.1, 0.1), "cm"))
  }else{
    if(v) cat("-- Interactive plot\n"); a <- list()
    if(v) cat("Bulding layout\n")
    for (i in seq_len(nrow(check_genes))) {
      m <- check_genes[i, ]
      a[[i]] <- list(
        x = m[["Fold"]],
        y = m[["FDR"]],
        text = m[["gene_name"]],
        xref = "x",
        yref = "y",
        showarrow = TRUE,
        arrowhead = 0.5,
        ax = 20,
        ay = -40
      )
    }
    titl <- paste0(titl, '\n', subt)
    p <- plot_ly(data = x, x = ~Fold, y = ~FDR, text = ~gene_name, mode = "markers", color = ~group) %>%
      layout(title = titl)
    if(isTRUE(legends)) p <- p %>% layout(annotations = a)
  }
  if(return_plot || interact) return(p) else{ cat("Plotting\n"); print(p)}
}
