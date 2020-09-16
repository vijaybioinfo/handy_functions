#!/usr/bin/R

### Gene Set Enrichment Analysis ### -------------------------------------------
gsea_tests <- gsea_liger <- function(
  res,
  gene_name = "gene_name",
  lfc.type = 'log2FoldChange',
  gsea_file = "~/asthma/info/gsea_lists.csv",
  method = "liger",
  myseed = 27,
  path = "./",
  plot_all = TRUE,
  v = FALSE
){
  depend <- c('liger', 'fgsea', 'parallel', 'mdgsa')
  load_packs(depend, install = FALSE, v = FALSE)
  ncores <- min(c(ceiling(detectCores() / 5), 3))
  prefix <- ifelse(grepl("\\/$", path), path, paste0(path, "_"))
  if(v) cat("Output at:", prefix)
  if(!dir.exists(prefix) && grepl("\\/$", prefix)){
    if(v) cat(" - created"); dir.create(prefix)
  }else if(!dir.exists(dirname(prefix))){
    if(v) cat(" - created"); dir.create(dirname(prefix))
  }; if(v) cat("\n")
  res <- res[!is.na(res[, lfc.type]), ]
  vals <- res[, lfc.type]
  guniverse <- unique(sub("'", "", res[, gene_name]))  # get universe
  names(vals) <- guniverse

  # get a gene set
  if(file.exists(gsea_file[[1]][1])){
    if(v) cat("Getting gene sets file\n")
    gset_list <- readfile(gsea_file, stringsAsFactors = F, header = TRUE, v = v)
    # l2add <- list(
    #   Singhania_il13_induced = c("POSTN", "SERPINB2", "CLCA1"), # sputum
    #   Singhania_il17_induced = c("CXCL1", "CXCL2", "CXCL3", "IL8", "CSF3"),
    #   Micosse_th9 = c("PPARG", "IL17R", "FAM169A", "LYL1", "TCEAL4", "RP11-1018N14.5", "SH3TC1", # this is the whole list
    #     "PXDC1", "PLA2G16", "IRX5", "RGL1", "IL6R", "EFCAB8", "FANK1", "NCKAP5L", "PRR5L",
    #     "YBX3", "FAAHP1", "NRP3", "TLR6", "ZEB2", "CD84", "CCR4", "LEF1", "STARD8", "HLF",
    #     "ANPEP", "CCL17", "GPR55", "PPP1R26", "PALLD", "ADGRB2", "SNED1", "STRP2", "PTGDR2",
    #     "HRASLS5", "SCAMP5", "RIMKLB", "GRASP", "CCL22"),
    #   Micosse_th9_th2 = c("PPARG", "IL17R", "IL6R", "CCR4", "CCL17", "PTGDR2", "CCL22") # TH2 genes
    # )
    # for(iname in names(l2add)){
    #   gset_list[[iname]] <- ""; gset_list[[iname]][1:length(l2add[[iname]])] <- l2add[[iname]]
    # }; #write.csv(gset_list, file = "gsea_lists.csv", row.names = FALSE)
  }else if(length(gsea_file)){
    if(v) cat("Given", class(gsea_file), "of genes\n")
    gset_list <- if(!is.list(gsea_file)) list(gene_list = gsea_file) else gsea_file
  }else{
    if(v) cat("Human Gene Ontology to HUGO Symbol list\n",
      "http://geneontology.org/page/download-go-annotations\n")
    gset_list <- org.Hs.GO2Symbol.list[1:5]
  }; gset_list <- lapply(gset_list, function(x){ # removing NAs and empty elements
    y <- x[!is.na(x)]; gsub("'| ", "", y[y != ""])
  }); if(is.null(names(gset_list))) names(gset_list) <- paste0("g", 1:length(gset_list))
  if(v) cat("Number of sets:", length(gset_list), "\n")
  if(v) str(gset_list[head(1:length(gset_list))])

  ## Correcting for multiple-hypothesis
  if(method == "liger"){
    fgseaRes <- bulk.gsea(vals, gset_list)
    colnames(fgseaRes) <- c("pval", "padj", "ES", "edge")
    fgseaRes <- cbind(pathway = rownames(fgseaRes), fgseaRes)
    fgseaRes$size = sapply(gset_list, function(x) sum(x %in% guniverse) )
    fgseaRes <- data.table(fgseaRes)
    myheight <- 9
  }
  if(method == "fgsea"){
    fgseaRes <- fgsea(
      pathways = gset_list,
      stats = vals,
      minSize = 3,
      maxSize = 500,
      nproc = ncores,
      nperm = 10000
    )
    myheight <- 7
    str(fgseaRes)
  }
  cat("Plotting top enrichements\n")
  topPathways <- c(
    fgseaRes[ES > 0, ][head(order(pval), n = 10), pathway],
    rev(fgseaRes[ES < 0, ][head(order(pval), n = 10), pathway])
  )
  if(length(topPathways) > 2){
    fname <- paste0(prefix, "a1_", method, '_top_up_and_down_pathways_', lfc.type, '.pdf')
    pdf(sub("\\/a1$", "_a1", fname), width = 8, ifelse(length(topPathways) >= 16, 16, 8))
    try(plotGseaTable(gset_list[topPathways], vals, fgseaRes, gseaParam = 1))
    graphics.off()
    tvar <- if("leadingEdge" %in% colnames(fgseaRes)) ncol(fgseaRes) - 1 else ncol(fgseaRes)
    write.csv(data.frame(fgseaRes[, 1:(tvar)]), file = sub("pdf", "csv", fname), row.names = FALSE)
  }
  if(!plot_all) return(fgseaRes)

  go_names <- suppressWarnings(getGOnames(names(gset_list), verbose = TRUE)); names(go_names) <- names(gset_list)
  go_names <- gsub(" activity|activity", "", sub("(^[A-z]{41,}).*", "\\1", go_names))
  if(v) cat("\n\n")
  for(indx in names(gset_list)){
    cat("----------------\n", indx, "\n")
    allgs <- sub("'", "", gset_list[[indx]]); allgs <- unique(allgs[allgs != ""])
    gs <- getfound(allgs, sour = guniverse, element = 'genes', v = v)
    if(!length(gs) > 1) next
    suffix <- paste0("_", method, "_", length(gs), "of", length(allgs))

    if(!is.na(go_names[indx])){
      if(v) cat(go_names[indx], "\n")
      indxname <- sub("^GO\\:", "", indx)
    }else{ indxname <- indx }

    fname <- paste0(prefix, indxname, "_", lfc.type, suffix, '.pdf')
    pdf(fname, height = myheight, width = 7)
    if(v) cat('Testing', lfc.type, '\n'); set.seed(myseed)
    pval <- paste0("adjusted P-value: ", round(fgseaRes[pathway == indx, "padj"], 5))
    if(method == "liger"){
      print(gsea(values = vals, geneset = gs, mc.cores = ncores, return.details = TRUE))
      if(!is.na(go_names[indx])) graphics::legend("bottomleft", legend = go_names[indx], bty = 'n', horiz = TRUE)
      graphics::legend("bottomright", legend = pval, bty = 'n', horiz = TRUE)
    }
    if(method == "fgsea"){
      p <- plotEnrichment(gs, vals) + labs(caption = pval)
      print(p)
    }
    graphics.off()
    # cat('- random gene set\n');set.seed(myseed); randomgs <- sample(names(vals), length(gs))
    # print(gsea(values = vals, geneset = randomgs, mc.cores = ncores, return.details = TRUE))
  }
  return(fgseaRes)
}

# Metrics for Ranking Genes
# Extracted from https://software.broadinstitute.org/gsea/doc/GSEAUserGuideFrame.html
# groups = make_list(annot, "celltypes", grouping = T)
gsea_metric <- function(
  groups,
  mat,
  metric = c('Signal2Noise', 'tTest', 'Ratio_of_Classes', 'Diff_of_Classes', 'log2_Ratio_of_Classes'),
  rnames = NULL,
  v = FALSE
) {
  metric <- match.arg(metric)
  if(is.null(rnames)) rnames <- rownames(mat)
  tvar <- unique(groups)
  if(v) cat(commas(tvar), "\n")
  groups <- sub(paste0("^", tvar[1], "$"), "X1", groups)
  groups <- sub(paste0("^", tvar[2], "$"), "X2", groups)
  tvar <- sub(paste0("^", tvar[1], "$"), "X1", sub(paste0("^", tvar[2], "$"), "X2", tvar))
  if(v) cat(commas(tvar), "\n")
  stattab <- get_stat_report(
    mat = mat,
    groups = groups,
    rnames = rnames,
    moments = c('mn', 'sd'),
    v = v
  ); str(stattab)
  stattab2 <- data.frame(sapply(colnames(stattab), function(x){
    y <- stattab[, x]
    if(grepl("_sd", x)){
      ttab <- as.matrix(data.frame(minsd = .2 * abs(stattab[, sub("_sd", "_mean", x)]), sd = y))
      y <- matrixStats::rowMaxs( ttab )
    }
    if(grepl("_mean", x)) y <- ifelse(y == 0, 1, y)
    y
  }), row.names = rownames(stattab)); stattab <- stattab2; rm(stattab2)
  dmean <- apply( stattab[, rev(grep("_mean", colnames(stattab)))], 1, diff )
  if(v) cat(metric, "\n")
  ymetric <- switch(metric,
    Signal2Noise = {
      if(v) cat("X1mean - X2mean\n")
      if(v) cat("---------------\n")
      if(v) cat("  X1sd + X2sd\n")
      dmean / apply( stattab[, grep("_sd", colnames(stattab))], 1, sum )
    }, tTest = {
      if(v) cat("          X1mean - X2mean\n")
      if(v) cat("  ------------------------------\n")
      if(v) cat("sqrt ((X1sd^2 / X1n) +( X2sd^2 / X2n))\n")
      ttab <- sweep(stattab[, grep("_sd", colnames(stattab))] ** 2, 2, table(groups), "/")
      dmean / sqrt( base::rowSums(ttab) )
    }, Ratio_of_Classes = {
      if(v) cat("X1mean / X2mean\n")
      stattab[, "X1_mean"] / stattab[, "X2_mean"]
    }, Diff_of_Classes = {
      if(v) cat("X1mean - X2mean\n"); dmean
    }, log2_Ratio_of_Classes = {
      if(v) cat("log2(X1mean / X2mean)\n")
      log2(stattab[, "X1_mean"] / stattab[, "X2_mean"])
    }
  )
  ymetric
}
