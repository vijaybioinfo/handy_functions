#!/usr/bin/R

##################
# Gene set lists #
##################

# This scripts gathers all gene sets used at Vijay Lab since 2017.
# They are all human [gene] symbols.

source("/home/ciro/scripts/handy_functions/devel/file_reading.R") # readfile
source("/home/ciro/scripts/handy_functions/devel/overlap.R") # overlap_list
source("/home/ciro/scripts/handy_functions/devel/plots.R") # make_title
source("/home/ciro/scripts/handy_functions/devel/filters.R") # getDEGenes
source("/home/ciro/scripts/handy_functions/devel/utilities.R") # vlist2df

### Initial signatures ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnames = c(
  "/home/ciro/simon/info/gsea_lists_extended.csv",
  "/home/ciro/vdv/vfajardo/general_data/cell_type_classification/module_definition/module_signatures/activation/activation_2020-02-24_NoFiltering/ModuleACTIVATIONFeaturesSignatureWithLFCTholds1.txt",
  "/home/ciro/large/covid19/results/dgea/CD8T24_R1n2N_sng_15p/comprs/activation/NCVvsCV/signature_genes/selected_genes.csv",
  "/home/ciro/covid19/info/signatures_hdm.csv"
)
# "~/asthma_pjs/info/trm_suptab2Lung_clarke.csv", # all genes taken; /home/ciro/asthma_airways/scripts/trm_madness.R
signatures_vijaylab <- lapply(fnames, readfile, stringsAsFactors = FALSE)
colnames(signatures_vijaylab[[2]]) = "tcell_activation_flu"
signatures_vijaylab <- unlist(signatures_vijaylab, recursive = FALSE)
signatures_vijaylab <- lapply(signatures_vijaylab, function(x){ y <- x[x != ""]; y[!is.na(y)] })
str(signatures_vijaylab)

# Adding James signatures
# cp ~/Documents/liai/cancer/amica/james_trm_table4.csv /Volumes/ciro/amica/info/
# # unique genes in TRM clusters. 2 is purple cluster and hast ITGAE
# cp ~/Documents/liai/cancer/amica/james_trm_table7.csv /Volumes/ciro/amica/info/
fnames = c(
  lung = "/home/ciro/asthma_pjs/info/deprecated/james_trm_table2.csv",
  tumor = "/home/ciro/amica/info/james_trm_table4.csv",
  sctrm = "/home/ciro/amica/info/james_trm_table7.csv"
)
res_list <- lapply(fnames, readfile, stringsAsFactors = FALSE, row.names = 1)
str(res_list, list.len = 10)
colnames(res_list[[1]])[1] <- "log2FoldChange"
colnames(res_list[[2]])[c(1,3)] <- c("log2FoldChange", "padj")
range(abs(res_list[[1]]$log2FoldChange))
range(abs(res_list[[2]]$log2FoldChange))
degs_list0 <- lapply(
  X = res_list[1:2],
  FUN = function(x){
    list(
      trm = getDEGenes(x, pv = 0.05, fc = 1, upreg = TRUE, v = TRUE),
      nontrm = getDEGenes(x, pv = 0.05, fc = 1, upreg = FALSE, v = TRUE)
    )
})
cat("All genes:", sum(sapply(degs_list0[[1]], length)) == nrow(res_list[[1]]), "\n")
cat("All genes:", sum(sapply(degs_list0[[2]], length)) == nrow(res_list[[2]]), "\n")
str(degs_list0)
degs_list1 = lapply(res_list[3], make_list, "Cluster.where.gene.is.enriched")
names(degs_list1[[1]]) <- paste0("cluster", names(degs_list1[[1]]))
str(degs_list1)
degs_list <- c(degs_list0, degs_list1)
degs_list <- unlist(degs_list, recursive = FALSE)
names(degs_list) <- gsub("(.*)\\.(.*)", "\\2_\\1_clarke", names(degs_list))
str(degs_list)
signatures_vijaylab = c(signatures_vijaylab, degs_list)

### COVID-19 CD8 paper ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signaturesl <- read.csv("/home/ciro/covid19/info/tcell_cd8story_signatures.csv", stringsAsFactor = FALSE)
signaturesl <- as.list(signaturesl)
signaturesl <- signaturesl[!grepl("tcell_cytotoxic_guo|arnon|dixhaust", names(signaturesl))]
## Converting to human
signaturesl <- lapply(signaturesl, function(x){ # removing NAs and empty elements
  y <- x[!is.na(x)]; gsub("'| ", "", y[y != ""])
});
str(signaturesl)
tem_signatures <- read.csv("/home/ciro/covid19/info/tcell_cd8story_signatures_2020_06_21.csv", stringsAsFactor = FALSE)
tem_signatures <- lapply(as.list(tem_signatures), function(x){ # removing NAs and empty elements
  y <- x[!is.na(x)]; gsub("'| ", "", y[y != ""])
});
signatures_cd8covid <- c(tem_signatures, signaturesl)
str(signatures_cd8covid)

### More signatures from Preethi ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fnames = list.files(path = '/home/ciro/preethi/info/tem_signatures_scrnaseq', pattern = 'txt', full.names = TRUE)
mynames <- make_title(gsub("&|,|.txt|genes|from|pathway", "", basename(fnames), ignore.case = TRUE))
mynames <- casefold(gsub(" ", "_", mynames))
signatures_tem <- suppressWarnings(lapply(fnames, readfile, stringsAsFactors = FALSE))
signatures_tem <- unlist(signatures_tem, recursive = FALSE)
names(signatures_tem) <- mynames
str(signatures_tem)

# ### From the internet ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# # extra_signatures = sigPathway::importGeneSets(fileNames = "GMT, GMX, GRP, and XML")
# # https://cran.r-project.org/web/packages/msigdbr/vignettes/msigdbr-intro.html
# urls = "http://www.gsea-msigdb.org/gsea/msigdb/download_geneset.jsp?geneSetName=KRAS.600.LUNG.BREAST_UP.V1_UP&fileType=gmx"
# msigdbr::msigdbr_species()
# h_gene_sets = msigdbr::msigdbr(species = "Homo sapiens")
# h_gene_sets
# selected_set = grep(pattern = "LUNG_UP", h_gene_sets[['gs_name']], ignore.case = TRUE)
# selected_set = h_gene_sets[['gs_name']] %in% "KRAS.LUNG_UP.V1_UP"
# table(h_gene_sets[['gs_name']][selected_set])
# extra_signatures = make_list(data.frame(h_gene_sets[selected_set, ]), "gs_name", "gene_symbol")

# ### check with previous object ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# signatures_vijaylab0 = readfile("/home/ciro/scripts/handy_functions/data/signatures_vijaylab.rdata")
# signatures_cd8covid0 = readfile("/home/ciro/scripts/handy_functions/data/signatures_cd8covid.rdata")
# signatures_tem0 = readfile("/home/ciro/scripts/handy_functions/data/signatures_tem.rdata")
# names(signatures_vijaylab0)[!names(signatures_vijaylab0) %in% names(signatures_vijaylab)]
# names(signatures_cd8covid0)[!names(signatures_cd8covid0) %in% names(signatures_cd8covid)]
# names(signatures_tem0)[!names(signatures_tem0) %in% names(signatures_tem)]
# # all cool!

### Merging ### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signatures_vijay_db = c(signatures_vijaylab, signatures_cd8covid, signatures_tem)
length(names(signatures_vijay_db)) == length(unique(names(signatures_vijay_db))) # no repeated names
signatures_vijay_db_ov = sapply(
  X = signatures_vijay_db,
  FUN = function(x){
    sapply(signatures_vijay_db, function(y) mean(x %in% y) )
})
diag(signatures_vijay_db_ov) = 0
compare_lists = function(x, y) list(missing_x = x[! x %in% y], missing_y = y[! y %in% x])
compare_lists(
  signatures_vijay_db[['hallmarkglycolisis_signature_gsea']],
  signatures_vijay_db[['hallmark_glycolysis']]
) # looks like it was updated?
allin = which(signatures_vijay_db_ov == 1, arr.ind = TRUE); allin # one :o
signatures_vijay_db_ov[allin[, 1], allin[, 2], drop = FALSE]
str(signatures_vijay_db[unlist(dimnames(signatures_vijay_db_ov[allin[, 1], allin[, 2], drop = FALSE]))])

ddf <- vlist2df(signatures_vijay_db); ddf[is.na(ddf)] <- ""; ddf = rbind(ddf[100000, ], ddf)
str(ddf)
# https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMX:_Gene_MatriX_file_format_.28.2A.gmx.29
write.table(
  x = ddf, file = "/home/ciro/scripts/handy_functions/data/signatures_vijay_db.gmx",
  row.names = FALSE, quote = FALSE, sep = "\t"
)
write.table(
  x = ddf, file = "/home/ciro/scripts/handy_functions/data/signatures_vijay_db.txt",
  row.names = FALSE, quote = FALSE, sep = "\t"
)
save(signatures_vijay_db, file = "/home/ciro/scripts/handy_functions/data/signatures_vijay_db.rdata")
