#!/usr/bin/R

#' Qlucore file generator
#'
#' This function creates a qlucore formated table.
#'
#' @param mat Expression matrix. It can be a Seurat object.
#' @param metadata Annotation data.frame.
#' @param rnames Features (genes).
#' @param cnames Sample/cells.
#' @param fname Path to the file you want to generate.
#' @param v Verbose
#' @keywords Qlucore
#'
#' @return Returns Qlucore format table [or writes it if fname is provided]
#'
#' @importFrom Seurat
#'
#' @export
#'
#' @examples
# ' qlucoref <- qlucore_format(mat = edata, metadata = mycells@meta.data)
#'

qlucore_format <- function(
  mat,
  metadata = NULL,
  rnames = NULL,
  cnames = NULL,
  fname = NULL,
  v = FALSE
){
  if(casefold(class(mat)) == 'seurat'){
    if(as.numeric(sub("(^...).*", "\\1", mat@version)) < 3) mat <- Seurat::UpdateSeuratObject(mat)
    if(v) cat('Seurat object found\n')
    metadata <- mat@meta.data
    cat("Assay", DefaultAssay(mat), "\n")
    mat <- GetAssayData(object = mat, slot = "data", assay = DefaultAssay(mat))
  }
  if(is.null(rnames)){
    rnames <- rownames(mat)
  }else{
    rnames <- rnames[rnames %in% rownames(mat)]
  }
  if(is.null(cnames)){
    cnames <- rownames(metadata)
  }else{
    cnames <- cnames[cnames %in% rownames(metadata)]
  }
  if(v){
    cat('Columns:', length(cnames), all(cnames %in% colnames(mat)), "\n")
    cat("Rows:", length(rnames), all(rnames %in% rownames(mat)), "\n")
  }
  mat <- mat[rnames, cnames]
  metadata <- metadata[cnames, ]
  # annotation with gap in Gene column
  qoldata <- cbind(Gene = rep("", ncol(metadata)), void = colnames(metadata), t(metadata))
  if(all(apply(metadata, 2, function(x) length(table(x)) != length(x) ))){ # check sample names
    qoldata <- rbind(Sample = c("", "Sample", rownames(metadata)), qoldata)
  }
  qoldata <- data.frame(qoldata, stringsAsFactors = F, check.names = F)
  if(v) cat('Building matrix\n')
  mat <- as.matrix(mat); mat <- round(mat, 2)
  if(v) cat(' Adding gap\n')
  qmat <- cbind(Gene = rownames(mat), void = rep("", nrow(mat)), mat) # expression with gap
  if(v) cat(' Proper format: data frame\n')
  qmat <- data.frame(qmat, stringsAsFactors = F, check.names = F) # between Gene names and matrix
  if(v) cat('Binding\n')
  qformat <- rbind(qoldata, rep("", ncol(qmat)), c("Gene_ID", rep("", ncol(qmat)-1)), qmat)
  colnames(qformat) <- rownames(qformat) <- NULL
  tvar <- grep("Gene_ID", qformat[, 1])
  if(v) print(qformat[c(1:5, tvar:(tvar + 2)), 1:4])
  if(!is.null(fname)){
    tvar <- list(c(".csv", ","), c(".txt", "\t"))[[1]]
    fname <- paste0(fname, "_qlucore_format_R", length(rnames), "xC", length(cnames), tvar[1])
    fname <- sub("\\/{2,}", "/", fname)
    write.table(qformat, file = fname, sep = tvar[2], col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  return(qformat)
}
