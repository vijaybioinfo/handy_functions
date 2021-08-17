#!/usr/bin/R

#' Qlucore file generator
#'
#' This function creates a qlucore formated table.
#'
#' @param edata Expression matrix. It can be a Seurat object.
#' @param metadata Annotation data.frame.
#' @param rnames Features (genes).
#' @param cnames Sample/cells.
#' @param fname Path to the file you want to generate.
#' @param verbose Show progress.
#' @keywords Qlucore
#'
#' @return Returns Qlucore format table [or writes it if fname is provided]
#'
#' @importFrom Seurat
#'
#' @export
#'
#' @examples
# ' qlucoref <- qlucore_format(edata = edata, metadata = mycells@meta.data)
#'

qlucore_format <- function(
  edata,
  metadata = NULL,
  rnames = NULL,
  cnames = NULL,
  fname = NULL,
  threshold = 0,
  verbose = FALSE
){
  if(casefold(class(edata)) == 'seurat'){
    if(verbose) cat('Seurat object found\n')
    if(grepl("^2.|^1.", edata@version))
      edata <- Seurat::UpdateSeuratObject(edata)
    metadata <- edata@meta.data
    if(verbose) cat("Assay", DefaultAssay(edata), "\n")
    edata <- edata@assays$RNA@data #Seurat::GetAssayData(object = edata)
  }
  if(is.null(rnames)) rnames <- rownames(edata)
  rnames <- show_found(x = rnames, rownames(edata), verbose = verbose)
  if(is.null(cnames)) cnames <- rownames(metadata)
  cnames <- show_found(x = cnames, rownames(metadata), verbose = verbose)
  edata <- edata[rnames, cnames]
  gmeans <- sort(Matrix::rowMeans(edata), decreasing = TRUE)
  rnames_sub <- names(gmeans)[gmeans > threshold]
  if(verbose){
    cat("Filter:", threshold, "\nRange of mean expression:\n")
    print(c(head(gmeans), tail(gmeans)))
    cat(length(rnames_sub), "/", nrow(edata), 'kept\n'); str(rnames_sub)
  }; edata <- edata[rnames_sub, ]
  metadata <- metadata[cnames, ]
  # annotation with gap in Gene column
  qoldata <- cbind(Gene = rep("", ncol(metadata)), void = colnames(metadata), t(metadata))
  if(all(apply(metadata, 2, function(x) length(table(x)) != length(x) ))){ # check sample names
    qoldata <- rbind(Sample = c("", "Sample", rownames(metadata)), qoldata)
  }
  qoldata <- data.frame(qoldata, stringsAsFactors = F, check.names = F)
  if(verbose) cat('Building matrix\n')
  edata <- as.matrix(edata); edata <- round(edata, 2)
  if(verbose) cat(' Adding gap\n')
  qmat <- cbind(Gene = rownames(edata), void = rep("", nrow(edata)), edata) # expression with gap
  if(verbose) cat(' Proper format: data frame\n')
  qmat <- data.frame(qmat, stringsAsFactors = F, check.names = F) # between Gene names and matrix
  if(verbose) cat('Binding\n')
  qformat <- rbind(qoldata, rep("", ncol(qmat)), c("Gene_ID", rep("", ncol(qmat)-1)), qmat)
  colnames(qformat) <- rownames(qformat) <- NULL
  tvar <- grep("Gene_ID", qformat[, 1])
  if(verbose) print(qformat[c(1:5, tvar:(tvar + 2)), 1:4])
  if(!is.null(fname)){
    tvar <- list(c(".csv", ","), c(".txt", "\t"))[[1]]
    fname <- paste0(fname, "_qlucore_format_R", length(rnames), "xC", length(cnames), tvar[1])
    fname <- sub("\\/{2,}", "/", fname)
    write.table(qformat, file = fname, sep = tvar[2], col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  return(qformat)
}
