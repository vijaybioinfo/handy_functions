#!/usr/bin/R

##################
# Seurat to loom #
##################

# This script transforms the Seurat object to loom format

#### Check arguments ####
if( grepl("3.5", sub(".*sion (.*) \\(.*", "\\1",version$version.string)) )
  .libPaths('~/R/newer_packs_library/3.5/')
library(optparse)
optlist <- list(
  make_option(
    opt_str = c("-m", "--mycellsf"), type = "character",
    help = "Object: the rdata object full path name."
  )
)
# Getting arguments from command line and setting their values to their respective variable names.
optparse <- OptionParser(option_list = optlist)
defargs <- parse_args(optparse)

setwd(dirname(defargs$mycellsf))
cat('Workingi in', getwd(), '\n')
cat('Loading object\n'); load(defargs$mycellsf)
mycells@reductions <- list()
mycells@graphs <- list()
cat('Transforming it to loom\n');
newname <- sub(".rdata", ".loom", defargs$mycellsf, ignore.case = TRUE)
if(!file.exists(newname)){
  mycellsloom <- Seurat::as.loom(mycells, filename = newname)
  mycellsloom$close_all()
}
cat("Done\n")
