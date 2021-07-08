#!/usr/bin/R

suppressPackageStartupMessages({
  library(openxlsx)
})

## create and add a style to the column headers
body_style <- function(format_type){
  createStyle(
    fontSize = 12, fontColour = "#000000",
    numFmt = format_type,
    border = "TopBottomLeftRight", borderColour = "#000000", borderStyle = "thin",
    fgFill = "#FFFFFF",
    halign = "center", valign = "center"
  )
}
body_section_style <- function(side = "right"){
  createStyle(
    fontSize = 12, fontColour = "#000000",
    numFmt = "GENERAL",
    border = side, borderColour = "#000000", borderStyle = "thick",
    fgFill = "#FFFFFF",
    halign = "center", valign = "center"
  )
}
feature_style <- function(fill = "#D1D1D1"){
  createStyle(
    fontSize = 12, fontColour = "#000000",
    numFmt = "TEXT",
    border = c("left", "right"), borderColour = "#000000", borderStyle = "thick",
    fgFill = fill,
    halign = "center", valign = "center",
    textDecoration = "italic"
  )
}
caption_style <- createStyle(
  fontName = "Arial", fontSize = 12, fontColour = "#000000", numFmt = "TEXT",
  valign = "center", textDecoration = "bold"
)
header_style <- function(fill = "#BFBFBF"){
  createStyle(
    fontSize = 12, fontColour = "#000000",
    numFmt = "TEXT",
    border = "TopBottomLeftRight", borderColour = "#000000", borderStyle = "thick",
    fgFill = fill,
    halign = "center", valign = "center",
    textDecoration = "bold",
    locked = FALSE
  )
}

# sinew::makeOxygen(supp_table)

#' @title Supplementary table.
#'
#' @description `supp_table` creates sheets from a list of tables.
#'
#' @details This is a function that takes a [named] list of data.frames
#' and adds them to an XLSX file. You can specify how columns are grouped.
#'
#' @param mytables data.frame or list of data.frames (each will be a sheet/tab).
#' Names will be used to name each sheet/tab ("description~TableName" is valid).
#' @param headers Headers description, Default: NULL. Example:
#' list("Columns group name" = c("hange" = "NUMBER", "^pvalue" = "SCIENTIFIC")).
#' @param title_name This will be appended to each table description,
#' Default: 'Supplementary table #.
#' @param workfile Workbook object to add tables to, Default: NULL.
#' @param select_features Function (return rownames), Default: NULL.
#' @param highlight_features Column to highlight (genes names?),
#' Default: TRUE (highlights the first column in each table); accepts character.
#' @param color_header Background color for headers, Default: NULL.
#' @param color_features Background color for feature/row names, Default: NULL.
#' @param filename File name, Default: NULL.
#' @param rename_columns Make columns unique (grouped columns), Default: TRUE.
#' @param round_num Number of decimal places to round to, Default: 2.
#' @param body_start Row to format elements after header, Default: 4.
#' @param verbose Show progress, Default: TRUE.
#'
#' @keywords EXCEL formats
#'
#' @return A Workbook object or writes to filename.xlsx
#'
#' @author Ciro Ramírez-Suástegui, \email{ksuasteguic@gmail.com}
#' @references
#' @seealso
#'
#' @importFrom openxlsx createWorkbook addWorksheet writeData addStyle setColWidths setRowHeights saveWorkbook
#'
#' @examples
#' \dontrun{
#'   supp_table(
#'     mytables = iris[, 5:1],
#'     filename = "supp_table"
#'   )
#'   supp_table(
#'     mytables = list(Flowers = iris[, 5:1], Cars = cbind(Model = rownames(mtcars), mtcars)),
#'     filename = "supp_table_composed"
#'   )
#'   supp_table(
#'     mytables = iris[, 5:1],
#'     headers = list(c(pecie = "TEXT"), Sepal = c(epal = "NUMBER"), Petal = c(etal = "NUMBER")),
#'     filename = "supp_table_grouped"
#'   )
#'   supp_table(
#'     mytables = list(Flowers = iris[, 5:1], Cars = cbind(Model = rownames(mtcars), mtcars)),
#'     headers = list(
#'       c(pecie = "TEXT"), Sepal = c(epal = "NUMBER"), Petal = c(etal = "NUMBER"),
#'       c(odel = "TEXT"), Specifications = c("mpg|cyl|disp|hp|drat|wt|qsec|vs|am|gear|carb" = "NUMBER")
#'     ),
#'     filename = "supp_table_composed_grouped"
#'   )
#' }
#'
#' @export
#'

supp_table <- function(
  mytables,
  headers = NULL,
  title_name = "Supplementary table #.",
  workfile = NULL,
  select_features = NULL,
  highlight_features = TRUE,
  color_header = "#BFBFBF",
  color_features = "#D1D1D1",
  filename = NULL,
  rename_columns = TRUE,
  round_num = 2,
  body_start = 4,
  verbose = TRUE
){

  if(is.null(workfile)){
    workfile <- createWorkbook(
      title = title_name
    )
  }

  if(!is.list(mytables) || is.data.frame(mytables)) mytables <- list(mytables)
  if(is.null(names(mytables)))
    names(mytables) <- paste0("Table ", 1:length(mytables))
  names(mytables) <- make.unique(names(mytables))
  if(verbose > 1) str(mytables)

  if(is.null(headers)){
    headers <- unlist(lapply(unname(mytables), function(x) lapply(x, class) ))
    headers <- headers[!duplicated(names(headers))]
    headers <- ifelse(headers == "numeric", "NUMBER", "GENERAL")
    headers[headers %in% c("character", "factor")] <- "TEXT"
  }
  if(!is.list(headers)) headers <- list(none = headers)
  names(headers)[names(headers) == ""] <- "none"
  colnamestype <- unlist(headers, use.names = FALSE)
  names(colnamestype) <- unlist(sapply(headers, names), use.names = FALSE)
  if(any(names(colnamestype) == ""))
    stop("Format instruction (headers) without column pattern")

  for(tname in names(mytables)){
    header_group = body_start_i = body_start
    if(verbose) cat(tname, "\n")
    ddf <- mytables[[tname]]
    if(!is.null(select_features)) ddf <- ddf[select_features(x = ddf), ]
    if(!is.null(ddf$gene_name)) ddf$gene_name <- gsub("'", "", ddf$gene_name)
    headers2write <- lapply(
      X = headers, # taking names (patterns) for each element in the list
      FUN = function(p){ # another apply to follow the order
        ps = names(p)[!names(p) %in% c("rev", "reverse", "exclude")]
        y <- lapply(ps,
          function(z) grep(pattern = z, x = colnames(ddf), value = TRUE)
        )
        y <- unlist(y, use.names = FALSE)
        if(any(names(p) %in% c("rev", "reverse"))) y <- rev(y)
        if(any(names(p) %in% "exclude"))
          y <- y[!grepl(pattern = p[["exclude"]], x = y)]
        return(y)
      }
    )
    taken_cols <- unname(unlist(headers2write));
    if(verbose){
      cat("Columns:", length(taken_cols), "\n")
      print(format(taken_cols, justify = "centre", trim = TRUE))
    }; if(length(taken_cols) == 0) warning("No columns found in ", tname, "\n")
    ddf <- ddf[, taken_cols, drop = FALSE]

    tvar <- which(colnamestype == "NUMBER")
    if(length(tvar) > 0){
      if(verbose) cat("Rounding numbers to", round_num, "decimals ")
      tvar <- paste0(names(tvar), collapse = "|")
      tvar <- grep(tvar, colnames(ddf), value = TRUE)
      if(verbose) cat("- in", length(tvar), "columns:")
      if(verbose > 1) print(format(tvar, justify = "centre", trim = TRUE))
      ddf[, tvar] <- round(ddf[, tvar], round_num); if(verbose) cat("\n")
    }

    tvar <- rename_columns
    if(isTRUE(tvar) || is.character(tvar) || is.list(tvar)){
      if(verbose) cat("Renaming columns\n")
      taken_new <- if(isTRUE(rename_columns)){ # parse names
        unname(unlist(lapply(headers2write, function(x){
          y <- sapply(strsplit(x = x, split = "_|\\."), c)
          # check which names are different across columns
          if(is.null(dim(y))) return(x)
          tvar <- apply(y, 1, function(z) length(unique(z)) ) == length(x)
          if(sum(tvar) == 0)
            tvar <- apply(y, 1, function(z) length(unique(z)) ) != length(x)
          apply(y[tvar, , drop = FALSE], 2, function(z) paste0(z, collapse = " ") )
        })))
      }else{
        if(!is.list(tvar)) rename_columns <- list(rename_columns)
        taken_new <- colnames(ddf)
        for(i in 1:length(rename_columns)){
          tvar <- ifelse(is.na(rename_columns[[i]][2]), "", rename_columns[[i]][2])
          taken_new <- gsub(rename_columns[[i]][1], tvar, taken_new)
        }; taken_new
      }; # taken_new <- make.unique(taken_new)
      if(verbose > 1) str(ddf); colnames(ddf) <- taken_new
    }; if(verbose) str(ddf)

    tabname = gsub(".*~", "", tname)
    if(verbose) cat("Tab name:", tabname, "\n")
    addWorksheet(wb = workfile, sheetName = tabname)
    writeData(
      wb = workfile, sheet = tabname,
      x = paste(title_name, gsub("~.*", "", tname)))
    heads <- sapply(headers2write, length);
    heads <- rep(names(heads), heads)
    header_group <- if(any(!heads == "none")){
      body_start_i-1
    }else body_start_i <- body_start_i-1
    if(verbose) cat("Merging headers titles\n")
    for(headname in unique(heads[heads!="none"])){
      headcols <- which(heads == headname)
      if(verbose){
        cat(" - ", heads[headcols[1]], ": ", min(headcols),
        " to ", max(headcols), "\n", sep = "")
      }
      writeData(
        wb = workfile, sheet = tabname, x = headname,
        startCol = min(headcols), startRow = header_group)
      mergeCells(wb = workfile, sheet = tabname,
        cols = headcols, rows = header_group)
    }
    if(verbose) cat("Starting:", body_start_i, "\n")
    writeData(
      wb = workfile,
      sheet = tabname,
      x = ddf,
      startRow = body_start_i
    )

    if(verbose) cat("Styles\n") # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hfeatures = if(is.character(highlight_features)){
      which(taken_cols %in% highlight_features)
    }else if(isTRUE(highlight_features)) 1 else highlight_features
    if(verbose) cat(" - title/caption\n")
    addStyle(wb = workfile, sheet = tabname,
      style = caption_style, rows = 1, cols = 1)
    if(verbose) cat(" - headers\n")
    for(taken in taken_cols){
      if(verbose > 1) cat("  *", taken, "> ")
      tvar <- which(sapply(headers2write, function(x) any(grepl(taken, x)) ))
      if(length(tvar) > 1){ print(tvar); print(headers2write)}
      tvar <- names(which(sapply(names(headers[[tvar]]), grepl, taken)))
      # tvar <- which(sapply(names(colnamestype), grepl, taken))
      # myformat <- colnamestype[[tvar]]
      myformat <- colnamestype[[which(names(colnamestype) %in% tvar)]]
      this_header = which(taken_cols %in% taken)
      if(verbose > 1) cat(myformat, "\n")
      addStyle(
        wb = workfile, sheet = tabname,
        style = body_style(format_type = myformat),
        rows = body_start_i:(nrow(ddf)+body_start_i),
        cols = this_header, gridExpand = TRUE
      ); tvar <- this_header %in% hfeatures && heads[this_header] == "none"
      tvar <- if(tvar) body_start_i else header_group
      addStyle(
        wb = workfile, sheet = tabname, style = header_style(fill = color_header),
        rows = tvar:body_start_i, cols = this_header, gridExpand = TRUE
      )
    }; if(verbose) cat(" - body by section\n")
    for(headname in heads){
      addStyle(
        wb = workfile, sheet = tabname, style = body_section_style(),
        rows = (body_start_i+1):(nrow(ddf)+body_start_i),
        cols = max(which(heads == headname)),
        gridExpand = TRUE, stack = TRUE
      )
    }; if(verbose) cat(" - body bottom\n")
    addStyle(
      wb = workfile, sheet = tabname, style = body_section_style(side = "bottom"),
      rows = (nrow(ddf)+body_start_i), cols = 1:length(taken_cols),
      gridExpand = TRUE, stack = TRUE
    )
    if(is.numeric(hfeatures)){
      if(verbose) cat(" - features (first column's rows)\n")
      addStyle(
        wb = workfile, sheet = tabname, style = feature_style(fill = color_features),
        rows = body_start_i:(nrow(ddf)+body_start_i), cols = hfeatures,
        gridExpand = TRUE, stack = TRUE
      )
    }

    if(verbose) cat("Setting widths\n") # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    setColWidths(wb = workfile, sheet = tabname, cols = 1:ncol(ddf), widths = 16)
    if(any(!heads == "none"))
      setRowHeights(wb = workfile, sheet = tabname, rows = header_group, heights = 36)
    setRowHeights(wb = workfile, sheet = tabname, rows = body_start_i, heights = 18)
  }
  if(!is.null(filename)){
    fname <- paste0(filename, ifelse(grepl("xlsx$", filename), "", ".xlsx"))
    if(verbose) cat("Writing to", fname, "\n")
    saveWorkbook(wb = workfile, file = fname, overwrite = TRUE)
  }else{ return(workfile) }
}
