#' Seurat bindings for h5Seurat files
#'
#' @importFrom hdf5r h5attr h5attr<- list.groups
#'
#' @name h5Seurat-bindings
#' @rdname h5Seurat-bindings
#'
NULL

#' @importFrom Seurat Cells
#' @inheritParams Seurat::Cells
#'
#' @aliases Cells
#'
#' @rdname h5Seurat-bindings
#' @method Cells h5Seurat
#' @export
#'
Cells.h5Seurat <- function(x) {
  if (!x$exists(name = 'cell.names')) {
    stop("Cannot find cell names in this h5Seurat file", call. = FALSE)
  }
  return(x[['cell.names']][])
}

#' @importFrom Seurat DefaultAssay
#' @inheritParams Seurat::DefaultAssay
#'
#' @aliases DefaultAssay
#'
#' @rdname h5Seurat-bindings
#' @method DefaultAssay h5Seurat
#' @export
#'
DefaultAssay.h5Seurat <- function(object, ...) {
  return(h5attr(x = object, which = 'active.assay'))
}

#' @importFrom Seurat DefaultAssay<-
#'
#' @rdname h5Seurat-bindings
#' @method DefaultAssay<- h5Seurat
#' @export
#'
"DefaultAssay<-.h5Seurat" <- function(object, ..., value) {
  if (!value %in% list.groups(object = object[['assays']])) {
    stop("", call. = FALSE)
  }
  h5attr(x = object, which = 'active.assay') <- value
  return(invisible(x = object))
}

#' @importFrom Seurat Idents
#' @inheritParams Seurat::Idents
#'
#' @aliases Idents
#'
#' @rdname h5Seurat-bindings
#' @method Idents h5Seurat
#' @export
#'
Idents.h5Seurat <- function(object, ...) {
  .NotYetImplemented()
}

#' @importFrom Seurat Project
#' @inheritParams Seurat::Project
#'
#' @aliases Project
#'
#' @rdname h5Seurat-bindings
#' @method Project h5Seurat
#' @export
#'
Project.h5Seurat <- function(object, ...) {
  return(h5attr(x = object, which = 'project'))
}

#' @importFrom Seurat Project<-
#'
#' @rdname h5Seurat-bindings
#' @method Project<- h5Seurat
#' @export
#'
"Project<-.h5Seurat" <- function(object, ..., value) {
  h5attr(x = object, which = 'project') <- value
  return(invisible(x = object))
}

# levels
# levels<-
# Loadings
# Misc
# Misc<-
# Project
# Project<-
# ReorderIdent
# RenameCells
# RenameIdents
# SetAssayData
# SetIdent
# StashIdent
# Stdev
# Tool
# Tool<-
# VariableFeatures
# VariableFeatures<-
# WhichCells
