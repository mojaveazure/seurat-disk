#' Seurat bindings for h5Seurat files
#'
#' @importFrom hdf5r h5attr list.groups
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
  object$attr_delete(attr_name = 'active.assay')
  object$create_attr(
    attr_name = 'active.assay',
    robj = value,
    dtype = GuessDType(x = value)
  )
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
  return(as.factor(x = object[['active.ident']]))
}

#' @importFrom Seurat IsGlobal
#' @inheritParams Seurat::IsGlobal
#'
#' @aliases IsGlobal
#'
#' @rdname h5Seurat-bindings
#' @method IsGlobal H5Group
#' @export
#'
IsGlobal.H5Group <- function(object) {
  return(
    object$attr_exists(attr_name = 'global') &&
      h5attr(x = object, which = 'global') == 1
  )
}

#' @importFrom Seurat Key
#' @inheritParams Seurat::Key
#'
#' @aliases Key
#'
#' @rdname h5Seurat-bindings
#' @method Key H5Group
#' @export
#'
Key.H5Group <- function(object, ...) {
  return(h5attr(x = object, which = 'key'))
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
  object$attr_delete(attr_name = 'project')
  object$create_attr(
    attr_name = 'project',
    robj = value,
    dtype = GuessDType(x = value)
  )
  return(invisible(x = object))
}

#' @importFrom Seurat Stdev
#' @inheritParams Seurat::Stdev
#'
#' @rdname h5Seurat-bindings
#' @method Stdev h5Seurat
#' @export
#'
Stdev.h5Seurat <- function(object, reduction = 'pca') {
  if (object[['reductions']]$exists(name = reduction)) {
    reduc.group <- object[['reductions']][[reduction]]
    if (reduc.group$exists(name = 'stdev')) {
      return(reduc.group[['stdev']][])
    }
  }
  return(numeric(length = 0L))
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
