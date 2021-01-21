#' @include loom.R
#'
NULL

#' Seurat binding for loom files
#'
#' @importFrom hdf5r h5attr list.groups
#'
#' @name loom-bindings
#' @rdname loom-bindings
#'
NULL

#' @importFrom Seurat DefaultAssay
#' @inheritParams Seurat::DefaultAssay
#'
#' @aliases DefaultAssay
#'
#' @rdname loom-bindings
#' @method DefaultAssay loom
#' @export
#'
DefaultAssay.loom <- function(object, ...) {
  if (Exists(x = object, name = H5Path('attrs', 'SEURAT_ASSAY'))) {
    return(object[['attrs']][['SEURAT_ASSAY']][])
  } else if (AttrExists(x = object, name = 'SEURAT_ASSAY')) {
    return(h5attr(x = object, which = 'SEURAT_ASSAY'))
  }
  warning(
    "Cannot find default assay information in this loom file",
    call. = FALSE,
    immediate. = TRUE
  )
  return(NULL)
}

#' @aliases dim
#'
#' @rdname loom-bindings
#' @method dim loom
#' @export
#'
dim.loom <- function(x) {
  return(rev(x = x[['matrix']]$dims))
}
