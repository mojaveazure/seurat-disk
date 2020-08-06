#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Loom-file Loading
#'
#' Load data from a loom file into a \code{\link[Seurat]{Seurat}} object
#'
#' @param file,x Name of loom file or a \code{loom} object to load data from
#' @param assay Name of assay to store expression data as; if \code{NULL}, will
#' search for an HDF5 attribute named \code{SEURAT_ASSAY} or an attribute
#' dataset named \code{/attrs/SEURAT_ASSAY} for assay name. If not found,
#' defaults to \dQuote{RNA}
#' @param cells Name of dataset in \code{/col_attrs} with cell names
#' @param features Name of dataset in \code{/row_attrs} with feature names
#' @param normalized Name of matrix in \code{/layers} to store normalized data
#' as; pass \dQuote{/matrix} to store \code{/matrix} as normalized data instead
#' of raw counts
#' @param scaled Name of dataset in \code{/layers} to store scaled data as
#' @param filter Keep only selected cells and/or features as specified by
#' \code{/col_attrs/Valid} and \code{/row_attrs/Valid}, respectively
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return A \code{\link[Seurat]{Seurat}} object
#'
#' @name LoadLoom
#' @rdname LoadLoom
#'
#' @inherit LoomLoading details
#'
#' @inheritSection LoomLoading Loom 0.1 Loading
#'
#' @inheritSection LoomLoading Loom 3.0.0 Loading
#'
#' @seealso
#' \href{http://linnarssonlab.org/loompy/conventions/index.html}{Loom
#' file conventions}
#'
#' @export
#'
LoadLoom <- function(
  file,
  assay = NULL,
  cells = 'CellID',
  features = 'Gene',
  normalized = NULL,
  scaled = NULL,
  filter = c('cells', 'features', 'all', 'none'),
  verbose = TRUE,
  ...
) {
  UseMethod(generic = 'LoadLoom', object = file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LoadLoom
#' @method LoadLoom character
#' @export
#'
LoadLoom.character <- function(file, ...) {
  loom <- loom$new(filename = file, mode = 'r')
  on.exit(expr = loom$close_all())
  return(LoadLoom(file = loom, ...))
}

#' @rdname LoadLoom
#' @method LoadLoom H5File
#' @export
#'
LoadLoom.H5File <- function(file, ...) {
  return(LoadLoom(
    file = as.loom(x = file),
    ...
  ))
}

#' @rdname LoadLoom
#' @method LoadLoom loom
#' @export
#'
LoadLoom.loom <- function(file, ...) {
  return(as.Seurat(x = file, ...))
}

#' @aliases as.Seurat
#'
#' @rdname LoadLoom
#' @method as.Seurat loom
#' @export
#'
as.Seurat.loom <- function(
  x,
  assay = NULL,
  cells = 'CellID',
  features = 'Gene',
  normalized = NULL,
  scaled = NULL,
  filter = c('cells', 'features', 'all', 'none'),
  verbose = TRUE,
  ...
) {
  version <- ClosestVersion(query = x$version(), targets = c('0.1.0', '3.0.0'))
  load.fxn <- switch(
    EXPR = version,
    '0.1.0' = LoadLoom0.1,
    '0.3.0' = LoadLoom3.0
  )
  object <- load.fxn(
      file = x,
      assay = assay,
      cells = cells,
      features = features,
      normalized = normalized,
      scaled = scaled,
      filter = filter,
      verbose = verbose
    )
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Loom-file Loading
#'
#' Version-specific loom-file loading functions
#'
#' @inheritParams LoadLoom
#'
#' @inherit LoadLoom return
#'
#' @name LoomLoading
#' @rdname LoomLoading
#'
#' @details
#' \code{LoadLoom} will try to automatically fill slots of a \code{Seurat}
#' object based on data presense or absence in a given loom file. This method
#' varies by loom specification version. For version-specific details, see
#' sections below
#'
#' @section Loom 0.1 Loading:
#' Loading data from loom files less than version 3.0.0 is not
#' currently supported
#'
#' @keywords internal
#'
LoadLoom0.1 <- function(
  file,
  assay = NULL,
  cells = 'CellID',
  features = 'Gene',
  normalized = NULL,
  scaled = NULL,
  filter = c('cells', 'features', 'all', 'none'),
  verbose = TRUE
) {
  .NotYetImplemented()
}

#' @name LoomLoading
#' @rdname LoomLoading
#'
#' @section Loom 3.0.0 Loading:
#' blah
#'
LoadLoom3.0 <- function(
  file,
  assay = NULL,
  cells = 'CellID',
  features = 'Gene',
  normalized = NULL,
  scaled = NULL,
  filter = c('cells', 'features', 'all', 'none'),
  verbose = TRUE
) {
  .NotYetImplemented()
}
