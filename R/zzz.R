#' @importFrom hdf5r h5const
#' @importFrom methods setOldClass
#'
NULL

#' @docType package
#' @name SeuratDisk-pacakge
#' @rdname SeuratDisk-package
#'
#' @section Package options:
#'
#' SeuratDisk uses the following options to control behaviour, users can configure
#' these with \code{\link[base]{options}}:
#'
#' \describe{
#'  \item{\code{SeuratDisk.loom.encoding}, \code{SeuratDisk.loom3.encoding}}{
#'  Character encoding for loom files. The \code{loom3} encoding (default UTF-8)
#'  is used for loom files >= v3.0.0; otherwise, the \code{loom} encoding
#'  (default ASCII) is used}
#'  \item{\code{SeuratDisk.loom.string_len}}{When the loom encoding is ASCII, set
#'  the character width (default 7L)}
#' }
#'
#' @aliases SeuratDisk
#'
"_PACKAGE"

default.options <- list(
  SeuratDisk.loom.encoding = h5const$H5T_CSET_ASCII,
  SeuratDisk.loom3.encoding = h5const$H5T_CSET_UTF8,
  SeuratDisk.loom.string_len = 7L
)

#' Get a class string with package information
#'
#' @param class Class name
#' @param packages A vector of packages to exclude from resulting class information
#'
#' @return A character vector with the class
#'
#' @importFrom methods getClass slot
#'
#' @keywords internal
#'
GetClass <- function(class, packages = 'Seurat') {
  .NotYetImplemented()
  class <- class[1]
  classdef <- getClass(Class = class)
  class <- paste(slot(object = classdef, name = 'package'), class, sep = ':')
  return(class)
}

#' Check to see if a matrix is empty
#'
#' Determine if a matrix is empty or not. A matrix is considered empty if it
#' satisfies one of the following conditions:
#' \itemize{
#'  \item The dimensions of the matrix are 0-by-0 (\code{all(dim(x) == 0)})
#'  \item The dimensions of the matrix are 1-by-1 (\code{all(dim(x) == 1)}) and
#'  the sole vlaue is \code{NA}
#' }
#' These two situations correspond to matrices generated with either
#' \code{new('matrix')} or \code{matrix()}
#'
#' @param x A matrix
#'
#' @return \code{TRUE} if the matrix is empty otherwise \code{FALSE}
#'
#' @keywords internal
#'
#' @examples
#' IsMatrixEmpty(new('matrix'))
#' IsMatrixEmpty(matrix())
#' IsMatrixEmpty(matrix(1:9, nrow = 3))
#'
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

.onLoad <- function(libname, pkgname) {
  # Make the classes defined in SeuratDisk compatible with S4 generics/methods
  setOldClass(Classes = c('scdisk', 'h5Seurat', 'loom'))
  # Set some default options
  op <- options()
  toset <- !names(x = default.options) %in% names(x = op)
  if (any(toset)) {
    options(default.options[toset])
  }
  invisible(x = NULL)
}
