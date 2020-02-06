#' @importFrom rlang %||%
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

modes <- list(
  'new' = c('w', 'w-', 'x'),
  'existing' = c('r', 'r+')
)

#' Find the closest version
#'
#' @param query A query version (\code{character} or \code{numeric_version})
#' @param targets A vector of target versions (\code{character} or
#' \code{numeric_version})
#' @param direction Which way should we check for closest version? Choose from:
#' \describe{
#'  \item{min}{Closest version less than or equal to \code{query}}
#'  \item{max}{Closest version greater than or equal to \code{query}}
#' }
#'
#' @return The version from \code{targets} that is closest to \code{query} as a
#' \code{character} vector
#'
#' @keywords internal
#'
ClosestVersion <- function(query, targets, direction = c('min', 'max')) {
  direction <- match.arg(arg = direction)
  query <- numeric_version(x = query)
  targets <- sort(x = numeric_version(x = targets))
  index <- suppressWarnings(expr = switch(
    EXPR = direction,
    'min' = max(which(x = targets <= query)),
    'max' = min(which(x = targets >= query))
  ))
  if (is.infinite(x = index)) {
    stop(
      "All target versions ",
      switch(EXPR = direction, 'min' = 'greater', 'max' = 'less'),
      " than query version (",
      as.character(x = query),
      ")",
      call. = FALSE
    )
  }
  return(as.character(x = targets[index]))
}

#' Enumerate a list or vector
#'
#' @param x A list or a vector
#'
#' @return A list of length \code{x} with the following named values:
#' \describe{
#'   \item{\code{name}}{The name of \code{x} at a given index}
#'   \item{\code{value}}{The value of \code{x} at the corresponding index}
#' }
#'
#' @note For any given index \code{i} in \code{x}, all attempts to use the name
#' of the value of \code{x} at \code{i} will be made. If there is no name
#' (eg. \code{nchar(x = names(x = x)[i]) == 0}), the index \code{i} will be used
#' in its stead
#'
#' @keywords internal
#'
Enumerate <- function(x) {
  indices <- seq_along(along.with = x)
  keys <- names(x = x) %||% as.character(x = indices)
  keys[nchar(x = keys) == 0] <- indices[nchar(x = keys) == 0]
  vals <- lapply(
    X = indices,
    FUN = function(i) {
      return(c('name' = keys[i], 'value' = unname(obj = x[i])))
    }
  )
  return(vals)
}

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
#' \dontrun{
#' SeuratDisk:::IsMatrixEmpty(new('matrix'))
#' SeuratDisk:::IsMatrixEmpty(matrix())
#' SeuratDisk:::IsMatrixEmpty(matrix(1:9, nrow = 3))
#' }
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
