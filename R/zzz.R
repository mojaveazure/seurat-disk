#' @importFrom rlang %||%
#' @importFrom methods setOldClass
#'
NULL

#' @docType package
#' @name SeuratDisk-pacakge
#' @rdname SeuratDisk-package
#'
# @section Package options:
#
# SeuratDisk uses the following options to control behaviour, users can configure
# these with \code{\link[base]{options}}:
#
# \describe{
#  \item{\code{SeuratDisk.loom.encoding}, \code{SeuratDisk.loom3.encoding}}{
#  Character encoding for loom files. The \code{loom3} encoding (default UTF-8)
#  is used for loom files >= v3.0.0; otherwise, the \code{loom} encoding
#  (default ASCII) is used}
#  \item{\code{SeuratDisk.loom.string_len}}{When the loom encoding is ASCII, set
#  the character width (default 7L)}
# }
#
# @aliases SeuratDisk
#
"_PACKAGE"


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

default.options <- list()

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global constants
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modes <- list(
  'new' = c('w', 'w-', 'x'),
  'existing' = c('r', 'r+')
)

version.regex <- '^\\d+(\\.\\d+){2}(\\.9\\d{3})?$'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal utility functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert a logical to an integer
#'
#' @param x A logical vector
#'
#' @return An integer vector
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::BoolToInt(x = c(TRUE, FALSE, NA))
#' }
#'
BoolToInt <- function(x) {
  x <- as.integer(x = x)
  x[which(x = is.na(x = x))] <- 2L
  return(x)
}

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
#' @examples
#' \donttest{
#' # Enumerate will use names if possible
#' x <- list(x = 1:3, y = letters[1:3], z = c('g1', 'g2', 'g3'))
#' SeuratDisk:::Enumerate(x = x)
#'
#' # If no (or missing) object names present, Enumerate will use the index
#' # number as the name
#' x <- unname(obj = x)
#' SeuratDisk:::Enumerate(x = x)
#' }
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
  class <- class[1]
  classdef <- getClass(Class = class)
  classpkg <- slot(object = classdef, name = 'package')
  if (classpkg %in% packages) {
    classpkg <- NULL
  }
  class <- paste(classpkg, class, sep = ':')
  return(gsub(pattern = '^:', replacement = '', x = class))
}

#' Guess an HDF5 Datatype
#'
#' Wrapper around \code{\link[hdf5r:guess_dtype]{hdf5r::guess_dtype}}, allowing
#' for the customization of string types rather than defaulting to
#' variable-length ASCII-encoded strings
#'
#' @inheritParams StringType
#' @inheritParams hdf5r::guess_dtype
#' @inheritDotParams hdf5r::guess_dtype
#'
#' @return An object of class \code{\link[hdf5r]{H5T}}
#'
#' @importFrom hdf5r guess_dtype H5T_COMPOUND
#'
#' @seealso \code{\link[hdf5r]{guess_dtype}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Characters can either be variable-width UTF8-encoded or
#' # fixed-width ASCII-encoded
#' SeuratDisk:::GuessDType(x = 'hello')
#' SeuratDisk:::GuessDType(x = 'hello', stype = 'ascii7')
#'
#' # Data frames are a compound type; character columns follow the same rules
#' # as character vectors
#' df <- data.frame(x = c('g1', 'g2', 'g3'), y = 1, 2, 3, stringsAsFactors = FALSE)
#' SeuratDisk:::GuessDType(x = df)
#' SeuratDisk:::GuessDType(x = df, stype = 'ascii7')
#'
#' # Logicals are turned into integers to ensure compatability with Python
#' # TRUE evaluates to 1, FALSE to 0, and NA to 2
#' SeuratDisk:::GuessDType(x = c(TRUE, FALSE, NA))
#' }
#'
GuessDType <- function(x, stype = 'utf8', ...) {
  dtype <- guess_dtype(x = x, ...)
  if (inherits(x = dtype, what = 'H5T_STRING')) {
    dtype <- StringType(stype = stype)
  } else if (inherits(x = dtype, what = 'H5T_COMPOUND')) {
    cpd.dtypes <- dtype$get_cpd_types()
    for (i in seq_along(along.with = cpd.dtypes)) {
      if (inherits(x = cpd.dtypes[[i]], what = 'H5T_STRING')) {
        cpd.dtypes[[i]] <- StringType(stype = stype)
      }
    }
    dtype <- H5T_COMPOUND$new(
      labels = dtype$get_cpd_labels(),
      dtypes = cpd.dtypes,
      size = dtype$get_size()
    )
  } else if (inherits(x = dtype, what = 'H5T_LOGICAL')) {
    dtype <- guess_dtype(x = BoolToInt(x = x), ...)
  }
  return(dtype)
}

#' Check the datatype of an HDF5 dataset
#'
#' @param x An HDF5 dataset (object of type \code{\link[hdf5r]{H5D}})
#' @param dtype A character vector of HDF5 datatype names, must be present in
#' \code{\link[hdf5r]{h5types}}
#'
#' @return A logical
#'
#' @importFrom hdf5r h5types
#'
#' @seealso \code{\link[hdf5r]{h5types}}
#'
#' @keywords internal
#'
IsDType <- function(x, dtype) {
  if (!inherits(x = x, what = 'H5D')) {
    stop("'IsDType' only works on HDF5 dataset", call. = FALSE)
  }
  dtypes <- unique(x = sapply(
    X = grep(pattern = '^H5T_', x = names(x = h5types), value = TRUE),
    FUN = function(i) {
      return(class(x = h5types[[i]])[1])
    },
    USE.NAMES = FALSE
  ))
  dtypes <- unique(x = c(dtypes, 'H5T_COMPOUND'))
  match.arg(arg = dtype, choices = dtypes, several.ok = TRUE)
  missing.dtypes <- setdiff(x = dtype, y = dtypes)
  if (length(x = missing.dtypes)) {
    dtype <- setdiff(x = dtype, y = missing.dtypes)
    if (!length(x = dtype)) {
      stop("None of the requested dtypes are valid HDF5 datatypes", call. = FALSE)
    } else {
      warning(
        "The following requested dtypes are not valid HDF5 datatypes: ",
        paste(missing.dtypes, sep = ", "),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  # dtype <- vapply(
  #   X = dtype,
  #   FUN = function(i) {
  #     return(class(x = h5types[[i]])[1])
  #   },
  #   FUN.VALUE = character(length = 1L),
  #   USE.NAMES = FALSE
  # )
  return(inherits(x = x$get_type(), what = dtype))
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
#' \donttest{
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

#' Make a space
#'
#' @param n Length space should be
#'
#' @return A space (' ') of length \code{n}
#'
#' @examples
#' \donttest{
#' SeuratDisk:::MakeSpace(n = 10)
#' cat('hello', SeuratDisk:::MakeSpace(n = 10), 'world\n', sep = '')
#' }
#'
#' @keywords internal
#'
MakeSpace <- function(n) {
  return(paste(rep_len(x = ' ', length.out = n), collapse = ''))
}

#' Generate an HDF5 string dtype
#'
#' @param stype Type of string encoding to use, choose from:
#' \describe{
#'  \item{utf8}{Variable-width, UTF-8}
#'  \item{ascii7}{Fixed-width (7 bits), ASCII}
#' }
#'
#' @return An \code{\link[hdf5r]{H5T_STRING}} object
#'
#' @importFrom hdf5r h5const H5T_STRING
#'
#' @keywords internal
#'
StringType <- function(stype = c('utf8', 'ascii7')) {
  stype <- match.arg(arg = stype)
  return(switch(
    EXPR = stype,
    'utf8' = H5T_STRING$new(size = Inf)$set_cset(cset = h5const$H5T_CSET_UTF8),
    'ascii7' = H5T_STRING$new(size = 7L)
  ))
}

#' Update slots in an object
#'
#' @param object An object to update
#'
#' @return \code{object} with the latest slot definitions
#'
#' @importFrom methods slotNames slot slot<-
#'
#' @keywords internal
#'
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(
        mode = class(x = xobj),
        length = 1L
      )
    }
  }
  return(object)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Loading handler
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
