#' @include TestObject.R
#' @include TestH5.R
#'
NULL

#' Check to see if a dataset, group, or attribute exists in an HDF5 file,
#' group, or dataset
#'
#' @param x An HDF5 \link[hdf5r:H5File]{file} or \link[hdf5r:H5Group]{group};
#' for \code{AttrExists}, may also be a \link[hdf5r:H5D]{dataset}
#' @param name Name of dataset, group, or attribute to test for
#'
#' @return \code{TRUE} if \code{name} exists in \code{x}, otherwise \code{FALSE}
#'
#' @name H5Exists
#' @rdname H5Exists
#'
#' @keywords internal
#'
AttrExists <- function(x, name) {
  if (!inherits(x = x, what = c('H5File', 'H5Group', 'H5D'))) {
    stop("'x' must be an HDF5 file, group, or dataset", call. = FALSE)
  }
  exists <- x$attr_exists(attr_name = name)
  if (isTRUE(x = exists)) {
    space <- x$attr_open(attr_name = name)$get_space()
    exists <- !(length(x = space$dims) > 0 && space$dims == 0)
  }
  return(exists)
}

#' Get the dimensions of an HDF5 dataset or sparse matrix
#'
#' @param x An HDF5 dataset or sparse matrix
#'
#' @return A vector with the dimensions of the dataset or sparse matrix. For
#' sparse matrices, if no dimensions are found in either the \dQuote{dims} or
#' \dQuote{shape} attribute, returns \code{c(NA_integer_, NA_integer_)}
#'
#' @seealso \code{\href{../doc/h5Seurat-spec.html}{vignette("h5Seurat-spec")}}
#'
#' @keywords internal
#'
Dims <- function(x) {
  UseMethod(generic = 'Dims', object = x)
}

#' @method Dims H5D
#'
Dims.H5D <- function(x) {
  return(x$dims)
}

#' @importFrom hdf5r h5attr
#'
#' @method Dims H5Group
#'
Dims.H5Group <- function(x) {
  if (!IsMatrix(x = x)) {
    stop("Not a matrix", call. = FALSE)
  }
  if (x$attr_exists(attr_name = 'dims')) {
    return(h5attr(x = x, which = 'dims'))
  } else if (x$attr_exists(attr_name = 'shape')) {
    return(rev(x = h5attr(x = x, which = 'shape')))
  }
  return(c(NA_integer_, NA_integer_))
}

#' @rdname H5Exists
#'
#' @keywords internal
#'
Exists <- function(x, name) {
  if (!inherits(x = x, what = c('H5File', 'H5Group'))) {
    stop("'x' must be an HDF5 file or group", call. = FALSE)
  }
  name <- unlist(x = strsplit(x = name[1], split = '/', fixed = TRUE))
  path <- character(length = 1L)
  exists <- TRUE
  for (i in seq_along(along.with = name)) {
    path <- paste(path, name[i], sep = '/')
    exists <- x$exists(name = path)
    if (isFALSE(x = exists)) {
      break
    }
  }
  return(exists)
}

#' Create an HDF5 object path
#'
#' @inheritParams base::paste
#'
#' @return A character vector with path ready for accessing data in an HDF5 file
#' or group
#'
#' @keywords internal
#'
H5Path <- function(..., collapse = NULL) {
  return(paste(..., sep = '/', collapse = collapse))
}

#' Is an HDF5 file or group writeable
#'
#' @param x An \code{\link[hdf5r]{H5File}} or \code{\link[hdf5r]{H5Group}}
#' object
#' @param error Throw an error
#'
#' @return ...
#'
#' @keywords internal
#'
Writeable <- function(x, error = TRUE) {
  mode <- as.character(x = x$get_file_id()$get_intent()) == 'H5F_ACC_RDONLY'
  if (isTRUE(x = mode)) {
    msg <- paste("File", x$get_filename(), "is not writeable")
    if (isTRUE(x = error)) {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, immediate. = TRUE, call. = FALSE)
    }
  }
  return(invisible(x = !mode))
}
