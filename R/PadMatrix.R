#' @include zzz.R
#' @include TestObject.R
#' @include TestH5.R
#' @importFrom methods setGeneric setMethod
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generic
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Pad a matrix
#'
#' @param ... ...
#'
#' @return ...
#'
#' @keywords internal
#'
setGeneric(
  name = 'PadMatrix',
  def = function(src, dest, dname, dims, index) {
    if (!inherits(x = dest, what = c('H5File', 'H5Group'))) {
      stop("'dest' must be an HDF5 file or group", call. = FALSE)
    } else if (!is.character(x = dname)) {
      stop("'dname' must be a character", call. = FALSE)
    }
    Writeable(x = dest)
    if (is.numeric(x = dims)) {
      if (all(dims == as.integer(x = dims))) {
        dims <- as.integer(x = dims)
      }
    }
    if (!is.integer(x = dims) || length(x = dims) != 2) {
      stop("'dims' must be a two-length integer value", call. = FALSE)
    }
    if (!(is.list(x = index) && !is.data.frame(x = index)) || length(x = index) != 2) {
      stop("'index' must be a two length list of integer values", call. = FALSE)
    }
    for (i in seq_along(along.with = index)) {
      if (is.numeric(x = index[[i]])) {
        if (all(index[[i]] == as.integer(x = index[[i]]))) {
          index[[i]] <- as.integer(x = index[[i]])
        }
        if (!is.integer(x = index[[i]]) || '') {
          ''
        }
      }
    }
    standardGeneric(f = 'PadMatrix')
  },
  signature = c('src')
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname PadMatrix
#'
setMethod(
  f = 'PadMatrix',
  signature = c(src = "H5D"),
  definition = function(src, dest, dname, dims) {
    if (isFALSE(x = IsMatrix(x = src))) {
      stop("Not a matrix dataset", call. = FALSE)
    } else if (dest$exists(name = dname)) {
      stop("A dataset named ", dname, " already exists", call. = FALSE)
    }
  }
)

#' @rdname PadMatrix
#'
setMethod(
  f = 'PadMatrix',
  signature = c(src = "H5Group"),
  definition = function(src, dest, dname, dims) {
    .NotYetImplemented()
  }
)
