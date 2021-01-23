#' @include TestObject.R
#' @importFrom methods setOldClass setMethod
#'
NULL

setOldClass(Classes = 'H5D')
setOldClass(Classes = 'H5Group')

#' Test HDF5 datasets and groups to see what kind of data they are
#'
#' @param x An HDF5 dataset or group
#'
#' @name TestH5
#' @rdname TestH5
#'
#' @seealso \code{\link{TestObject}}
#'
#' @keywords internal
#'
NULL

#' @rdname TestH5
#'
setMethod(
  f = 'IsDataFrame',
  signature = c('x' = 'H5D'),
  definition = function(x) {
    return(inherits(x = x$get_type(), what = 'H5T_COMPOUND'))
  }
)

#' @rdname TestH5
#'
setMethod(
  f = 'IsDataFrame',
  signature = c('x' = 'H5Group'),
  definition = function(x) {
    if (AttrExists(x = x, name = 's4class')) {
      return(FALSE)
    }
    check <- sapply(
      X = names(x = x),
      FUN = function(i) {
        if (inherits(x = x[[i]], what = 'H5Group')) {
          if (IsFactor(x = x[[i]])) {
            return(x[[i]][['values']]$dims)
          }
        } else {
          return(ifelse(test = x[[i]]$dims == 1, yes = NA_real_, no = x[[i]]$dims))
        }
        return(NA_real_)
      },
      USE.NAMES = FALSE
    )
    return(!any(is.na(x = check)) && length(x = unique(x = check)) == 1)
  }
)

#' @rdname TestH5
#'
setMethod(
  f = 'IsFactor',
  signature = c('x' = 'H5D'),
  definition = function(x) {
    return(FALSE)
  }
)

#' @rdname TestH5
#'
setMethod(
  f = 'IsFactor',
  signature = c('x' = 'H5Group'),
  definition = function(x) {
    check <- x$exists(name = 'levels') &&
      inherits(x = x[['levels']], what = 'H5D') &&
      IsDType(x = x[['levels']], dtype = 'H5T_STRING')
    if (isTRUE(x = check)) {
      check <- x$exists(name = 'values') &&
        inherits(x = x[['values']], what = 'H5D') &&
        IsDType(x = x[['values']], dtype = 'H5T_INTEGER')
    }
    if (all(check)) {
      check <- vapply(
        X = c('levels', 'values'),
        FUN = function(i) {
          return(length(x = x[[i]]$dims) == 1)
        },
        FUN.VALUE = logical(length = 1L)
      )
    }
    if (all(check)) {
      check <- length(x = na.omit(object = unique(x = x[['values']]$read()))) <= x[['levels']]$dims
    }
    return(all(check))
  }
)

#' @rdname TestH5
#'
setMethod(
  f = 'IsList',
  signature = c('x' = 'H5Group'),
  definition = function(x) {
    return(!IsDataFrame(x = x) && !IsFactor(x = x) && !IsMatrix(x = x))
  }
)

#' @rdname TestH5
#'
setMethod(
  f = 'IsLogical',
  signature = c('x' = 'H5D'),
  definition = function(x) {
    check <- AttrExists(x = x, name = 's3class')
    if (isTRUE(x = check)) {
      check <- x$attr_open('s3class')$read() == 'logical'
    }
    if (isTRUE(x = check)) {
      check <- IsDType(x = x, dtype = 'H5T_INTEGER')
    }
    if (isTRUE(x = check)) {
      # check <- length(x = Dims(x = x)) == 1
      check <- length(x = x$dims) == 1
    }
    if (isTRUE(x = check)) {
      check <- all(unique(x = x$read()) %in% seq.int(from = 0, to = 2))
    }
    return(check)
  }
)

#' @rdname TestH5
#'
setMethod(
  f = 'IsMatrix',
  signature = c('x' = 'H5D'),
  definition = function(x) {
    return(length(x = x$dims) == 2)
  }
)

#' @rdname TestH5
#'
setMethod(
  f = 'IsMatrix',
  signature = c('x' = 'H5Group'),
  definition = function(x) {
    slots <- c('indices', 'indptr', 'data')
    check <- vapply(
      X = slots,
      FUN = function(i) {
        return(x$exists(name = i) && inherits(x = x[[i]], what = 'H5D'))
      },
      FUN.VALUE = logical(length = 1L),
      USE.NAMES = FALSE
    )
    return(all(check))
  }
)
