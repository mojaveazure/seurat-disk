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
#' @param src A source matrix
#' @param dest Destination HDF5 file or group for the padded matrix
#' @param dname Destination name for the padded matrix
#' @param dims A two-length integer vector with the number of rows and number of
#' columns in the padded matrix
#' @param index A two-length list of integer vectors describing the rows and
#' columns that \code{src} exists in
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
    # for (i in seq_along(along.with = index)) {
    #   if (is.numeric(x = index[[i]])) {
    #     if (all(index[[i]] == as.integer(x = index[[i]]))) {
    #       index[[i]] <- as.integer(x = index[[i]])
    #     }
    #     if (!is.integer(x = index[[i]]) || '') {
    #       ''
    #     }
    #   }
    # }
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
  definition = function(src, dest, dname, dims, index) {
    if (isFALSE(x = IsMatrix(x = src))) {
      stop("Not a matrix dataset", call. = FALSE)
    } else if (dest$exists(name = dname)) {
      stop("A dataset named ", dname, " already exists", call. = FALSE)
    }
    pmtx <- dest$create_dataset(
      name = dname,
      dtype = src$get_type(),
      space = H5S$new(dims = dims, maxdims = dims)
    )
    MARGIN <- GetMargin(dims = src$dims, MARGIN = 'largest')
    chunk.points <- ChunkPoints(
      dsize = src$dims[MARGIN],
      csize = src$chunk_dims[MARGIN]
    )
    na <- NA
    class(x = na) <- class(x = src[1, 1])
    na.write <- dims.write <- dims.read <- vector(mode = 'list', length = 2L)
    dims.read[[-MARGIN]] <- seq.default(from = 1, to = src$dims[-MARGIN])
    dims.write[[-MARGIN]] <- index[[-MARGIN]]
    # Write the NAs
    na.write[[MARGIN]] <- setdiff(
      x = seq.default(from = 1, to = dims[MARGIN]),
      y = index[[MARGIN]]
    )
    na.write[[-MARGIN]] <- setdiff(
      x = seq.default(from = 1, to = dims[-MARGIN]),
      y = index[[-MARGIN]]
    )
    for (i in seq_along(along.with = na.write)) {
      if (!length(x = na.write[[i]])) {
        next
      }
      na.dims <- vector(mode = 'list', length = 2L)
      na.dims[[i]] <- na.write[[i]]
      na.dims[[-i]] <- seq.default(from = 1, to = dims[-i])
      pmtx$write(
        args = na.dims,
        value = na
      )
    }
    # Write the data
    for (i in seq.default(from = 1, to = nrow(x = chunk.points))) {
      dims.read[[MARGIN]] <- seq.default(
        from = chunk.points[i, 'start'],
        to = chunk.points[i, 'end']
      )
      dims.write[[MARGIN]] <- index[[MARGIN]][dims.read[[MARGIN]]]
      pmtx$write(
        args = dims.write,
        value = src$read(args = dims.read)
      )
    }
    return(invisible(x = NULL))
  }
)

#' @rdname PadMatrix
#'
setMethod(
  f = 'PadMatrix',
  signature = c(src = "H5Group"),
  definition = function(src, dest, dname, dims, index) {
    .NotYetImplemented()
  }
)
