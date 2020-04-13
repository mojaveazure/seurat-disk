#' @include zzz.R
#' @include TestObject.R
#' @include TestH5.R
#' @include sparse_matrix.R
#'
NULL

#' Transpose a matrix
#'
#' @param x A matrix to transpose
#' @param ... Arguments passed to other methods
#'
#' @name Transpose
#' @rdname Transpose
#'
#' @export
#'
Transpose <- function(x, ...) {
  UseMethod(generic = 'Transpose', object = x)
}

#' @importFrom methods slot
#' @importFrom Matrix sparseMatrix
#'
#' @return \code{\link[Matrix]{dgCMatrix}} method: returns a
#' \code{\link[Matrix]{dgCMatrix}} with the data of \code{x} transposed
#'
#' @rdname Transpose
#' @method Transpose dgCMatrix
#' @export
#'
Transpose.dgCMatrix <- function(x, ...) {
  i.order <- order(slot(object = x, name = 'i'))
  return(sparseMatrix(
    i = PointerToIndex(p = slot(object = x, name = 'p'))[i.order],
    p = IndexToPointer(j = slot(object = x, name = 'i') + 1),
    x = slot(object = x, name = 'x')[i.order],
    dims = rev(x = dim(x = x)),
    dimnames = rev(x = dimnames(x = x)),
    giveCsparse = TRUE
  ))
}

#' @param dest ...
#' @param dname ...
#' @param overwrite ...
#' @param verbose Show progress updates
#'
#' @return \code{\link[hdf5r]{H5D}} and \code{\link[hdf5r]{H5Group}} methods:
#' Invisibly returns \code{NULL}
#'
#' @importFrom hdf5r H5S
#'
#' @rdname Transpose
#' @method Transpose H5D
#' @export
#'
Transpose.H5D <- function(
  x,
  dest = GetParent(x = x),
  dname = paste0('t_', basename(path = x$get_obj_name())),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  if (isFALSE(x = IsMatrix(x = x))) {
    stop("Not a matrix dataset", call. = FALSE)
  } else if (!inherits(x = dest, what = c('H5File', 'H5Group'))) {
    stop("'dest' must be an HDF5 file or group", call. = FALSE)
  }
  Writeable(x = dest)
  if (isTRUE(x = dest$exists(name = dname))) {
    if (isTRUE(x = overwrite)) {
      dest$link_delete(name = dname)
    } else {
      stop("'", dname, "' already exists", call. = FALSE)
    }
  }
  if (isTRUE(x = verbose)) {
    pb <- PB()
  }
  mtx.dims <- rev(x = x$dims)
  tmtx <- dest$create_dataset(
    name = dname,
    dtype = x$get_type(),
    space = H5S$new(dims = mtx.dims, maxdims = mtx.dims)
  )
  MARGIN <- GetMargin(dims = x$dims, MARGIN = 'largest')
  chunk.points <- ChunkPoints(
    dsize = x$dims[MARGIN],
    csize = x$chunk_dims[MARGIN]
  )
  dims <- vector(mode = 'list', length = 2L)
  dims[[-MARGIN]] <- seq.default(from = 1, to = x$dims[-MARGIN])
  for (i in seq.default(from = 1, to = nrow(x = chunk.points))) {
    dims[[MARGIN]] <- seq.default(
      from = chunk.points[i, 'start'],
      to = chunk.points[i, 'end']
    )
    tmtx$write(
      args = dims[c(MARGIN, c(1, 2)[-MARGIN])],
      value = t(x = x$read(args = dims))
    )
    if (verbose) {
      setTxtProgressBar(pb = pb, value = i / nrow(x = chunk.points))
    }
  }
  if (isTRUE(x = verbose)) {
    close(con = pb)
  }
  return(invisible(x = NULL))
}

#' @rdname Transpose
#' @method Transpose H5Group
#' @export
#'
Transpose.H5Group <- function(
  x,
  dest = GetParent(x = x),
  dname = paste0('t_', basename(path = x$get_obj_name())),
  overwrite = FALSE,
  ...
) {
  if (!IsMatrix(x = x)) {
    stop("Not a sparse matrix", call. = FALSE)
  } else if (!inherits(x = dest, what = c('H5File', 'H5Group'))) {
    stop("'dest' must be an HDF5 file or group", call. = FALSE)
  } else if (as.character(x = dest$get_file_id()$get_intent()) == 'H5F_ACC_RDONLY') {
    stop("read only", call. = FALSE)
  }
  if (dest$exists(name = dname)) {
    if (overwrite) {
      dest$link_delete(name = dname)
    } else {
      stop("exists", call. = FALSE)
    }
  }
  nfile <- dest$create_group(name = dname)
  i.order <- order(x[['indices']]$read())
  nfile$create_dataset(
    name = 'indices',
    robj = PointerToIndex(p = x[['indptr']]$read())[i.order] - 1,
    dtype = GuessDType(x = integer(length = 1L))
  )
  nfile$create_dataset(
    name = 'indptr',
    robj = IndexToPointer(j = x[['indices']]$read() + 1),
    dtype = GuessDType(x = integer(length = 1L))
  )
  nfile$create_dataset(
    name = 'data',
    robj = x[['data']]$read()[i.order],
    dtype = GuessDType(x = x[['data']][1])
  )
  dims <- vapply(
    X = c('dims', 'h5sparse_shape'),
    FUN = x$attr_exists,
    FUN.VALUE = logical(length = 1L),
    USE.NAMES = TRUE
  )
  if (any(dims)) {
    dims <- names(x = which(x = dims))[1]
    obj.dims <- x$attr_open(attr_name = dims)$read()
    nfile$create_attr(
      attr_name = dims,
      robj = rev(x = obj.dims),
      dtype = GuessDType(x = obj.dims)
    )
  }
  return(invisible(x = NULL))
}
