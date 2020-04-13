#' Convert sparse matrix pointers to indices and vice versa
#'
#' @param j A vector of sparse matrix colum indices
#'
#' @return \code{IndexToPointer}: A vector of index pointers (p)
#'
#' @name SparsePointers
#' @rdname SparsePointers
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' dat <- dat <- c(0, 0, 1, 4, 0, 2, 0, 9, 0)
#' smat <- Matrix::Matrix(data = dat, nrow = 3, sparse = TRUE)
#' j <- SeuratDisk:::PointerToIndex(p = smat@p)
#' Matrix::sparseMatrix(i = smat@i + 1, j = j, x = smat@x)
#' p <- SeuratDisk:::IndexToPointer(j = j)
#' Matrix::sparseMatrix(i = smat@i + 1, p = p, x= smat@x)
#' }
IndexToPointer <- function(j) {
  p <- vector(mode = 'integer', length = max(j) + 1)
  index <- seq.int(from = 2, to = length(x = p))
  for (i in seq_along(along.with = index)) {
    p[index[i]] <- sum(j <= i)
  }
  return(p)
}

#' @param p A vector of sparse matrix pointers
#'
#' @return \code{PointerToIndex}: A vector of column (j) indices
#'
#' @rdname SparsePointers
#'
#' @keywords internal
#'
#' @source \code{PointerToIndex} came from
#' \href{https://stackoverflow.com/questions/20008200/r-constructing-sparse-matrix}{StackOverflow}
#' @author \code{PointerToIndex} was written by
#' \href{https://stackoverflow.com/users/980833/josh-obrien}{Josh O'Brien on StackOverflow}
#'
PointerToIndex <- function(p) {
  dp <- diff(x = p)
  j <- rep.int(x = seq_along(along.with = dp), times = dp)
  return(j)
}
