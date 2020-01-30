#' @include scdisk.R
#'
NULL

#' A class for connections to loom files
#'
#' @docType class
#' @name loom-class
#' @rdname loom-class
#' @aliases loom
#' @format An \code{\link[R6]{R6Class}} object
#' @seealso \code{\link[hdf5r]{H5File}}
#'
#' @importFrom R6 R6Class
#' @importFrom hdf5r H5File
#'
#' @export
#'
loom <- R6Class(
  classname = 'loom',
  inherit = scdisk,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
  ),
  private = list(
    # Methods
    validate = function(validate = TRUE) {
      invisible(x = NULL)
    }
  )
)
