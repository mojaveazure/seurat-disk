#' @importFrom methods setGeneric
#'
NULL

#' Test an object's class
#'
#' Generic functions for testing if an object satisfies class membership
#'
#' @name TestObject
#' @rdname TestObject
#'
#' @keywords internal
#'
NULL

#' @rdname TestObject
#'
setGeneric(
  name = 'IsDataFrame',
  def = function(x) {
    standardGeneric(f = 'IsDataFrame')
  },
  signature = c('x')
)

#' @rdname TestObject
#'
setGeneric(
  name = 'IsFactor',
  def = function(x) {
    standardGeneric(f = 'IsFactor')
  },
  signature = c('x')
)

#' @rdname TestObject
#'
setGeneric(
  name = 'IsList',
  def = function(x) {
    standardGeneric(f = 'IsList')
  },
  signature = c('x')
)

#' @rdname TestObject
#'
setGeneric(
  name = 'IsMatrix',
  def = function(x) {
    standardGeneric(f = 'IsMatrix')
  },
  signature = c('x')
)
