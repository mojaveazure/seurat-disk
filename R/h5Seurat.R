#' @include scdisk.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definitions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' A class for connections to h5Seurat files
#'
#' @docType class
#' @name h5Seurat-class
#' @rdname h5Seurat-class
#' @aliases h5Seurat
#' @format An \code{\link[R6]{R6Class}} object
#' @seealso \code{\link[hdf5r]{H5File}}
#'
#' @importFrom R6 R6Class
#' @importFrom hdf5r H5File
#' @importFrom utils packageVersion
#'
#' @export
#'
h5Seurat <- R6Class(
  classname = 'h5Seurat',
  inherit = scdisk,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
  ),
  private = list(
    versions = c('3.1.2'),
    # Methods
    create = function(version, verbose = TRUE) {
      if (self$mode == 'r') {
        stop(private$errors(type = 'mode'), call. = FALSE)
      }
      version <- ClosestVersion(query = version, targets = private$versions)
      switch(
        EXPR = version,
        '3.1.2' = {
          if (verbose) {
            message("Creating h5Seurat file for version ", version)
          }
          for (group in c('assays', 'commands', 'graphs', 'misc', 'reductions', 'tools')) {
            self$create_group(name = group)
          }
          attrs <- c(
            'active.assay' = '',
            'project' = 'SeuratDiskProject',
            'version' = '3.1.2'
          )
          for (i in seq_along(along.with = attrs)) {
            self$create_attr(attr_name = names(x = attrs)[i], robj = attrs[i])
          }
        },
        stop("Unknown version ", version, call. = FALSE)
      )
      return(invisible(x = self))
    },
    validate = function(verbose = TRUE, ...) {
      if (self$mode %in% modes$new) {
        private$create(
          version = packageVersion(pkg = 'Seurat'),
          verbose = verbose
        )
      }
      if (verbose) {
        message("Validating h5Seurat file")
      }
      if (!self$attr_exists(attr_name = 'version')) {
        stop("Invalid h5Seurat file: cannot find attribute 'version'", call. = FALSE)
      }
      version <- self$attr_open(attr_name = 'version')$read()
      switch(
        EXPR = ClosestVersion(query = version, targets = private$versions),
        '3.1.2' = {
          ''
        },
        stop("Unknown version ", version, call. = FALSE)
      )
      return(invisible(x = NULL))
    }
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods for Seurat-defined generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat character
#' @export
#'
LoadH5Seurat.character <- function(file, ...) {
  hfile <- h5Seurat$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(LoadH5Seurat(file = hfile, ...))
}

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat H5File
#' @export
#'
LoadH5Seurat.H5File <- function(file, ...) {
  return(LoadH5Seurat(file = as.h5Seurat(x = file), ...))
}

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat h5Seurat
#' @export
#'
LoadH5Seurat.h5Seurat <- function(file, ...) {
  .NotYetImplemented()
}
