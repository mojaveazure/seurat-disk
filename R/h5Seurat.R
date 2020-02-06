#' @include scdisk.R
#'
NULL

h5seurat.validate <- list(
  '3.1.2' = function(hfile, verbose = TRUE) {
    valid <- TRUE
    .NotYetImplemented()
    return(valid)
  }
)

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
    # Methods
    validate = function(verbose = TRUE, ...) {
      message("Validating h5Seurat file")
      if (self$mode %in% modes$new) {
        ''
      }
      return(invisible(x = NULL))
    }
  )
)

#' @rdname as.h5Seurat
#' @method as.h5Seurat H5File
#' @export
#'
as.h5Seurat.H5File <- function(x, ...) {
  return(h5Seurat$new(filename = x$filename, mode = x$mode))
}

#' @param filename Name of file to save \code{x} to
#' @param overwrite Overwrite \code{filename} if present?
#' @param verbose Show progress updates
#'
#' @importFrom tools file_ext
#' @importFrom Seurat Project
#'
#' @rdname as.h5Seurat
#' @method as.h5Seurat Seurat
#' @export
#'
as.h5Seurat.Seurat <- function(
  x,
  filename = paste0(Project(object = x), '.h5seurat'),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  if (!grepl(pattern = '\\.h5seurat', x = file_ext(x = filename), ignore.case = TRUE)) {
    filename <- paste0(filename, '.h5seurat')
  }
  if (file.exists(filename)) {
    if (overwrite) {
      warning(
        "Overwriting previous file ",
        filename,
        call. = FALSE,
        immediate. = TRUE
      )
      file.remove(filename)
    } else {
      stop("H5Seurat file at ", filename, " already exists", call. = FALSE)
    }
  }
  hfile <- h5Seurat$new(filename = filename, mode = 'w')
  return(hfile)
}

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

#' Save an object as an h5seurat file
#'
#' @inheritParams as.h5Seurat
#'
#' @return Invisibly returns h5seurat file name
#'
#' @export
#'
SaveH5Seurat <- function(
  object,
  filename = 'object.h5seurat',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  h5seurat <- as.h5Seurat(
    x = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose
  )
  h5seurat$close_all()
  return(invisible(x = NULL))
}
