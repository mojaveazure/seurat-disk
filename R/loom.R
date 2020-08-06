#' @include scdisk.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    #' @description Get version information
    #' @return A \code{\link[base]{numeric_version}} object with the loom
    #' specification version information
    version = function() {
      loom.version <- 'LOOM_SPEC_VERSION'
      if (AttrExists(x = self, loom.version)) {
        version <- h5attr(x = self, which = loom.version)
      } else if (Exists(x = self, name = H5Path('attrs', loom.version))) {
        version <- self[[H5Path('attrs', loom.version)]][]
      } else {
        stop("Cannot find version information in this loom file", call. = FALSE)
      }
      return(numeric_version(x = version))
    },
    #' @description Add a timestamp to a dataset or group as an HDF5 attribute
    #' @param name Name of dataset or group to add timestamp to; if \code{NULL},
    #' timestamps the file as a whole
    #' @return Invisilby returns the object
    timestamp = function(name = NULL) {
      super$timestamp(
        name = name,
        attr = 'last_modified',
        tz = 'UTC',
        format = TSFormats(type = 'loom')
      )
      return(invisible(x = self))
    },
    #' @description Retrieve a timestamp from a dataset or group
    #' @param name Name of dataset or group to retrieve timestamp from; if
    #' \code{NULL}, retrieves timestamp from at the file-level
    #' @param locale Change the timestamp of to the timezone of the locale
    #' @return A character with the timestamp
    last.modified = function(name = NULL, locale = FALSE) {
      ds <- if (is.null(x = name)) {
        self
      } else {
        if (!Exists(x = self, name = name)) {
          stop("Cannot find a member named '", name, "'", call. = FALSE)
        }
        self[[name]]
      }
      if (!AttrExists(x = ds, name = 'last_modified')) {
        warning("No timestamp found", call. = FALSE, immediate. = TRUE)
        return(NA_character_)
      }
      time <- h5attr(x = ds, which = 'last_modified')
      time <- unlist(x = strsplit(x = time, split = '.', fixed = TRUE))[1]
      if (!grepl(pattern = 'Z$', x = time)) {
        time <- paste0(time, 'Z')
      }
      return(FormatTime(time = time, locale = locale))
    }
  ),
  private = list(
    # Methods
    validate = function(validate = TRUE) {
      invisible(x = NULL)
    }
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Check that a dataset is a proper loom matrix
#'
#' @param lfile A \code{\link{loom}} or \code{\link[hdf5r]{H5File}} object
#' @param name Name of matrix to check
#' @param dims If provided, ensure \code{lfile[[name]]} has these dimensions;
#' should be a two-dimensional numberic vector with ncells/ncol as the first
#' value and nfeature/nrow as the second
#'
#' @return If all checks pass succesfully, invisibly returns \code{name}
#'
#' @keywords internal
#'
CheckMatrix <- function(lfile, name, dims = NULL) {
  if (!inherits(x = lfile, what = "H5File")) {
    stop("Not an HDF5 file", call. = FALSE)
  }
  if (!Exists(x = lfile, name = name)) {
    stop("Cannot find '", name, "' in this file", call. = FALSE)
  } else if (!inherits(x = lfile[[name]], what = 'H5D') || !IsMatrix(x = lfile[[name]])) {
    stop("'", name, "' is not a matrix dataset", call. = FALSE)
  }
  if (is.numeric(x = dims) && length(x = dims) >= 2) {
    if (!all(Dims(x = lfile[[name]]) == dims[1:2])) {
      stop("'", name, "' does not match the provided dimensions", call. = FALSE)
    }
  }
  return(invisible(x = name))
}

#' Validate Loom Files
#'
#' Functions for validating loom files before connection
#'
#' @param ... ...
#'
#' @return ...
#'
#' @name ValiateLoom
#' @rdname ValidateLoom
#'
#' @keywords internal
#'
#' @seealso
#' \href{http://linnarssonlab.org/loompy/format/index.html}{Loom file format}
#'
LoomValidate0.1 <- function(lfile, verbose = TRUE) {
  .NotYetImplemented()
}

#' @name ValidateLoom
#' @rdname ValidateLoom
#'
LoomValidate3.0.0 <- function(lfile, verbose = TRUE) {
  # Check /matrix
  CheckMatrix(lfile = lfile, name = 'matrix')
  matrix.shape <- Dims(x = lfile[['matrix']])
  # Check /layers
  if (lfile$exists(name = 'layers')) {
    if (!inherits(x = lfile[['layers']], what = 'H5Group')) {
      stop("'layers' is not a group", call. = FALSE)
    }
    for (layer in names(x = lfile[['layers']])) {
      CheckMatrix(
        lfile = lfile,
        name = H5Path('layers', layer),
        dims = matrix.shape
      )
    }
  }
  # Check /attrs/LOOM_SPEC_VERSION
  if (!Exists(x = lfile, name = 'attrs/LOOM_SPEC_VERSION')) {
    stop("Cannot find the loom version", call. = FALSE)
  } else if (FALSE) {
    ''
  }
  # Check /row_attrs
  # Check /col_attrs
  # Check /row_graphs
  # Check /col_graphs
  return(invisible(x = TRUE))
}
