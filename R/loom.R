#' @include scdisk.R
#'
NULL

loom.validate <- list(
  '0.1.0' = function(lfile, verbose = TRUE) {
    .NotYetImplemented()
  },
  '3.0.0' = function(lfile, verbose = TRUE) {
    .NotYetImplemented()
  }
)

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
