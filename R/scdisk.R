#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' A disk-based object for single-cell analysis
#'
#' @docType class
#' @name scdisk-class
#' @rdname scdisk-class
#' @aliases scdisk
#' @format An \code{\link[R6]{R6Class}} object
#' @seealso \code{\link[hdf5r]{H5File}}
#'
#' @importFrom R6 R6Class
#' @importFrom hdf5r H5File
#'
#' @export
#'
scdisk <- R6Class(
  classname = 'scdisk',
  inherit = H5File,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
    #' @description Create a new \code{scdisk} object
    #' @param filename Name of on-disk file to connect to
    #' @param mode How to open the file, choose from:
    #' \describe{
    #'  \item{a}{Create new or open existing file, allow read and write}
    #'  \item{r}{Open existing file, allow read only}
    #'  \item{r+}{Open existing file, allow read and write}
    #'  \item{w}{Create new file (deleting any existing one), allow read and write}
    #'  \item{w-, x}{Create new file (error if exists), allow read and write}
    #' }
    #' @param validate Validate the file upon connection
    #' @param ... Extra arguments passed to validation routine
    initialize = function(
      filename = NULL,
      mode = c('a', 'r', 'r+', 'w', 'w-', 'x'),
      validate = TRUE,
      ...
      ) {
      mode <- match.arg(arg = mode)
      if (!file.exists(filename) && !mode %in% c('a', 'w', 'w-', 'x')) {
        stop("Cannot find file ", filename, call. = FALSE)
        }
      super$initialize(filename = filename, mode = mode)
      private$validate(validate = validate, ...)
    },
    #' @description Handle the loss of reference to this \code{scdisk} object
    finalizer = function() {
      self$close_all(close_self = TRUE)
    },
    #' @description Generate chunk points for a dataset
    #' @param dataset Name of dataset
    #' @param MARGIN Direction to chunk in; defaults to largest dimension of dataset
    #' @param csize Size of chunk; defaults to hdf5r-suggested chunk size
    #' @return A matrix where each row is a chunk, column 1 is start points,
    #' column 2 is end points
    chunk.points = function(
      dataset,
      MARGIN = getOption(x = 'SeuratDisk.chunking.MARGIN', default = 'largest'),
      csize = NULL
    ) {
      if (!self$exists(name = dataset)) {
        stop(
          "Cannot find dataset ",
          dataset,
          " in this ",
          class(x = self),
          " file",
          call. = FALSE)
      } else if (!inherits(x = self[[dataset]], what = 'H5D')) {
        stop("'dataset' must be an HDF5 dataset", call. = FALSE)
      }
      MARGIN <- GetMargin(dims = self[[dataset]]$dims, MARGIN = MARGIN)
      csize <- csize %||% self[[dataset]]$chunk_dims[MARGIN]
      return(ChunkPoints(dsize = self[[dataset]]$dims[MARGIN], csize = csize))
    }
  ),
  private = list(
    # Methods
    # @description Prebuilt error messages to reduce code duplication
    # @param type Type of error message to produce, choose from
    # \describe{
    #  \item{mode}{Cannot modify a file that has been opened as read-only}
    #  \item{ambiguous}{Unable to uniquely locate a dataset}
    # }
    # @return A character vector with the error message requested
    # @keywords internal
    errors = function(type = c('mode', 'ambiguous')) {
      type <- match.arg(arg = type)
      return(switch(
        EXPR = type,
        'mode' = paste(
          'Cannot modify a',
          class(x = self)[1],
          'file in read-only mode'
        ),
        'ambiguous' = 'Cannot identify the dataset provided, found too many like it; please be more specific'
      ))
    },
    # @description Check to see if a name is a dataset or group
    # @param name Name of object in HDF5 file
    # @param type Type object is, choose from
    # \describe{
    #  \item{H5D}{A HDF5 dataset}
    #  \item{H5Group}{An HDF5 group}
    # }
    # @return \code{TRUE} if \code{name} exists and is of type \code{type},
    # otherwise \code{FALSE}
    is.data = function(name, type = c('H5D', 'H5Group')) {
      type <- match.arg(arg = type)
      return(self$exists(name = name) && inherits(x = self[[name]], what = type))
    },
    # @description Validate ...
    # @param ... Ignored for \code{scdisk} method
    # @note The validation routine should be overwritten by subclasses of
    # \code{scdisk}; each validation routine must take at least one argument:
    # \code{validate} (a logical). Extra arguments are allowed.
    validate = function(...) {
      if (class(x = self)[1] == 'scdisk') {
        stop("Cannot create an scdisk object directly", call. = FALSE)
      }
      warning(
        "No validation method present for ",
        class(x = self)[1],
        " files",
        call. = FALSE,
        immediate. = TRUE
      )
    }
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# External Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Get and Register \code{\link{scdisk}} Subclasses
#'
#' Mechanisms for registration of \code{\link{scdisk}} subclass generators for
#' use in functions that rely on the class definition instead of an object.
#'
#' @inheritParams IsSCDisk
#'
#' @return \code{GetSCDisk}: if \code{r6class} is \code{NULL}, then a vector of
#' all registered \code{scdisk} subclasses; otherwise, a
#' \link[R6::R6Class]{generator} for the requested \code{scdisk} subclass
#'
#' @name RegisterSCDisk
#' @rdname RegisterSCDisk
#'
#' @details
#' While \code{scdisk}-subclassed objects (eg. \code{\link{h5Seurat}} objects)
#' follow traditional inheritance patterns (can be determined through
#' \code{\link{inherits}}), the class definitions and object generators do not.
#' These functions provide a simple mechanism for adding and getting the
#' defintions of \code{scdisk} subclasses for functions that utilize the object
#' generators or other aspects of the class definition (such as
#' \code{\link{Convert}})
#'
#' To register a subclass of \code{scdisk}, simply add a call to
#' \code{RegisterSCDisk} in your \link[base::ns-hooks]{load hook}
#'
#' \preformatted{
#' .onLoad <- function(libname, pkgname) {
#'   RegisterSCDisk(classgen)
#'   # Other code to be run on load
#' }
#' }
#'
#' @export
#'
#' @examples
#' GetSCDisk()
#' GetSCDisk("h5Seurat")
#'
GetSCDisk <- function(r6class = NULL) {
  if (is.null(x = r6class)) {
    return(names(x = scdisk.types))
  }
  classname <- grep(
    pattern = paste0('^', r6class, '$'),
    x = names(x = scdisk.types),
    ignore.case = TRUE,
    value = TRUE
  )
  if (!length(x = classname)) {
    stop("Unknown file type: ", r6class, call. = FALSE)
  }
  return(scdisk.types[[classname]])
}

#' @name RegisterSCDisk
#' @rdname RegisterSCDisk
#'
#' @return \code{RegisterSCDisk}: adds \code{r6class} to the internal subclass
#' registry and invisibly returns \code{NULL}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' RegisterSCDisk(h5Seurat)
#' }
#'
RegisterSCDisk <- function(r6class) {
  if (isTRUE(x = IsSCDisk(r6class = scdisk))) {
    r6pkg <- environmentName(env = r6class$parent_env)
    scpkg <- environmentName(env = scdisk$parent_env)
    # Ensure unique scdisk classes
    if (r6class$classname %in% names(x = scdisk.types) && r6pkg != scpkg) {
      stop(
        r6class$classname, " already registered by package ",
        r6pkg,
        call. = FALSE
      )
    }
    scdisk.types[[r6class$classname]] <- r6class
  } else {
    warning(
      r6class$classname,
      " does not inherit from scdisk",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  return(invisible(x = NULL))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Does an R6 class inherit from scdisk
#'
#' @param r6class An \link[R6::R6Class]{R6 class generator} or a character name
#' of an R6 class generator
#'
#' @return If \code{r6class} inherits from scdisk, returns \code{TRUE};
#' otherwise, returns \code{FALSE}
#'
#' @importFrom R6 is.R6Class
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::IsSCDisk("H5File")
#' SeuratDisk:::IsSCDisk("scdisk")
#' SeuratDisk:::IsSCDisk("h5Seurat")
#' }
#'
IsSCDisk <- function(r6class) {
  if (is.character(x = r6class)) {
    r6class <- eval(expr = as.symbol(x = r6class))
  }
  if (is.R6Class(x = r6class)) {
    if (identical(x = r6class, y = scdisk)) {
      return(TRUE)
    }
    return(IsSCDisk(r6class = r6class$get_inherit()))
  }
  return(FALSE)
}
