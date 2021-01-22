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
    #' @description Add an attribute
    #' @param x Object to add as an attribute
    #' @param name Name to store attribute as
    #' @param type Type of attribute to add
    add_attribute = function(x, name, type = c('global', 'row', 'col')) {
      if (is.list(x = x)) {
        stop("Cannot add lists as attributes", call. = FALSE)
      }
      type <- match.arg(arg = type)
      name <- basename(path = name)
      if (type == 'global') {
        WriteAttribute(x = x, name = name, lfile = self, stype = private$stype())
      } else {
        SizedAttribute(
          x = x,
          name = name,
          lfile = self,
          type = type,
          stype = private$stype()
        )
      }
      return(invisible(x = self))
    },
    #' @description Add a graph
    #' @param x ...
    #' @param name ...
    #' @param type ...
    #' @param verbose ...
    add_graph = function(x, name, type = c('col', 'row'), verbose = TRUE) {
      type <- match.arg(arg = type)
      name <- basename(path = name)
      WriteGraph(
        x = x,
        name = name,
        lfile = self,
        type = type,
        verbose = verbose
      )
      return(invisible(x = self))
    },
    #' @description Add a layer to this loom file
    #' @param x An object to save as a layer
    #' @param name Name to store layer as
    #' @param transpose ...
    #' @param verbose ...
    #' @return Invisibly returns \code{NULL}
    add_layer = function(x, name, transpose = TRUE, verbose = TRUE) {
      dims <- if (isTRUE(x = transpose)) {
        dim(x = self)
      } else {
        rev(x = dim(x = self))
      }
      if (!identical(x = dim(x = x), y = dims)) {
        stop(
          "'x' must be a matrix-like object with ",
          dims[1],
          " rows  and ",
          dims[2],
          " columns",
          call. = FALSE
        )
      }
      layers <- c('layers', '/layers')
      if (!dirname(path = name) %in% layers) {
        name <- H5Path('/layers', name)
      }
      if (!dirname(path = name) %in% layers) {
        stop("Subgroups are not allowed in layers")
      }
      if (isTRUE(x = verbose)) {
        message("Adding layer ", basename(path = name))
      }
      WriteMatrix(
        x = x,
        name = name,
        lfile = self,
        transpose = transpose,
        verbose = verbose
      )
      return(invisible(x = self))
    },
    #' @description Get version information
    #' @return A \code{\link[base]{numeric_version}} object with the loom
    #' specification version information
    version = function() {
      loom.version <- 'LOOM_SPEC_VERSION'
      version <- if (AttrExists(x = self, loom.version)) {
        h5attr(x = self, which = loom.version)
      } else if (Exists(x = self, name = H5Path('attrs', loom.version))) {
        self[[H5Path('attrs', loom.version)]][]
      } else {
        version <- ifelse(
          test = Exists(x = self, name = '/attrs'),
          yes = '3.0.0',
          no = '0.1.0'
        )
        warning(
          "Cannot find version information in this loom file, assuming to be ",
          version,
          call. = FALSE,
          immediate. = TRUE
        )
        # NULL
        version
      }
      return(numeric_version(x = version))
    },
    #' @description Add a timestamp to a dataset or group as an HDF5 attribute
    #' @param name Name of dataset or group to add timestamp to; if \code{NULL},
    #' timestamps the file as a whole
    #' @return Invisibly returns the object
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
    create = function() {
      if (self$mode == 'r') {
        stop(private$errors(type = 'mode'), call. = FALSE)
      }
      groups <- c(
        'layers',
        'attrs',
        paste0(c('row', 'col'), '_attrs'),
        paste0(c('row', 'col'), '_graphs')
      )
      for (group in groups) {
        self$create_group(name = group)
      }
      suppressWarnings(expr = self$add_attribute(
        x = '3.0.0',
        name = 'LOOM_SPEC_VERSION'
      ))
      return(invisible(x = NULL))
    },
    stype = function() {
      version <- try(
        expr = ClosestVersion(
          query = self$version(),
          targets = c('0.1.0', '3.0.0')
        ),
        silent = TRUE
      )
      return(switch(
        EXPR = version,
        '0.1.0' = 'ascii7',
        '3.0.0' = 'utf8',
        {
          warning(
            "Unknown loom version ",
            version,
            ", using UTF-8 as string type",
            call. = FALSE,
            immediate. = TRUE
          )
          'utf8'
        }
      ))
    },
    validate = function(validate = TRUE) {
      if (self$mode %in% modes$new) {
        private$create()
      }
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
#' should be a two-dimensional numeric vector with ncells/ncol as the first
#' value and nfeature/nrow as the second
#'
#' @return If all checks pass successfully, invisibly returns \code{name}
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

methods::setGeneric(
  name = 'SizedAttribute',
  def = function(x, name, lfile, type, stype, verbose = TRUE) {
    if (!inherits(x = lfile, what = 'loom')) {
      stop("'lfile' must be a loom object", call. = FALSE)
    }
    standardGeneric(f = 'SizedAttribute')
  },
  signature = c('x')
)

methods::setMethod(
  f = 'SizedAttribute',
  signature = c('x' = 'vector'),
  definition = function(x, name, lfile, type, stype, verbose = TRUE) {
    if (is.list(x = x)) {
      stop("Cannot add lists as attributes", call. = FALSE)
    }
    dfun <- switch(EXPR = type, 'row' = nrow, 'col' = ncol)
    if (!length(x = x) == dfun(x = lfile)) {
      stop("Wrong attribute length", call. = FALSE)
    }
    if (isTRUE(x = verbose)) {
      message("Adding ", type, " attribute ", name)
    }
    lfile$create_dataset(
      name = H5Path(paste0(type, '_attrs'), name),
      robj = x,
      dtype = GuessDType(x = x, stype = stype)
    )
  }
)

methods::setGeneric(
  name = 'WriteAttribute',
  def = function(x, name, lfile, stype, transpose = TRUE, verbose = TRUE) {
    if (!inherits(x = lfile, what = 'loom')) {
      stop("'lfile' must be a loom object", call. = FALSE)
    }
    standardGeneric(f = 'WriteAttribute')
  },
  signature = c('x')
)

methods::setMethod(
  f = 'WriteAttribute',
  signature = c('x' = 'matrix'),
  definition = function(
    x,
    name,
    lfile,
    stype,
    transpose = TRUE,
    verbose = TRUE
  ) {
    if (isTRUE(x = transpose)) {
      x <- t(x = x)
    }
    lfile$create_dataset(
      name = H5Path('/attrs', name),
      robj = x,
      dtype = GuessDType(x = x[1, 1], stype = stype)
    )
    return(invisible(x = NULL))
  }
)

methods::setMethod(
  f = 'WriteAttribute',
  signature = c('x' = 'vector'),
  definition = function(
    x,
    name,
    lfile,
    stype,
    transpose = TRUE,
    verbose = TRUE
  ) {
    if (is.list(x = x)) {
      stop("Cannot add lists as attributes", call. = FALSE)
    }
    lfile$create_dataset(
      name = H5Path('/attrs', name),
      robj = x,
      dtype = GuessDType(x = x, stype = stype)
    )
    return(invisible(x = NULL))
  }
)

methods::setGeneric(
  name = 'WriteGraph',
  def = function(x, name, lfile, type, verbose = TRUE) {
    if (!inherits(x = lfile, what = 'loom')) {
      stop("'lfile' must be a loom object", call. = FALSE)
    }
    standardGeneric(f = 'WriteGraph')
  },
  signature = c('x')
)

#' @importClassesFrom Matrix dgCMatrix
#'
methods::setMethod(
  f = 'WriteGraph',
  signature = c('x' = 'dgCMatrix'),
  definition = function(x, name, lfile, type, verbose = TRUE) {
    graph <- lfile$create_group(name = H5Path(paste0(type, '_graphs'), name))
    graph$create_dataset(
      name = 'a',
      robj = methods::slot(object = x, name = 'i'),
      dtype = GuessDType(x = 1L)
    )
    graph$create_dataset(
      name = 'b',
      robj = PointerToIndex(p = methods::slot(object = x, name = 'p')) - 1,
      dtype = GuessDType(x = 1L)
    )
    graph$create_dataset(
      name = 'w',
      robj = methods::slot(object = x, name = 'x'),
      dtype = GuessDType(x = 1)
    )
    return(invisible(x = NULL))
  }
)
