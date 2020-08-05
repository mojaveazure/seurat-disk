#' @importFrom rlang %||%
#' @importFrom hdf5r H5T_COMPOUND
#' @importFrom methods setOldClass
#'
NULL

#' @docType package
#' @name SeuratDisk-package
#' @rdname SeuratDisk-package
#'
#' @section Package options:
#'
#' SeuratDisk uses the following options to control behaviour, users can configure
#' these with \code{\link[base]{options}}:
#'
#' \describe{
#'  \item{\code{SeuratDisk.dtypes.logical_to_int}}{
#'   When writing \link[base]{logical} vectors, coerce to integer types to
#'   ensure compatibility across languages (see \code{\link{BoolToInt}} for
#'   more details)
#'  }
#'  \item{\code{SeuratDisk.dtypes.dataframe_as_group}}{
#'   When writing \link[base]{data.frame}s, always write out as a group
#'   regardless of factor presence
#'  }
#'  \item{\code{SeuratDisk.chunking.MARGIN}}{
#'   Default direction for chunking datasets; choose from:
#'   \describe{
#'    \item{largest}{Chunk along the largest dimension of a dataset}
#'    \item{smallest}{Chunk along the smallest dimension}
#'    \item{first}{Chunk along the first dimension}
#'    \item{last}{Chunk along the last dimension}
#'   }
#'  }
#' }
#'
#' @aliases SeuratDisk
#'
"_PACKAGE"

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Options
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

default.options <- list(
  "SeuratDisk.dtypes.logical_to_int" = TRUE,
  "SeuratDisk.dtypes.dataframe_as_group" = TRUE,
  "SeuratDisk.chunking.MARGIN" = c("largest", "smallest", "first", "last")
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global constants
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

modes <- list(
  'new' = c('w', 'w-', 'x'),
  'existing' = c('r', 'r+')
)

version.regex <- '^\\d+(\\.\\d+){2}(\\.9\\d{3})?$'

scdisk.types <- new.env()

spatial.version <- '3.1.5.9900'

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal utility functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert a logical to an integer
#'
#' Unlike most programming languages, R has three possible \link[base]{logical}
#' (boolean) values: \code{TRUE}, \code{FALSE}, and \code{\link[base]{NA}};
#' moreover, the \code{NA} value has representations in other data types, such
#' as \code{NA_integer_}, \code{NA_real_}, and \code{NA_character_}. Simply
#' writing out the logical values to an HDF5 file would cause issues when trying
#' to read the data in to another language, such as Python. To encode these three
#' logical values for other languages, we can encode the logicals as integers:
#' \itemize{
#'  \item \code{FALSE} becomes \code{0L}
#'  \item \code{TRUE} becomes \code{1L}
#'  \item \code{NA} becomes \code{2L}
#' }
#' This encoding scheme allows other languages to handle \code{NA}s in their own
#' manner while preserving all three logicals for R
#'
#' @param x A logical vector
#'
#' @return An integer vector
#'
#' @seealso \link[base]{integer} \link[base]{logical} \code{\link[base]{NA}}
#' \code{\link{WriteH5Seurat}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::BoolToInt(x = c(TRUE, FALSE, NA))
#' }
#'
BoolToInt <- function(x) {
  x <- as.integer(x = x)
  x[which(x = is.na(x = x))] <- 2L
  return(x)
}

#' Generate chunk points
#'
#' @param dsize Size of data being chunked
#' @param csize Size of chunk; if \code{NA}, assumes single chunk
#'
#' @return A matrix where each row is a chunk, column 1 is start points, column
#' 2 is end points
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::ChunkPoints(100, 3)
#' SeuratDisk:::ChunkPoints(100, NA)
#' }
#'
ChunkPoints <- function(dsize, csize) {
  if (is.na(x = csize)) {
    return(matrix(
      data = c(1, dsize),
      ncol = 2,
      dimnames = list(NULL, c('start', 'end'))
    ))
  }
  return(t(x = vapply(
    X = seq.default(from = 1L, to = ceiling(dsize / csize)),
    FUN = function(i) {
      return(c(
        start = (csize * (i - 1L)) + 1L,
        end = min(csize * i, dsize)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  )))
}

#' Find the closest version
#'
#' API changes happen at set versions, and knowing how a current running version
#' relates to versions introducing API changes is important.
#' \code{ClosestVersion} approximages both \dQuote{rounding down} (eg. to
#' determine minimum version with new API addition) and \dQuote{rounding up}
#' (eg. to determine maximum version before API deletion) for semantic versions.
#'
#' @param query A query version (\code{\link[base]{character}} or
#' \code{\link[base]{numeric_version}})
#' @param targets A vector of target versions (\code{\link[base]{character}} or
#' \code{\link[base]{numeric_version}})
#' @param direction Which way should we check for closest version? Choose from:
#' \describe{
#'  \item{min}{Closest version less than or equal to \code{query}}
#'  \item{max}{Closest version greater than or equal to \code{query}}
#' }
#' @param inclusive Perform an inclusive comparison (eg. \code{>=} or \code{<=}
#' versus to \code{>} or \code{<}) for \dQuote{rounding}
#'
#' @return The version from \code{targets} that is closest to \code{query} as a
#' \code{\link[base]{character}} vector
#'
#' @seealso \code{\link[base]{numeric_version}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::ClosestVersion('3.1.0', targets = c('3.0.0', '1.4.9', '4.3.2'))
#' SeuratDisk:::ClosestVersion('3.1.0', targets = c('3.0.0', '1.4.9', '4.3.2'), direction = 'max')
#' }
#'
ClosestVersion <- function(
  query,
  targets,
  direction = c('min', 'max'),
  inclusive = direction == 'min'
) {
  direction <- match.arg(arg = direction)
  query <- numeric_version(x = query)
  targets <- sort(x = numeric_version(x = targets))
  switch(
    EXPR = direction,
    'min' = {
      compare <- ifelse(test = inclusive, yes = `<=`, no = `<`)
      collapse <- max
    },
    'max' = {
      compare <- ifelse(test = inclusive, yes = `>=`, no = `>`)
      collapse <- min
    }
  )
  index <- suppressWarnings(expr = collapse(which(x = compare(
    e1 = targets,
    e2 = query
  ))))
  if (is.infinite(x = index)) {
    stop(
      "All target versions ",
      switch(EXPR = direction, 'min' = 'greater', 'max' = 'less'),
      " than query version (",
      as.character(x = query),
      ")",
      call. = FALSE
    )
  }
  return(as.character(x = targets[index]))
}

#' Convert an HDF5 compound dataset to a group
#'
#' @param src An HDF5 dataset (\code{\link[hdf5r]{H5D}}) of type
#' \code{\link[hdf5r]{H5T_COMPOUND}}
#' @param dest An HDF5 file (\code{\link[hdf5r]{H5File}}) or group
#' (\code{\link[hdf5r]{H5Group}})
#' @param dst.name Name of group in \code{dest}
#' @param order Name of HDF5 attribute to store column order as
#' @param index Integer values of which values to pull; defaults to all values
#' @param overwrite Overwrite existing group \code{dst.name} in \code{dest}
#'
#' @return Invisibly returns \code{NULL}
#'
#'
#' @keywords internal
#'
CompoundToGroup <- function(
  src,
  dest,
  dst.name = basename(path = src$get_obj_name()),
  order = c('colnames', 'column-order'),
  index = NULL,
  overwrite = FALSE
) {
  order <- match.arg(arg = order)
  if (!IsDType(src, 'H5T_COMPOUND')) {
    stop("'src' must be an HDF5 compound dataset", call. = FALSE)
  } else if (!inherits(x = dest, what = c('H5File', 'H5Group'))) {
    stop("'dest' must be a HDF5 file or group", call. = FALSE)
  }
  if (dest$exists(name = dst.name)) {
    if (overwrite) {
      dest$link_delete(name = dst.name)
    } else {
      stop(dst.name, " already exists in the destination", call. = FALSE)
    }
  }
  index <- index %||% seq.default(from = 1, to = src$dims)
  group <- dest$create_group(name = dst.name)
  cpd <- src$get_type()
  for (i in seq_along(along.with = cpd$get_cpd_labels())) {
    name <- cpd$get_cpd_labels()[i]
    dtype <- cpd$get_cpd_types()[[i]]
    group$create_dataset(
      name = name,
      robj = unlist(
        x = src$read_low_level(mem_type = H5T_COMPOUND$new(
          labels = name,
          dtypes = dtype
        )),
        use.names = FALSE
      )[index],
      dtype = dtype
    )
  }
  group$create_attr(
    attr_name = order,
    robj = cpd$get_cpd_labels(),
    dtype = GuessDType(x = cpd$get_cpd_labels())
  )
  return(invisible(x = NULL))
}

#' Determine a filetype based on its extension
#'
#' @param file Name of file
#'
#' @return The extension, all lowercase
#'
#' @importFrom tools file_ext
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::FileType('pbmc3k.h5Seurat')
#' SeuratDisk:::FileType('h5ad')
#' }
#'
FileType <- function(file) {
  ext <- file_ext(x = file)
  ext <- ifelse(test = nchar(x = ext), yes = ext, no = basename(path = file))
  return(tolower(x = ext))
}

#' Get a class string with package information
#'
#' S4 classes are useful in the context of their defining package (benefits of
#' stricter typing). In order to ensure class information is properly retained
#' in HDF5 files, S4 class names are written as \dQuote{package:classname} with
#' certain exceptions (eg. S4 classes defined by
#' \link[Seurat:Seurat-package]{Seurat})
#'
#' @param class Class name
#' @param packages A vector of packages to exclude from resulting class
#' information
#'
#' @return A character vector with the class
#'
#' @importFrom methods getClass slot
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::GetClass('Seurat')
#' SeuratDisk:::GetClass('Matrix')
#' }
#'
GetClass <- function(class, packages = 'Seurat') {
  class <- class[1]
  classdef <- getClass(Class = class)
  classpkg <- slot(object = classdef, name = 'package')
  if (classpkg %in% packages) {
    classpkg <- NULL
  }
  class <- paste(classpkg, class, sep = ':')
  return(gsub(pattern = '^:', replacement = '', x = class))
}

#' Determine the margin to use for a dataset
#'
#' @param dims Dimensions of a dataset
#' @param MARGIN Either an integer value contained within
#' \code{1:length(x = dims)} or one of the possible values of
#' \code{\link[SeuratDisk]{SeuratDisk.chunking.MARGIN}}
#'
#' @return An integer value with the \code{MARGIN}
#'
#' @seealso \code{\link[SeuratDisk]{SeuratDisk.chunking.MARGIN}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::GetMargin(c(4, 10))
#' }
#'
GetMargin <- function(dims, MARGIN = getOption(x = 'SeuratDisk.chunking.MARGIN')) {
  if (isFALSE(x = is.numeric(x = MARGIN))) {
    MARGIN <- tryCatch(
      expr = match.arg(
        arg = MARGIN,
        choices = default.options[['SeuratDisk.chunking.MARGIN']]
      ),
      error = function(err) {
        warning(err$message, call. = FALSE, immediate. = TRUE)
        return(default.options[['SeuratDisk.chunking.MARGIN']][1])
      }
    )
    MARGIN <- switch(
      EXPR = MARGIN,
      'largest' = which.max(x = dims),
      'smallest' = which.min(x = dims),
      'first' = 1L,
      'last' = length(x = dims)
    )
  }
  if (isFALSE(x = MARGIN %in% seq.int(from = 1, to = length(x = dims)))) {
    stop("'MARGIN' must be within the dimensions of the dataset", call. = FALSE)
  }
  return(MARGIN)
}

#' Get the parent of an HDF5 dataset or group
#'
#' @param x An HDF5 dataset or group
#'
#' @return An \code{\link[hdf5r]{H5File}} or \code{\link[hdf5r]{H5Group}} object
#'
#' @keywords internal
#'
GetParent <- function(x) {
  dname <- dirname(path = x$get_obj_name())
  dest <- if (dname == '/') {
    x$get_file_id()
  } else {
    x$get_file_id()[[dname]]
  }
  return(dest)
}

#' Guess an HDF5 Datatype
#'
#' Wrapper around \code{\link[hdf5r:guess_dtype]{hdf5r::guess_dtype}}, allowing
#' for the customization of string types rather than defaulting to
#' variable-length ASCII-encoded strings. Also encodes logicals as
#' \code{\link[hdf5r]{H5T_INTEGER}} instead of \code{\link[hdf5r]{H5T_LOGICAL}}
#' to ensure cross-language compatibility (controlled via
#' \link[=SeuratDisk-package]{package options})
#'
#' @inheritParams StringType
#' @inheritParams hdf5r::guess_dtype
#' @inheritDotParams hdf5r::guess_dtype
#'
#' @return An object of class \code{\link[hdf5r]{H5T}}
#'
#' @importFrom hdf5r guess_dtype
#'
#' @seealso \code{\link[hdf5r]{guess_dtype}} \code{\link{BoolToInt}}
#' \code{\link{StringType}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' # Characters can either be variable-width UTF8-encoded or
#' # fixed-width ASCII-encoded
#' SeuratDisk:::GuessDType(x = 'hello')
#' SeuratDisk:::GuessDType(x = 'hello', stype = 'ascii7')
#'
#' # Data frames are a compound type; character columns follow the same rules
#' # as character vectors
#' df <- data.frame(x = c('g1', 'g2', 'g3'), y = 1, 2, 3, stringsAsFactors = FALSE)
#' SeuratDisk:::GuessDType(x = df)
#' SeuratDisk:::GuessDType(x = df, stype = 'ascii7')
#'
#' # Logicals are turned into integers to ensure compatibility with Python
#' # TRUE evaluates to 1, FALSE to 0, and NA to 2
#' SeuratDisk:::GuessDType(x = c(TRUE, FALSE, NA))
#' }
#'
GuessDType <- function(x, stype = 'utf8', ...) {
  dtype <- guess_dtype(x = x, ...)
  if (inherits(x = dtype, what = 'H5T_STRING')) {
    dtype <- StringType(stype = stype)
  } else if (inherits(x = dtype, what = 'H5T_COMPOUND')) {
    cpd.dtypes <- dtype$get_cpd_types()
    for (i in seq_along(along.with = cpd.dtypes)) {
      if (inherits(x = cpd.dtypes[[i]], what = 'H5T_STRING')) {
        cpd.dtypes[[i]] <- StringType(stype = stype)
      }
    }
    dtype <- H5T_COMPOUND$new(
      labels = dtype$get_cpd_labels(),
      dtypes = cpd.dtypes,
      size = dtype$get_size()
    )
  } else if (inherits(x = dtype, what = 'H5T_LOGICAL')) {
    if (getOption(x = "SeuratDisk.dtypes.logical_to_int", default = TRUE)) {
      dtype <- guess_dtype(x = BoolToInt(x = x), ...)
    }
  }
  return(dtype)
}

#' Check the datatype of an HDF5 dataset
#'
#' Effectively, an implementation of \code{\link[methods]{is}} for HDF5 datasets;
#' useful to ensure HDF5 validity for specific file structures
#'
#' @param x An HDF5 dataset (object of type \code{\link[hdf5r]{H5D}})
#' @param dtype A character vector of HDF5 datatype names, must be present in
#' \code{\link[hdf5r]{h5types}}
#'
#' @return A logical
#'
#' @importFrom hdf5r h5types
#'
#' @seealso \code{\link[hdf5r]{h5types}}
#'
#' @keywords internal
#'
IsDType <- function(x, dtype) {
  if (!inherits(x = x, what = 'H5D')) {
    stop("'IsDType' only works on HDF5 dataset", call. = FALSE)
  }
  dtypes <- unique(x = sapply(
    X = grep(pattern = '^H5T_', x = names(x = h5types), value = TRUE),
    FUN = function(i) {
      return(class(x = h5types[[i]])[1])
    },
    USE.NAMES = FALSE
  ))
  dtypes <- unique(x = c(dtypes, 'H5T_COMPOUND'))
  match.arg(arg = dtype, choices = dtypes, several.ok = TRUE)
  missing.dtypes <- setdiff(x = dtype, y = dtypes)
  if (length(x = missing.dtypes)) {
    dtype <- setdiff(x = dtype, y = missing.dtypes)
    if (!length(x = dtype)) {
      stop("None of the requested dtypes are valid HDF5 datatypes", call. = FALSE)
    } else {
      warning(
        "The following requested dtypes are not valid HDF5 datatypes: ",
        paste(missing.dtypes, sep = ", "),
        call. = FALSE,
        immediate. = TRUE
      )
    }
  }
  return(inherits(x = x$get_type(), what = dtype))
}

#' Check to see if a matrix is empty
#'
#' Determine if a matrix is empty or not. A matrix is considered empty if it
#' satisfies one of the following conditions:
#' \itemize{
#'  \item The dimensions of the matrix are 0-by-0 (\code{all(dim(x) == 0)})
#'  \item The dimensions of the matrix are 1-by-1 (\code{all(dim(x) == 1)}) and
#'  the sole vlaue is \code{NA}
#' }
#' These two situations correspond to matrices generated with either
#' \code{new('matrix')} or \code{matrix()}
#'
#' @param x A matrix
#'
#' @return \code{TRUE} if the matrix is empty otherwise \code{FALSE}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::IsMatrixEmpty(new('matrix'))
#' SeuratDisk:::IsMatrixEmpty(matrix())
#' SeuratDisk:::IsMatrixEmpty(matrix(1:9, nrow = 3))
#' }
#'
IsMatrixEmpty <- function(x) {
  matrix.dims <- dim(x = x)
  matrix.na <- all(matrix.dims == 1) && all(is.na(x = x))
  return(all(matrix.dims == 0) || matrix.na)
}

#' Make a space
#'
#' Generate a blank space \code{n} characters long; useful for aligning text to
#' be printed to console
#'
#' @param n Length space should be
#'
#' @return A space (' ') of length \code{n}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::MakeSpace(n = 10)
#' cat('hello', SeuratDisk:::MakeSpace(n = 10), 'world\n', sep = '')
#' }
#'
MakeSpace <- function(n) {
  return(paste(rep_len(x = ' ', length.out = n), collapse = ''))
}

#' Add names for unnamed or partially named objects
#'
#' @param x An object that can be named
#' @param prefix A prefix to be added to each name
#'
#' @return \code{x} with unnamed values named
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' a <- list(1, b = 2, 3)
#' SeuratDisk:::PadNames(a)
#' }
#'
PadNames <- function(x, prefix = 'index') {
  if (length(x = x)) {
    xnames <- names(x = x) %||% paste0(
      prefix,
      seq.int(from = 1L, to = length(x = x))
    )
    missing <- which(x = !nchar(x = xnames))
    if (length(x = missing)) {
      xnames[missing] <- paste0(prefix, missing)
    }
    names(x = x) <- xnames
  }
  return(x)
}

#' Create a progress bar
#'
#' Progress bars are useful ways of getting updates on how close a task is to
#' completion. However, they can get in the way of RMarkdown documents with
#' lots of unnecesssary printing. \code{PB} is a convenience function that
#' creates progress bars with the following defaults
#' \itemize{
#'  \item \code{char = '='}
#'  \item \code{style = 3}
#'  \item \code{file = stderr()}
#' }
#'
#' @return An object of class \code{\link[utils]{txtProgressBar}}
#'
#' @importFrom utils txtProgressBar
#'
#' @seealso \code{\link[utils]{txtProgressBar}} \code{\link[base]{stderr}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' pb <- SeuratDisk:::PB()
#' for (i in 1:10) {
#'   utils::setTxtProgressBar(pb, i / 10)
#' }
#' close(pb)
#' }
PB <- function() {
  return(txtProgressBar(char = '=', style = 3, file = stderr()))
}

#' Generate a random string of characters
#'
#' @param length Length (\code{\link[base]{nchar}}) of string to generate
#' @param ... Extra parameters passed to \code{\link[base]{sample}}
#'
#' @return A random string of characters of length (\code{\link[base]{nchar}})
#' of \code{length}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::RandomName()
#' }
#'
RandomName <- function(length = 5L, ...) {
  # CheckDots(..., fxns = "sample")
  return(paste(sample(x = letters, size = length, ...), collapse = ""))
}

#' Create a scalar space
#'
#' @return An object of type \code{\link[hdf5r:H5S]{H5S}} denoting a scalar HDF5
#' space
#'
#' @keywords internal
#'
Scalar <- function() {
  return(H5S$new(type = 'scalar'))
}

#' Generate an HDF5 string dtype
#'
#' Presets for encoding variations of \code{\link[hdf5r]{H5T_STRING}}; used to
#' generate HDF5 datatype specifications with specific string encodings
#'
#' @param stype Type of string encoding to use, choose from:
#' \describe{
#'  \item{utf8}{Variable-width, UTF-8}
#'  \item{ascii7}{Fixed-width (7 bits), ASCII}
#' }
#'
#' @return An \code{\link[hdf5r]{H5T_STRING}} object
#'
#' @importFrom hdf5r h5const H5T_STRING
#'
#' @seealso \code{\link[hdf5r]{H5T_STRING}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::StringType()
#' SeuratDisk:::StringType('ascii7')
#' }
#'
StringType <- function(stype = c('utf8', 'ascii7')) {
  stype <- match.arg(arg = stype)
  return(switch(
    EXPR = stype,
    'utf8' = H5T_STRING$new(size = Inf)$set_cset(cset = h5const$H5T_CSET_UTF8),
    'ascii7' = H5T_STRING$new(size = 7L)
  ))
}

#' Update a Seurat key
#'
#' Attempts to validate a string to use as a Seurat key. Valid keys must match
#' the regular expression \code{^[[:alnum:]]+_$}; if \code{key} fails this
#' regular expression, an attempt to modify it to said key will be made by
#' removing all non-alphanumeric characters, collapsing the resulting vector,
#' and appending \dQuote{_}. If this stil fails, a random string of lowercase
#' characters will be generated, followed by \dQuote{_}, to be used as the key
#'
#' @param key A key to validate and update
#'
#' @return \code{key}, updated if invalid
#'
#' @seealso \code{\link[Seurat]{Key}} \code{\link{RandomName}}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::UpdateKey("RNA_")
#' SeuratDisk:::UpdateKey("potato")
#' SeuratDisk:::UpdateKey("*@)")
#' }
#'
UpdateKey <- function(key) {
  if (grepl(pattern = "^[[:alnum:]]+_$", x = key)) {
    return(key)
  } else {
    new.key <- regmatches(x = key, m = gregexpr(pattern = "[[:alnum:]]+",text = key))
    new.key <- paste0(paste(unlist(x = new.key), collapse = ""), "_")
    if (new.key == "_") {
      new.key <- paste0(RandomName(length = 3), "_")
    }
    return(new.key)
  }
}

#' Update slots in an object
#'
#' @param object An object to update
#'
#' @return \code{object} with the latest slot definitions
#'
#' @importFrom methods slotNames slot slot<-
#'
#' @keywords internal
#'
UpdateSlots <- function(object) {
  object.list <- sapply(
    X = slotNames(x = object),
    FUN = function(x) {
      return(tryCatch(
        expr = slot(object = object, name = x),
        error = function(...) {
          return(NULL)
        }
      ))
    },
    simplify = FALSE,
    USE.NAMES = TRUE
  )
  object.list <- Filter(f = Negate(f = is.null), x = object.list)
  object.list <- c('Class' = class(x = object)[1], object.list)
  object <- do.call(what = 'new', args = object.list)
  for (x in setdiff(x = slotNames(x = object), y = names(x = object.list))) {
    xobj <- slot(object = object, name = x)
    if (is.vector(x = xobj) && !is.list(x = xobj) && length(x = xobj) == 0) {
      slot(object = object, name = x) <- vector(
        mode = class(x = xobj),
        length = 1L
      )
    }
  }
  return(object)
}

#' Write an attribute to an HDF5 file, group, or dataset
#'
#' @param h5 An HDF5 \link[hdf5r:H5File]{file}, \link[hdf5r:H5Group]{group}, or
#' \link[hdf5r:H5D]{dataset}
#' @param name Name to store attribute as
#' @param robj An object to write out
#' @param dtype Data type of attribute
#' @param scalar Is this a scalar or simple (vectorized) attribute?
#' @param overwrite Overwrite the attribute if it already exists
#' @param ... Extra paramters passed to \code{\link[hdf5r:H5S]{H5S$new}}
#'
#' @return Invisibly returns \code{NULL}
#'
#' @importFrom hdf5r H5S
#'
#' @keywords internal
#'
WriteAttribute <- function(
  h5,
  name,
  robj,
  dtype = GuessDType(x = robj),
  scalar = length(x = robj) == 1,
  overwrite = FALSE,
  ...
) {
  if (!inherits(x = h5, what = c('H5File', 'H5Group', 'H5D'))) {
    stop("'h5' must be an HDF5 file, group, or dataset", call. = FALSE)
  }
  if (h5$attr_exists(attr_name = name)) {
    if (overwrite) {
      h5$attr_delete(attr_name = name)
    } else {
      stop("Attribute ", name, " already exists", call. = FALSE)
    }
  }
  if (is.logical(x = robj) && getOption(x = "SeuratDisk.dtypes.logical_to_int", default = TRUE)) {
    robj <- BoolToInt(x = robj)
  }
  space.type <- ifelse(test = isTRUE(x = scalar), yes = 'scalar', no = 'simple')
  dims <- if (space.type == 'scalar') {
    NULL
  } else {
    dim(x = robj) %||% length(x = robj)
  }
  h5$create_attr(
    attr_name = name,
    robj = robj,
    dtype = dtype,
    space = H5S$new(type = space.type, dims = dims, ...)
  )
  return(invisible(x = NULL))
}

#' Get the proper HDF5 connection mode for writing depending on overwrite status
#'
#' @param overwrite Overwrite a file
#'
#' @return \code{w} if \code{overwrite} else \code{w-}
#'
#' @keywords internal
#'
#' @examples
#' \donttest{
#' SeuratDisk:::WriteMode(TRUE)
#' SeuratDisk:::WriteMode(FALSE)
#' }
#'
WriteMode <- function(overwrite = FALSE) {
  return(ifelse(test = overwrite, yes = 'w', no = 'w-'))
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Loading handler
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

.onLoad <- function(libname, pkgname) {
  # Make the classes defined in SeuratDisk compatible with S4 generics/methods
  # setOldClass(Classes = c('scdisk', 'h5Seurat', 'loom'))
  setOldClass(Classes = c('scdisk', 'h5Seurat'))
  RegisterSCDisk(r6class = h5Seurat)
  # Set some default options
  op <- options()
  toset <- !names(x = default.options) %in% names(x = op)
  if (any(toset)) {
    options(default.options[toset])
  }
  invisible(x = NULL)
}
