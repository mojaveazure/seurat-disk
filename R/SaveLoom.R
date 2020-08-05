#' @importFrom Seurat Project
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Save a \code{\link[Seurat]{Seurat}} object to a loom file
#'
#' @inheritParams SaveH5Seurat
#'
#' @return \code{SaveLoom}: Invisibly returns \code{filename}
#'
#' @name SaveLoom
#' @rdname SaveLoom
#'
#' @export
#'
SaveLoom <- function(object, filename, overwrite = FALSE, verbose = TRUE, ...) {
  UseMethod(generic = 'SaveLoom', object = object)
}

#' @return \code{as.loom}: A \code{\link{loom}} object
#'
#' @name SaveLoom
#' @rdname SaveLoom
#'
#' @export
#'
as.loom <- function(x, ...) {
  UseMethod(generic = 'as.loom', object = x)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname SaveLoom
#' @method SaveLoom default
#' @export
#'
SaveLoom.default <- function(
  object,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  tryCatch(
    expr = object <- as.Seurat(object = object, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        paste(class(x = object), collapse = ', '),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = object), '.loom')
  }
  return(invisible(x = SaveLoom(
    object = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )))
}

#' @rdname SaveLoom
#' @method SaveLoom Seurat
#' @export
#'
SaveLoom.Seurat <- function(
  object,
  filename = paste0(Project(object = object), '.loom'),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  loom <- as.loom(
    x = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  loom$close_all()
  return(invisible(x = loom$filename))
}

#' @rdname SaveLoom
#' @method as.loom default
#' @export
#'
as.loom.default <- function(x, filename, overwrite = FALSE, verbose = TRUE) {
  tryCatch(
    expr = x <- as.Seurat(object = x, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        paste(class(x = x), collapse = ', '),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = x), '.loom')
  }
  return(as.loom(
    object = x,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname SaveLoom
#' @method as.loom H5File
#' @export
#'
as.loom.H5File <- function(x, ...) {
  .NotYetImplemented()
}

#' @rdname SaveLoom
#' @method as.loom Seurat
#' @export
#'
as.loom.Seurat <- function(
  x,
  filename = paste0(Project(object = x), '.loom'),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  .NotYetImplemented()
}
