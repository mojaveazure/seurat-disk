#' @include zzz.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load data from a loom file into a \code{\link[Seurat]{Seurat}} object
#'
#' @inheritParams LoadH5Seurat
#'
#' @return A \code{\link[Seurat]{Seurat}} object
#'
#' @name LoadLoom
#' @rdname LoadLoom
#'
#' @export
#'
LoadLoom <- function(file, ...) {
  UseMethod(generic = 'LoadLoom', object = file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LoadLoom
#' @method LoadLoom character
#' @export
#'
LoadLoom.character <- function(file, ...) {
  loom <- loom$new(filename = file, mode = 'r')
  on.exit(expr = loom$close_all())
  return(LoadLoom(file = loom, ...))
}

#' @rdname LoadLoom
#' @method LoadLoom H5File
#' @export
#'
LoadLoom.H5File <- function(file, ...) {
  return(LoadLoom(
    file = as.loom(x = file),
    ...
  ))
}

#' @rdname LoadLoom
#' @method LoadLoom loom
#' @export
#'
LoadLoom.loom <- function(file, ...) {
  return(as.Seurat(x = file, ...))
}

#' @aliases as.Seurat
#'
#' @rdname LoadLoom
#' @method as.Seurat loom
#' @export
#'
as.Seurat.loom <- function(x, ...) {
  .NotYetImplemented()
}
