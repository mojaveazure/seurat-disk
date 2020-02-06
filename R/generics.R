#' Save an object as an H5Seurat file
#'
#' @param x An object
#' @param ... Arguments passed to other methods
#'
#' @return An \code{\link{h5Seurat}} object
#'
#' @export
#'
as.h5Seurat <- function(x, ...) {
  UseMethod(generic = 'as.h5Seurat', object = x)
}

#' Save an object as a loom file
#'
#' @inheritParams as.h5Seurat
#'
#' @return A \code{\link{loom}} object
#'
#' @export
#'
as.loom <- function(
  x,
  filename,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  UseMethod(generic = 'as.loom', object = x)
}

#' Load a saved \code{Seurat} object from an h5Seurat file
#'
#' @param file Name of h5Seurat or connected h5Seurat file to load
#' @param assays One of:
#' \itemize{
#'   \item A character vector with names of assays
#'   \item A character vector with one or more of \code{counts}, \code{data},
#'   \code{scale.data} describing which slots of \strong{all assays} to load
#'   \item A named list where each entry is either the name of an assay or a vector
#'   describing which slots (described above) to take from which assay
#'   \item \code{NULL} for all assays
#' }
#' @param reductions One of:
#' \itemize{
#'   \item A character vector with names of reductions
#'   \item \code{NULL} for all reductions
#'   \item \code{NA} or \code{FALSE} for no reductions
#' }
#' \strong{Note}: Only reductions associated with an assay loaded in \code{assays}
#' @param graphs One of:
#' \itemize{
#'   \item A character vector with names of graphs
#'   \item \code{NULL} for all graphs
#'   \item \code{NA} or \code{FALSE} for no reductions
#' }
#' @param meta.data Load object metadata?
#' @param commands Load command information
#' @param misc Load miscellaneous data?
#' @param tools Load tool-specific information?
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return If \code{type} is info, information about the data contained within the
#' h5Seurat file. Otherwise, a \code{Seurat} object with the data requested
#'
#' @export
#'
LoadH5Seurat <- function(file, ...) {
  UseMethod(generic = 'LoadH5Seurat', object = file)
}
