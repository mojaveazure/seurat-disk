#' @include generics.R
#'
NULL

#' @rdname as.h5Seurat
#' @method as.h5Seurat H5File
#' @export
#'
as.h5Seurat.H5File <- function(x, ...) {
  return(h5Seurat$new(
    filename = x$filename,
    mode = ifelse(test = x$mode == 'r', yes = 'r', no = 'r+')
  ))
}

#' @param filename Name of file to save \code{x} to
#' @param overwrite Overwrite \code{filename} if present?
#' @param verbose Show progress updates
#'
#' @importFrom tools file_ext
#' @importFrom Seurat Project Assays Reductions DefaultAssay
#' Idents Command Misc Tool
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
  if (!grepl(pattern = '^h5seurat$', x = file_ext(x = filename), ignore.case = TRUE)) {
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
  # Add Assays
  for (assay in Assays(object = x)) {
    WriteH5Group(
      x = x[[assay]],
      name = assay,
      hgroup = hfile[['assays']],
      verbose = verbose
    )
  }
  # Add DimReducs
  for (reduc in Reductions(object = x)) {
    WriteH5Group(
      x = x[[reduc]],
      name = reduc,
      hgroup = hfile[['reductions']],
      verbose = verbose
    )
  }
  # Add Graphs
  graphs <- Filter(
    f = function(g) {
      return(inherits(x = x[[g]], what = 'Graph'))
    },
    x = names(x = x)
  )
  for (graph in graphs) {
    WriteH5Group(
      x = x[[graph]],
      name = graph,
      hgroup = hfile[['graphs']],
      verbose = verbose
    )
  }
  # Add attributes for project, default assay, and version
  Project(object = hfile) <- Project(object = x)
  DefaultAssay(object = hfile) <- DefaultAssay(object = x)
  h5attr(x = hfile, which = 'version') <- as.character(x = slot(
    object = x,
    name = 'version'
  ))
  # TODO: Add Images
  if (package_version(x = h5attr(x = hfile, which = 'version')) >= package_version(x = '3.2.0')) {
    warning(
      "Support for spatial image data is not yet implemented",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  # Add metadata, cell names, and identity classes4
  hfile[['meta.data']] <- x[[]]
  hfile[['cell.names']] <- colnames(x = x)
  hfile[['active.ident']] <- Idents(object = x)
  # Add SeuratCommands
  # browser()
  for (cmd in Command(object = x)) {
    WriteH5Group(
      x = x[[cmd]],
      name = cmd,
      hgroup = hfile[['commands']],
      verbose = verbose
    )
  }
  # Add miscellaneous data
  if (length(x = Misc(object = x))) {
    WriteH5Group(
      x = Misc(object = x),
      name = 'misc',
      hgroup = hfile,
      verbose = verbose
    )
  }
  # Add tool data
  for (tool in Tool(object = x)) {
    WriteH5Group(
      x = Tool(object = x, slot = tool),
      name = tool,
      hgroup = hfile[['tools']],
      verbose = verbose
    )
  }
  return(hfile)
}

#' Save an object as an h5seurat file
#'
#' @inheritParams as.h5Seurat
#' @param filename Name of file to save \code{object} to
#'
#' @return Invisibly returns h5seurat file name
#'
#' @seealso \code{\link{as.h5Seurat}}
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
