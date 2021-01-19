#' @include WriteH5Group.R
#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Save a \code{Seurat} object to an h5Seurat file
#'
#' @param object,x An object
#' @param filename Name of file to save the object to
#' @param overwrite Overwrite \code{filename} if present
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return \code{SaveH5Seurat}: Invisbly returns \code{filename}
#'
#' @name SaveH5Seurat
#' @rdname SaveH5Seurat
#'
#' @export
#'
SaveH5Seurat <- function(
                         object,
                         filename,
                         overwrite = FALSE,
                         verbose = TRUE,
                         ...) {
  UseMethod(generic = "SaveH5Seurat", object = object)
}

#' @return \code{as.h5Seurat}: An \code{\link{h5Seurat}} object
#'
#' @rdname SaveH5Seurat
#'
#' @export
#'
as.h5Seurat <- function(x, ...) {
  UseMethod(generic = "as.h5Seurat", object = x)
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom Seurat as.Seurat Project
#'
#' @rdname SaveH5Seurat
#' @method SaveH5Seurat default
#' @export
#'
SaveH5Seurat.default <- function(
                                 object,
                                 filename,
                                 overwrite = FALSE,
                                 verbose = TRUE,
                                 ...) {
  tryCatch(
    expr = object <- as.Seurat(object = object, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        paste(class(x = object), collapse = ", "),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = object), ".h5Seurat")
  }
  return(invisible(x = SaveH5Seurat(
    object = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )))
}

#' @importFrom Seurat Project
#'
#' @rdname SaveH5Seurat
#' @method SaveH5Seurat Seurat
#' @export
#'
SaveH5Seurat.Seurat <- function(
                                object,
                                filename = paste0(Project(object = object), ".h5Seurat"),
                                overwrite = FALSE,
                                verbose = TRUE,
                                ...) {
  h5seurat <- as.h5Seurat(
    x = object,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  h5seurat$close_all()
  return(invisible(x = h5seurat$filename))
}

#' @importFrom Seurat as.Seurat Project
#'
#' @rdname SaveH5Seurat
#' @method as.h5Seurat default
#' @export
#'
as.h5Seurat.default <- function(
                                x,
                                filename,
                                overwrite = FALSE,
                                verbose = TRUE,
                                ...) {
  tryCatch(
    expr = x <- as.Seurat(x = x, verbose = verbose, ...),
    error = function(...) {
      stop(
        "Unable to coerce an object of class ",
        paste(class(x = x), collapse = ", "),
        " to a Seurat object",
        call. = FALSE
      )
    }
  )
  if (missing(x = filename)) {
    filename <- paste0(Project(object = x), ".h5Seurat")
  }
  return(SaveH5Seurat(
    object = x,
    filename = filename,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname SaveH5Seurat
#' @method as.h5Seurat H5File
#' @export
#'
as.h5Seurat.H5File <- function(x, ...) {
  return(h5Seurat$new(
    filename = x$filename,
    mode = ifelse(test = Writeable(x = x, error = FALSE), yes = "r", no = "r+")
  ))
}

#' @importFrom tools file_ext
#' @importFrom Seurat Project Assays Reductions DefaultAssay<- DefaultAssay
#' Idents Command Misc Tool 
#' @importFrom hdf5r h5attr
#'
#' @rdname SaveH5Seurat
#' @method as.h5Seurat Seurat
#' @export
#'
as.h5Seurat.Seurat <- function(
                               x,
                               filename = paste0(Project(object = x), ".h5seurat"),
                               overwrite = FALSE,
                               verbose = TRUE,
                               ...) {
  if (!grepl(pattern = "^h5seurat$", x = file_ext(x = filename), ignore.case = TRUE)) {
    filename <- paste0(filename, ".h5seurat")
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
  hfile <- h5Seurat$new(filename = filename, mode = "w")
  # Set the version
  object.version <- as.character(x = slot(object = x, name = "version"))
  tryCatch(
    expr = hfile$set.version(version = object.version),
    error = function(err) {
      file.remove(hfile$filename)
      test.message <- paste0(
        "All target versions greater than query version (",
        object.version,
        ")"
      )
      if (err$message == test.message) {
        stop(
          "Object too old to save, please update your Seurat object to at least v3.1.2 using UpdateSeuratObject",
          call. = FALSE
        )
      }
      stop(err$message, call. = FALSE)
    }
  )
  # Add Assays
  for (assay in Assays(object = x)) {
    hfile <- WriteH5Group(
      x = x[[assay]],
      name = assay,
      group = "/assays",
      hfile = hfile,
      verbose = verbose
    )
  }
  # Add DimReducs
  for (reduc in Reductions(object = x)) {
    hfile <- WriteH5Group(
      x = x[[reduc]],
      name = reduc,
      group = "/reductions",
      hfile = hfile,
      verbose = verbose
    )
  }
  # Add Neighbors
  neighbors <- names(x = slot(object = x, name = "neighbors"))
  for (neighbor in neighbors) {
    hfile <- WriteH5Group(
      x = x[[neighbor]],
      name = neighbor,
      group = "/neighbors",
      hfile = hfile,
      verbose = verbose
    )
  }
  # Add Graphs
  graphs <- Filter(
    f = function(g) {
      return(inherits(x = x[[g]], what = "Graph"))
    },
    x = names(x = x)
  )
  for (graph in graphs) {
    hfile <- WriteH5Group(
      x = x[[graph]],
      name = graph,
      group = "/graphs",
      hfile = hfile,
      verbose = verbose
    )
  }
  # Add attributes for project and default assay
  h5attr(x = hfile, which = "project") <- Project(object = x)
  h5attr(x = hfile, which = "active.assay") <- DefaultAssay(object = x)
  
  # Add Images
  if (package_version(x = object.version) >= package_version(x = spatial.version)) {
    # Older versions of Seurat don't have Images, call directly instead
    for (image in Seurat::Images(object = x)) {
      if (verbose) {
        message("Adding image ", image)
      }
      hfile <- WriteH5Group(
        x = x[[image]],
        name = image,
        group = "/images",
        hfile = hfile,
        verbose = verbose
      )
    }
  }
  # Add metadata, cell names, and identity classes
  hfile <- WriteH5Group(
    x = x[[]],
    name = "meta.data",
    group = "/",
    hfile = hfile,
    verbose = verbose
  )
  hfile <- WriteH5Group(
    x = colnames(x = x),
    name = "cell.names",
    group = "/",
    hfile = hfile,
    verbose = verbose
  )
  hfile <- WriteH5Group(
    x = Idents(object = x),
    name = "active.ident",
    group = "/",
    hfile = hfile,
    verbose = verbose
  )
  # Add SeuratCommands
  for (cmd in Command(object = x)) {
    hfile <- WriteH5Group(
      x = x[[cmd]],
      name = cmd,
      group = "/commands",
      hfile = hfile,
      verbose = verbose
    )
  }
  # Add miscellaneous data
  for (misc in names(x = Misc(object = x))) {
    hfile <- WriteH5Group(
      x = Misc(object = x, slot = misc),
      name = misc,
      group = "/misc",
      hfile = hfile,
      verbose = verbose
    )
  }
  # Add tool data
  for (tool in Tool(object = x)) {
    hfile <- WriteH5Group(
      x = Tool(object = x, slot = tool),
      name = tool,
      group = "/tools",
      hfile = hfile,
      verbose = verbose
    )
  }
  return(hfile)
}
