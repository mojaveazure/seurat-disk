#' @include AssembleObject.R
#' @include GetObject.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Append data from an h5Seurat file to a preexisting
#' \code{\link[Seurat]{Seurat}} object
#'
#' @inheritParams LoadH5Seurat
#' @param object A \code{\link[Seurat]{Seurat}} object to append data to
#' @param assays One of:
#' \itemize{
#'   \item A character vector with names of assays
#'   \item A character vector with one or more of \code{counts}, \code{data},
#'   \code{scale.data} describing which slots of \strong{all assays} to load
#'   \item A named list where each entry is either the name of an assay or a vector
#'   describing which slots (described above) to take from which assay
#'   \item \code{NULL} for all assays
#'   \item \code{FALSE} for no assays
#' }
#' @param extras Extra information to load; supports any combination of the
#' following values:
#' \describe{
#'  \item{\dQuote{commands}}{Load command logs. If \code{overwrite = TRUE},
#'  replaces existing command logs}
#' }
#' @param overwrite Overwrite existing data in \code{object} with data from
#' \code{file}
#'
#' @return \code{object} with the extra data requested
#'
#' @export
#'
AppendData <- function(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = 'commands',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  if (!inherits(x = object, what = 'Seurat')) {
    stop("'object' must be a Seurat object", call. = FALSE)
  }
  UseMethod(generic = 'AppendData', object = file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname AppendData
#' @method AppendData character
#' @export
#'
AppendData.character <- function(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = 'commands',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  if (!file.exists(file)) {
    stop("Cannot find h5Seurat file ", file, call. = FALSE)
  }
  hfile <- h5Seurat$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(AppendData(
    file = hfile,
    object = object,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @rdname AppendData
#' @method AppendData H5File
#' @export
#'
AppendData.H5File <- function(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = 'commands',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  return(AppendData(
    file = as.h5Seurat(x = file),
    object = object,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    overwrite = overwrite,
    verbose = verbose,
    ...
  ))
}

#' @importFrom Seurat Cells Assays Reductions Command
#'
#' @rdname AppendData
#' @method AppendData h5Seurat
#' @export
#'
AppendData.h5Seurat <- function(
  file,
  object,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  extras = 'commands',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  index <- file$index()
  if (!all(Cells(x = file) == Cells(x = object))) {
    stop(
      "Mismatched cells between the h5Seurat file and the Seurat object",
      call. = FALSE
    )
  }
  extras <- match.arg(arg = extras, several.ok = TRUE)
  obj.all <- all(vapply(
    X = c(assays, reductions, graphs),
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  ))
  # Load Assays
  if (!isFALSE(x = assays)) {
    assays <- GetAssays(assays = assays, index = index)
    if (!overwrite) {
      assays <- assays[!names(x = assays) %in% Assays(object = object)]
    } else if (verbose) {
      overwritten <- intersect(x = names(x = assays), y = Assays(object = object))
      if (length(x = overwritten)) {
        message(
          "Overwriting the following assays: ",
          paste(overwritten, collapse = ', ')
        )
      }
    }
    if (length(x = assays) && verbose) {
      message("Adding data for ", length(x = assays), " assays")
    }
    for (assay in names(x = assays)) {
      object[[assay]] <- AssembleAssay(
        assay = assay,
        file = file,
        slots = assays[[assay]],
        verbose = verbose
      )
    }
  }
  # Load DimReducs
  reductions <- GetDimReducs(
    reductions = reductions,
    index = index,
    assays = Assays(object = object)
  )
  if (!overwrite) {
    reductions <- setdiff(x = reductions, y = Reductions(object = object))
  } else if (verbose) {
    overwritten <- intersect(x = reductions, y = Reductions(object = object))
    if (length(x = overwritten)) {
      message(
        "Overwriting the following dimensional reductions: ",
        paste(overwritten, collapse = ', ')
      )
    }
  }
  if (length(x = reductions) && verbose) {
    message("Adding data for ", length(x = reductions), " dimensional reductions")
  }
  for (reduc in reductions) {
    object[[reduc]] <- AssembleDimReduc(
      reduction = reduc,
      file = file,
      verbose = verbose
    )
  }
  # Load Graphs
  graphs <- GetGraphs(
    graphs = graphs,
    index = index,
    assays = Assays(object = object)
  )
  object.graphs <- Filter(
    f = function(x) {
      return(inherits(x = object[[x]], what = 'Graph'))
    },
    x = names(x = graphs)
  )
  if (!overwrite) {
    graphs <- setdiff(x = graphs, y = object.graphs)
  } else if (verbose) {
    overwritten <- intersect(x = graphs, object.graphs)
    if (length(x = overwritten)) {
      message(
        "Overwriting the following nearest-neighbor graphs: ",
        paste(overwritten, collapse = ', ')
      )
    }
  }
  if (length(x = graphs)) {
    message("Adding data for ", length(x = graphs), " nearest-neighbor graphs")
  }
  for (graph in graphs) {
    object[[graph]] <- AssembleGraph(
      graph = graph,
      file = file,
      verbose = verbose
    )
  }
  # TODO: Load Neighbors
  # TODO: Load SpatialImages
  # Load SeuratCommands
  if ('commands' %in% extras) {
    commands <- GetCommands(index = index, assays = Assays(object = object))
    if (!overwrite) {
      commands <- setdiff(x = commands, y = Command(object = object))
    } else if (verbose) {
      overwritten <- intersect(x = commands, y = Command(object = object))
      if (length(x = overwritten)) {
        message(
          "The following command logs will be overwritten: ",
          paste(overwritten, collapse = ', ')
        )
      }
    }
    if (length(x = commands) && verbose) {
      message("Adding ", length(x = commands), " command logs")
    }
    for (cmd in commands) {
      object[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = file,
        verbose = verbose
      )
    }
  }
  # TODO: Load meta.data
  # TODO: Load misc
  # TODO: Load tools
  # Load no.assay information
  if (obj.all && !is.null(x = index$no.assay)) {
    if (verbose) {
      message("Adding data that was not associated with an assay")
    }
    for (graph in index$no.assay$graphs) {
      object[[graph]] <- AssembleGraph(
        graph = graph,
        file = file,
        verbose = verbose
      )
    }
    if ('commands' %in% extras) {
      for (cmd in index$no.assay$commands) {
        object[[cmd]] <- AssembleSeuratCommand(
          cmd = cmd,
          file = file,
          verbose = verbose
        )
      }
    }
  }
  return(object)
}
