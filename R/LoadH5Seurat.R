#' @include zzz.R
#' @include h5Seurat.R
#' @include GetObject.R
#' @include AssembleObject.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Load a saved \code{Seurat} object from an h5Seurat file
#'
#' @param file,x Name of h5Seurat or connected h5Seurat file to load
#' @param assays One of:
#' \itemize{
#'  \item A character vector with names of assays
#'  \item A character vector with one or more of \code{counts}, \code{data},
#'  \code{scale.data} describing which slots of \strong{all assays} to load
#'  \item A named list where each entry is either the name of an assay or a vector
#'  describing which slots (described above) to take from which assay
#'  \item \code{NULL} for all assays
#' }
#' @param reductions One of:
#' \itemize{
#'  \item A character vector with names of reductions
#'  \item \code{NULL} for all reductions
#'  \item \code{NA} for \link[Seurat:IsGlobal]{global} reductions
#'  \item \code{FALSE} for no reductions
#' }
#' \strong{Note}: Only reductions associated with an assay loaded in
#' \code{assays} or marked as \link[Seurat:IsGlobal]{global} will be loaded
#' @param graphs One of:
#' \itemize{
#'  \item A character vector with names of graphs
#'  \item \code{NULL} for all graphs
#'  \item \code{FALSE} for no graphs
#' }
#' \strong{Note}: Only graphs associated with an assay loaded in \code{assays}
#' will be loaded
#' @param images One of:
#' \itemize{
#'  \item A character vector with names of images
#'  \item \code{NULL} for all images
#'  \item \code{NA} for \link[Seurat:IsGlobal]{global} images
#'  \item \code{FALSE} for no images
#' }
#' @param meta.data Load object metadata
#' @param commands Load command information \cr
#' \strong{Note}: only commands associated with an assay loaded in
#' \code{assays} will be loaded
#' @param misc Load miscellaneous data
#' @param tools Load tool-specific information
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return A \code{Seurat} object with the data requested
#'
#' @export
#'
LoadH5Seurat <- function(file, ...) {
  UseMethod(generic = 'LoadH5Seurat', object = file)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat character
#' @export
#'
LoadH5Seurat.character <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  hfile <- h5Seurat$new(filename = file, mode = 'r')
  on.exit(expr = hfile$close_all())
  return(LoadH5Seurat(
    file = hfile,assays = assays,
    reductions = reductions,
    graphs = graphs,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @rdname LoadH5Seurat
#' @method LoadH5Seurat H5File
#' @export
#'
LoadH5Seurat.H5File <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  return(LoadH5Seurat(
    file = as.h5Seurat(x = file),
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @importFrom Seurat as.Seurat
#'
#' @rdname LoadH5Seurat
#' @method LoadH5Seurat h5Seurat
#' @export
#'
LoadH5Seurat.h5Seurat <- function(
  file,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = is.null(x = assays),
  tools = is.null(x = assays),
  verbose = TRUE,
  ...
) {
  return(as.Seurat(
    x = file,
    assays = assays,
    reductions = reductions,
    graphs = graphs,
    images = images,
    meta.data = meta.data,
    commands = commands,
    misc = misc,
    tools = tools,
    verbose = verbose,
    ...
  ))
}

#' @importClassesFrom Seurat Seurat
#' @importFrom methods slot<-
#' @importFrom Seurat as.Seurat DefaultAssay Cells
#' Idents<- Idents Project<- Project
#' AddMetaData
#'
#' @aliases as.Seurat
#'
#' @rdname LoadH5Seurat
#' @method as.Seurat h5Seurat
#' @export
#'
as.Seurat.h5Seurat <- function(
  x,
  assays = NULL,
  reductions = NULL,
  graphs = NULL,
  images = NULL,
  meta.data = TRUE,
  commands = TRUE,
  misc = TRUE,
  tools = TRUE,
  verbose = TRUE,
  ...
) {
  index <- x$index()
  obj.all <- all(vapply(
    X = c(assays, reductions, graphs),
    FUN = is.null,
    FUN.VALUE = logical(length = 1L)
  ))
  # Load Assays
  assays <- GetAssays(assays = assays, index = index)
  if (!DefaultAssay(object = index) %in% names(x = assays)) {
    active.assay <- names(x = assays)[1]
    warning(
      "Default assay not requested, using ",
      active.assay,
      " instead",
      call. = FALSE,
      immediate. = TRUE
    )
  } else {
    active.assay <- DefaultAssay(object = index)
  }
  assay.objects <- vector(mode = 'list', length = length(x = assays))
  names(x = assay.objects) <- names(x = assays)
  for (assay in names(x = assays)) {
    assay.objects[[assay]] <- AssembleAssay(
      assay = assay,
      file = x,
      slots = assays[[assay]],
      verbose = verbose
    )
  }
  default.assay <- list(assay.objects[[active.assay]])
  names(x = default.assay) <- active.assay
  object <- new(
    Class = 'Seurat',
    assays = default.assay,
    active.assay = active.assay,
    meta.data = data.frame(row.names = Cells(x = x)),
    version = package_version(x = x$version())
  )
  for (assay in names(x = assay.objects)) {
    if (assay != active.assay) {
      object[[assay]] <- assay.objects[[assay]]
    }
  }
  # Load DimReducs
  reductions <- GetDimReducs(
    reductions = reductions,
    index = index,
    assays = assays
  )
  for (reduc in reductions) {
    if (verbose) {
      message("Adding reduction ", reduc)
    }
    reduction <- AssembleDimReduc(
      reduction = reduc,
      file = x,
      verbose = verbose
    )
    if (isTRUE(x = getOption(x = 'SeuratDisk.dimreducs.allglobal', default = FALSE))) {
      slot(object = reduction, name = 'global') <- TRUE
    }
    object[[reduc]] <- reduction
  }
  # Load Graphs
  graphs <- GetGraphs(graphs = graphs, index = index, assays = assays)
  for (graph in graphs) {
    if (verbose) {
      message("Adding graph ", graph)
    }
    object[[graph]] <- AssembleGraph(graph = graph, file = x, verbose = verbose)
  }
  # Load SpatialImages
  if (packageVersion(pkg = 'Seurat') >= numeric_version(x = spatial.version)) {
    images <- GetImages(images = images, index = index, assays = assays)
    for (image in images) {
      if (verbose) {
        message("Adding image ", image)
      }
      object[[image]] <- AssembleImage(
        image = image,
        file = x,
        verbose = verbose
      )
    }
  }
  # Load SeuratCommands
  if (commands) {
    if (verbose) {
      message("Adding command information")
    }
    cmds <- GetCommands(index = index, assays = assays)
    cmdlogs <- vector(mode = 'list', length = length(x = cmds))
    names(x = cmdlogs) <- cmds
    for (cmd in cmds) {
      cmdlogs[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = x,
        verbose = verbose
      )
    }
    slot(object = object, name = 'commands') <- cmdlogs
  }
  # Load meta.data
  if (meta.data) {
    if (verbose) {
      message("Adding cell-level metadata")
    }
    md <- as.data.frame(x = x[['meta.data']], row.names = Cells(x = x))
    if (ncol(x = md)) {
      object <- AddMetaData(object = object, metadata = md)
    }
  }
  # Set cell identities and object project
  Idents(object = object) <- Idents(object = x)
  Project(object = object) <- Project(object = x)
  # Load misc
  if (misc) {
    if (verbose) {
      message("Adding miscellaneous information")
    }
    slot(object = object, name = 'misc') <- as.list(x = x[['misc']])
  }
  # Load tools
  if (tools) {
    if (verbose) {
      message("Adding tool-specific results")
    }
    slot(object = object, name = 'tools') <- as.list(x = x[['tools']])
  }
  # Load no.assay information
  if (obj.all && !is.null(x = index$no.assay)) {
    if (verbose) {
      message("Adding data that was not associated with an assay")
    }
    for (graph in index$no.assay$graphs) {
      object[[graph]] <- AssembleGraph(
        graph = graph,
        file = x,
        verbose = verbose
      )
    }
    for (cmd in index$no.assay$commands) {
      object[[cmd]] <- AssembleSeuratCommand(
        cmd = cmd,
        file = x,
        verbose = verbose
      )
    }
  }
  return(object)
}
