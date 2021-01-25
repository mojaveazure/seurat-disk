#' @include zzz.R
#' @importFrom hdf5r h5attr
#' @importFrom methods slot<- new
#' @importFrom Seurat Cells Key<- Key Cells
#'
NULL

#' Assemble an object from an h5Seurat file
#'
#' @param assay,reduction,graph,image,neighbor,cmd Name of assay, reduction,
#' graph, image, neighbor, or command to load
#' @param file A connected h5Seurat file to pull the data from
#' @param verbose Show progress updates
#'
#' @name AssembleObject
#' @rdname AssembleObject
#'
#' @keywords internal
#'
NULL

#' @param slots Optional vector of assay slots to load, defaults to all slots
#' present in assay
#'
#' @importFrom Seurat CreateAssayObject GetAssayData SetAssayData AddMetaData
#' VariableFeatures<-
#'
#' @return \code{AssembleAssay}: An \code{Assay} object
#'
#' @rdname AssembleObject
#'
#' @aliases AssembleAssay
#'
AssembleAssay <- function(assay, file, slots = NULL, verbose = TRUE) {
  index <- file$index()
  if (!assay %in% names(x = index)) {
    stop("Cannot find assay ", assay, " in this h5Seurat file", call. = FALSE)
  }
  slots.assay <- names(x = Filter(f = isTRUE, x = index[[assay]]$slots))
  slots <- slots %||% slots.assay
  slots <- match.arg(arg = slots, choices = slots.assay, several.ok = TRUE)
  if (!any(c('counts', 'data') %in% slots)) {
    stop("At least one of 'counts' or 'data' must be loaded", call. = FALSE)
  }
  assay.group <- file[['assays']][[assay]]
  features <- FixFeatures(features = assay.group[['features']][])
  # Add counts if not data, otherwise add data
  if ('counts' %in% slots && !'data' %in% slots) {
    if (verbose) {
      message("Initializing ", assay, " with counts")
    }
    counts <- as.matrix(x = assay.group[['counts']])
    rownames(x = counts) <- features
    colnames(x = counts) <- Cells(x = file)
    obj <- CreateAssayObject(counts = counts, min.cells = -1, min.features = -1)
  } else {
    if (verbose) {
      message("Initializing ", assay, " with data")
    }
    data <- as.matrix(x = assay.group[['data']])
    rownames(x = data) <- features
    colnames(x = data) <- Cells(x = file)
    obj <- CreateAssayObject(data = data)
  }
  Key(object = obj) <- Key(object = assay.group)
  # Add remaining slots
  for (slot in slots) {
    if (IsMatrixEmpty(x = GetAssayData(object = obj, slot = slot))) {
      if (verbose) {
        message("Adding ", slot, " for ", assay)
      }
      dat <- as.matrix(x = assay.group[[slot]])
      colnames(x = dat) <- Cells(x = file)
      rownames(x = dat) <- if (slot == 'scale.data') {
        FixFeatures(features = assay.group[['scaled.features']][])
      } else {
        features
      }
      obj <- SetAssayData(object = obj, slot = slot, new.data = dat)
    }
  }
  # Add meta features
  if (assay.group$exists(name = 'meta.features')) {
    if (verbose) {
      message("Adding feature-level metadata for ", assay)
    }
    meta.data <- as.data.frame(
      x = assay.group[['meta.features']],
      row.names = features
    )
    if (ncol(x = meta.data)) {
      obj <- AddMetaData(
        object = obj,
        metadata = meta.data
      )
    }
  }
  # Add variable feature information
  if (assay.group$exists(name = 'variable.features')) {
    if (verbose) {
      message("Adding variable feature information for ", assay)
    }
    VariableFeatures(object = obj) <- assay.group[['variable.features']][]
  }
  # Add miscellaneous information
  if (assay.group$exists(name = 'misc')) {
    if (verbose) {
      message("Adding miscellaneous information for ", assay)
    }
    slot(object = obj, name = 'misc') <- as.list(x = assay.group[['misc']])
  }
  if (assay.group$attr_exists(attr_name = 's4class')) {
    classdef <- unlist(x = strsplit(
      x = h5attr(x = assay.group, which = 's4class'),
      split = ':'
    ))
    pkg <- classdef[1]
    cls <- classdef[2]
    formal <- methods::getClassDef(Class = cls, package = pkg, inherits = FALSE)
    missing <- setdiff(
      x = slotNames(x = formal),
      y = slotNames(x = methods::getClass(Class = 'Assay'))
    )
    missing <- intersect(x = missing, y = names(x = assay.group))
    missing <- sapply(
      X = missing,
      FUN = function(x) {
        return(as.list(x = assay.group[[x]]))
      },
      simplify = FALSE
    )
    obj <- c(SeuratObject::S4ToList(object = obj), missing)
    attr(x = obj, which = 'classDef') <- paste(classdef, collapse = ':')
    obj <- SeuratObject::ListToS4(x = obj)
  }
  return(obj)
}

#' @importClassesFrom Seurat JackStrawData
#' @importFrom Seurat CreateDimReducObject Stdev IsGlobal JS<-
#'
#' @rdname AssembleObject
#'
AssembleDimReduc <- function(reduction, file, verbose = TRUE) {
  index <- file$index()
  index.check <- vapply(
    X = setdiff(x = names(x = index), y = c('global', 'no.assay')),
    FUN = function(x) {
      return(reduction %in% names(x = index[[x]]$reductions))
    },
    FUN.VALUE = logical(length = 1L)
  )
  if (!any(index.check)) {
    stop(
      "Cannot find reduction ",
      reduction,
      " in this h5Seurat file",
      call. = FALSE
    )
  } else if (sum(index.check) > 1) {
    stop("Multiple reductions named ", reduction, call. = FALSE)
  }
  assay <- names(x = which(x = index.check))
  reduc.group <- file[['reductions']][[reduction]]
  key <- Key(object = reduc.group)
  # Pull cell embeddings
  if (index[[assay]]$reductions[[reduction]][['cell.embeddings']]) {
    if (verbose) {
      message("Adding cell embeddings for ", reduction)
    }
    embeddings <- as.matrix(x = reduc.group[['cell.embeddings']])
    rownames(x = embeddings) <- Cells(x = file)
    colnames(x = embeddings) <- paste0(key, 1:ncol(x = embeddings))
  } else {
    if (verbose) {
      warning(
        "No cell embeddings for ",
        reduction,
        call. = FALSE,
        immediate. = TRUE
      )
    }
    embeddings <- new(Class = 'matrix')
  }
  # Pull feature loadings
  if (index[[assay]]$reductions[[reduction]][['feature.loadings']]) {
    if (verbose) {
      message("Adding feature loadings for ", reduction)
    }
    loadings <- as.matrix(x = reduc.group[['feature.loadings']])
    rownames(x = loadings) <- reduc.group[['features']][]
    colnames(x = loadings) <- paste0(key, 1:ncol(x = loadings))
  } else {
    loadings <- new(Class = 'matrix')
  }
  # Pull projected loadings
  if (index[[assay]]$reductions[[reduction]][['feature.loadings.projected']]) {
    if (verbose) {
      message("Adding projected loadings for ", reduction)
    }
    projected <- as.matrix(x = reduc.group[['feature.loadings.projected']])
    rownames(x = projected) <- reduc.group[['projected.features']][]
    colnames(x = projected) <- paste0(key, 1:ncol(x = projected))
  } else {
    projected <- new(Class = 'matrix')
  }
  # Build the object
  obj <- CreateDimReducObject(
    embeddings = embeddings,
    loadings = loadings,
    projected = projected,
    assay = assay,
    stdev = Stdev(object = file, reduction = reduction),
    key = key,
    global = IsGlobal(object = reduc.group)
  )
  # Add misc
  if (reduc.group$exists(name = 'misc')) {
    if (verbose) {
      message("Adding miscellaneous information for ", reduction)
    }
    slot(object = obj, name = 'misc') <- as.list(x = reduc.group[['misc']])
  }
  # Add jackstraw
  if (index[[assay]]$reductions[[reduction]][['jackstraw']]) {
    if (verbose) {
      message("Loading JackStraw data for ", reduction)
    }
    js <- new(Class = 'JackStrawData')
    for (slot in names(x = reduc.group[['jackstraw']])) {
      JS(object = js, slot = slot) <- as.matrix(x = reduc.group[['jackstraw']][[slot]])
    }
    JS(object = obj) <- js
  }
  return(obj)
}

#' @importFrom Seurat as.sparse as.Graph
#'
#' @rdname AssembleObject
#'
AssembleGraph <- function(graph, file, verbose = TRUE) {
  index <- file$index()
  obj <- as.sparse(x = file[['graphs']][[graph]])
  rownames(x = obj) <- colnames(x = obj) <- Cells(x = file)
  obj <- as.Graph(x = obj)
  if (file[['graphs']][[graph]]$attr_exists(attr_name = 'assay.used')) {
    assay <- h5attr(x = file[['graphs']][[graph]], which = 'assay.used')
    if (graph %in% index[[assay]]$graphs) {
      DefaultAssay(object = obj) <- assay
    }
  }
  return(obj)
}

#' @rdname AssembleObject
#'
AssembleImage <- function(image, file, verbose = TRUE) {
  index <- file$index()
  obj <- as.list(x = file[['images']][[image]], row.names = Cells(x = file))
  if (file[['images']][[image]]$attr_exists(attr_name = 'assay')) {
    assay <- h5attr(x = file[['images']][[image]], which = 'assay')
    if (image %in% index[[assay]]$images) {
      DefaultAssay(object = obj) <- assay
    }
  }
  return(obj)
}

#' @importClassesFrom Seurat Neighbor
#'
#' @rdname AssembleObject
#'
AssembleNeighbor <- function(neighbor, file, verbose = TRUE) {
  neighbor.group <- file[['neighbors']][[neighbor]]
  obj <- new(
    Class = 'Neighbor',
    nn.idx =  as.matrix(x = neighbor.group[["nn.idx"]]),
    nn.dist = as.matrix(x = neighbor.group[["nn.dist"]]),
    cell.names =  as.matrix(x = neighbor.group[["cell.names"]])[,1]
  )
  return(obj)
}

#' @importClassesFrom Seurat SeuratCommand
#'
#' @rdname AssembleObject
#'
AssembleSeuratCommand <- function(cmd, file, verbose = TRUE) {
  index <- file$index()
  index.check <- vapply(
    X = setdiff(x = names(x = index), y = 'global'),
    FUN = function(x) {
      return(cmd %in% index[[x]]$commands)
    },
    FUN.VALUE = logical(length = 1L)
  )
  if (!any(index.check)) {
    stop("Cannot find command ", cmd, " in this h5Seurat file", call. = FALSE)
  } else if (sum(index.check) > 1) {
    stop("Multiple commands named ", cmd, call. = FALSE)
  }
  cmd.group <- file[['commands']][[cmd]]
  cmdlog <- new(
    Class = 'SeuratCommand',
    name = h5attr(x = cmd.group, which = 'name'),
    time.stamp = as.POSIXct(x = h5attr(x = cmd.group, which = 'time.stamp')),
    call.string = h5attr(x = cmd.group, 'call.string'),
    params = as.list(x = cmd.group)
  )
  assay <- names(x = which(x = index.check))
  if (assay != 'no.assay') {
    slot(object = cmdlog, name = 'assay.used') <- assay
  }
  return(cmdlog)
}
