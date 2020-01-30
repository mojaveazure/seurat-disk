#' @importFrom hdf5r h5attr h5attr<-
#' @importFrom Seurat GetAssayData Key VariableFeatures Misc Embeddings Loadings
#' DefaultAssay Stdev JS
#' @importFrom methods setOldClass setClassUnion setGeneric setMethod
#' slotNames slot tryNew
#'
NULL

#' Write data to an HDF5 group
#'
#' @param x An object
#' @param name Name to save data as
#' @param h5group An HDF5 file or group (\code{H5File} or \code{H5Group} objects
#' from hdf5r)
#' @param verbose Show progress updates
#'
#' @return Invisibly returns \code{NULL}
#'
#' @export
#'
setGeneric(
  name = 'WriteH5Group',
  def = function(x, name, hgroup, verbose = TRUE) {
    if (!inherits(x = hgroup, what = c('H5File', 'H5Group'))) {
      stop(
        "'hgroup' must be an HDF5 file or group object from hdf5r",
        call. = FALSE
      )
    }
    standardGeneric(f = 'WriteH5Group')
  },
  signature = c('x')
)

setMethod(
  f = 'WriteH5Group',
  signature = c(x = 'ANY'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    if (isS4(x)) {
      xgroup <- hgroup$create_group(name = name)
      h5attr(x = xgroup, which = 's4class') <- ''
      for (i in slotNames(x = x)) {
        WriteH5Group(
          x = slot(object = x, name = i),
          name = i,
          hgroup = xgroup,
          verbose = verbose
        )
      }
    } else if (!is.null(x = x)) {
      hgroup[[name]] <- x
    }
    return(invisible(x = NULL))
  }
)

#' @importClassesFrom Seurat Assay
#'
setMethod(
  f = 'WriteH5Group',
  signature = c('x' = 'Assay'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    xgroup <- hgroup$create_group(name = name)
    # Write out expression data
    # TODO: determine if empty matrices should be present
    for (i in c('counts', 'data', 'scale.data')) {
      dat <- GetAssayData(object = x, slot = i)
      if (!IsMatrixEmpty(x = dat)) {
        if (verbose) {
          message("Adding ", i, " for ", name)
        }
        WriteH5Group(x = dat, name = i, hgroup = xgroup, verbose = verbose)
      }
      # For scale.data, ensure we have the features used
      if (i == 'scale.data') {
        WriteH5Group(
          x = rownames(x = dat),
          name = 'scaled.features',
          hgroup = xgroup,
          verbose = verbose
        )
      }
    }
    # Write out feature names
    WriteH5Group(
      x = rownames(x = x),
      name = 'features',
      hgroup = xgroup,
      verbose = verbose
    )
    # WRite out the key
    h5attr(x = xgroup, which = 'key') <- Key(object = x)
    # Write out variable features
    if (length(x = VariableFeatures(object = x))) {
      if (verbose) {
        message("Adding variable features for ", name)
      }
      WriteH5Group(
        x = VariableFeatures(object = x),
        name = 'variable.features',
        hgroup = xgroup,
        verbose = verbose
      )
    } else if (verbose) {
      message("No variable features found for ", name)
    }
    # Write out meta.features
    if (ncol(x = x[[]])) {
      if (verbose) {
        message("Adding feature-level metadata for ", name)
      }
      WriteH5Group(
        x = x[[]],
        name = 'meta.features',
        hgroup = xgroup,
        verbose = verbose
      )
    } else if (verbose) {
      message("No feature-level metadata found for ", name)
    }
    # Write out miscellaneous data
    WriteH5Group(
      x = Misc(object = x),
      name = 'misc',
      hgroup = xgroup,
      verbose = verbose
    )
    # Write out other slots for extended assay objects
    if (class(x = x)[1] != 'Assay') {
      h5attr(x = xgroup, which = 's4class') <- GetClass(class = class(x = x))
      slots.extended <- setdiff(
        x = slotNames(x = x),
        y = slotNames(x = tryNew(Class = 'Assay'))
      )
      for (slot in slots.extended) {
        if (verbose) {
          message("Writing out ", slot, " for ", name)
        }
        WriteH5Group(
          x = slot(object = x, name = slot),
          name = slot,
          hgroup = xgroup,
          verbose = verbose
        )
      }
    }
    return(invisible(x = NULL))
  }
)

setMethod(
  f = 'WriteH5Group',
  signature = c('x' = 'data.frame'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    hgroup[[name]] <- x
    return(invisible(x = NULL))
  }
)

setMethod(
  f = 'WriteH5Group',
  signature = c('x' = 'dgCMatrix'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    xgroup <- hgroup$create_group(name = name)
    xgroup[['indices']] <- slot(object = x, name = 'i')
    xgroup[['indptr']] <- slot(object = x, name = 'p')
    xgroup[['data']] <- slot(object = x, name = 'x')
    return(invisible(x = NULL))
  }
)

#' @importClassesFrom Seurat DimReduc
#'
setMethod(
  f = 'WriteH5Group',
  signature = c('x' = 'DimReduc'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    xgroup <- hgroup$create_group(name = name)
    # Add cell embeddings
    if (verbose) {
      message("Adding cell embeddings for ", name)
    }
    WriteH5Group(
      x = Embeddings(object = x),
      name = 'cell.embeddings',
      hgroup = xgroup,
      verbose = verbose
    )
    # Add feature loadings
    for (i in c('feature.loadings', 'feature.loadings.projected')) {
      projected <- grepl(pattern = 'projected', x = i)
      type <- ifelse(test = projected, yes = 'projected loadings', no = 'loadings')
      if (!IsMatrixEmpty(x = Loadings(object = x, projected = projected))) {
        if (verbose) {
          message("Adding ", type, " for ", name)
        }
        loadings <- Loadings(object = x, projected = projected)
        WriteH5Group(x = loadings, name = i, hgroup = xgroup, verbose = verbose)
        WriteH5Group(
          x = rownames(x = loadings),
          name = ifelse(test = projected, yes = 'projected.features', no = 'features'),
          hgroup = xgroup,
          verbose = verbose
        )
      } else if (verbose) {
        message("No ", type, " for ", name)
      }
    }
    # Add assay and key
    hdf5r::h5attr(x = hgroup, which = 'active.assay') <- DefaultAssay(object = x)
    hdf5r::h5attr(x = hgroup, which = 'key') <- Key(object = x)
    # Add standard deviations
    if (length(x = Stdev(object = x)) > 0) {
      if (verbose) {
        message("Adding standard deviations for ", name)
      }
      WriteH5Group(
        x = Stdev(object = x),
        name = 'stdev',
        hgroup = xgroup,
        verbose = verbose
      )
    } else if (verbose) {
      message("No standard deviations for ", name)
    }
    # Add JackStraw
    if (as.logical(x = JS(object = x))) {
      if (verbose) {
        message("Adding JackStraw information for ", name)
      }
      WriteH5Group(
        x = JS(object = x),
        name = 'jackstraw',
        hgroup = xgroup,
        verbose = verbose
      )
    } else if (verbose) {
      message("No JackStraw data for ", name)
    }
    # Add misc
    WriteH5Group(
      x = slot(object = x, name = 'misc'),
      name = 'misc',
      hgroup = xgroup,
      verbose = verbose
    )
    return(invisible(x = NULL))
  }
)

setMethod(
  f = 'WriteH5Group',
  signature = c('x' = 'factor'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    xgroup <- hgroup$create_group(name = name)
    # Write the integer values out
    WriteH5Group(
      x = as.integer(x = x),
      name = 'values',
      hgroup = xgroup,
      verbose = verbose
    )
    # Write the levels out
    WriteH5Group(
      x = levels(x = x),
      name = 'levels',
      hgroup = xgroup,
      verbose = verbose
    )
    return(invisible(x = NULL))
  }
)

setMethod(
  f = 'WriteH5Group',
  signature = c('x' = 'list'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    xgroup <- hgroup$create_group(name = name)
    for (i in Enumerate(x = x)) {
      WriteH5Group(x = i$value, name = i$name, hgroup = xgroup)
    }
    return(invisible(x = NULL))
  }
)

#' @importClassesFrom Seurat SeuratCommand
#'
setMethod(
  f = 'WriteH5Group',
  signature = c('x' = 'SeuratCommand'),
  definition = function(x, name, hgroup, verbose = TRUE) {
    # Write out params
    WriteH5Group(
      x = Filter(
        f = Negate(f = is.function),
        x = as.list(x = x)
      ),
      name = name,
      hgroup = hgroup,
      verbose = verbose
    )
    # Add extra information as HDF5 attributes
    for (slot in slotNames(x = x)) {
      if (slot == 'params') {
        next
      }
      h5attr(x = hgroup[[name]], which = slot) <- as.character(x = slot(
        object = x,
        name = slot
      ))
    }
    return(invisible(x = NULL))
  }
)
