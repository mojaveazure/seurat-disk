#' @include zzz.R
#' @importFrom methods setOldClass setClassUnion setGeneric setMethod slotNames
#' slot tryNew
#' @importFrom Seurat GetAssayData Key VariableFeatures Misc Embeddings Loadings
#' DefaultAssay IsGlobal Stdev JS
#'
NULL

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Internal utility functions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write lists and other data to an HDF5 dataset
#'
#' @inheritParams WriteH5Group
#'
#' @return Invisibly returns \code{NULL}
#'
#' @keywords internal
#'
BasicWrite <- function(x, name, group, hfile = hfile, verbose = TRUE) {
  if (is.data.frame(x = x)) {
    WriteH5Group(x = x, name = name, group = group, hfile = hfile, verbose = verbose)
  } else if (is.list(x = x)) {
    x <- PadNames(x = x)
    if (group != "/") {
      hgroup <- hfile[[group]]
    } else {
      hgroup <- hfile
    }
    xgroup <- hgroup$create_group(name = name)
    for (i in seq_along(along.with = x)) {
      WriteH5Group(
        x = x[[i]],
        name = names(x = x)[i],
        group = paste0(group, "/", name),
        hfile = hfile,
        verbose = verbose
      )
    }
    if (!is.null(x = names(x = x)) && length(x = names(x = x))) {
      xgroup$create_attr(
        attr_name = "names",
        robj = intersect(x = names(x = x), y = names(x = xgroup)),
        dtype = GuessDType(x = names(x = x)[1])
      )
    }
    if (!all(class(x = x) == "list")) {
      xgroup$create_attr(
        attr_name = "s3class",
        robj = class(x = x),
        dtype = GuessDType(x = class(x = x)[1])
      )
    }
  } else if (!is.null(x = x)) {
    hgroup$create_dataset(name = name, robj = x, dtype = GuessDType(x = x))
  }
  return(invisible(x = NULL))
}

#' Write a SpatialImage object to an HDF5 dataset
#'
#' @inheritParams WriteH5Group
#'
#' @return Invisibly returns \code{NULL}
#'
#' @keywords internal
#'
ImageWrite <- function(x, name, group, hfile, verbose = TRUE) {
  if (!inherits(x = x, what = "SpatialImage")) {
    stop(
      "'ImageWrite' work only for SpatialImage-derived objects",
      call. = FALSE
    )
  }
  if (group != "/") {
    hgroup <- hfile[[group]]
  } else {
    hgroup <- hfile
  }
  xgroup <- hgroup$create_group(name = name)
  # Add assay, globality, and class information
  xgroup$create_attr(
    attr_name = "assay",
    robj = DefaultAssay(object = x),
    dtype = GuessDType(x = DefaultAssay(object = x))
  )
  xgroup$create_attr(
    attr_name = "global",
    robj = BoolToInt(x = IsGlobal(object = x)),
    dtype = GuessDType(x = IsGlobal(object = x))
  )
  xgroup$create_attr(
    attr_name = "s4class",
    robj = GetClass(class = class(x = x)[1]),
    dtype = GuessDType(x = GetClass(class = class(x = x)[1]))
  )
  # Write out slots other than assay
  slots <- setdiff(x = slotNames(x = x), y = c("assay", "global"))
  for (slot in slots) {
    WriteH5Group(
      x = slot(object = x, name = slot),
      name = slot,
      group = paste0(group, "/", name),
      hfile = hfile,
      verbose = verbose
    )
  }
  return(invisible(x = NULL))
}

#' Write a sparse matrix to an HDF5 dataset
#'
#' @inheritParams WriteH5Group
#'
#' @return Invisibly returns \code{NULL}
#'
#' @keywords internal
#'
#' @import HDF5Array
#'
SparseWrite <- function(x, name, group, hfile, verbose = TRUE) {
  if (group != "/") {
    hgroup <- hfile[[group]]
  } else {
    hgroup <- hfile
  }
  filename <- hfile$get_filename()
  hfile$close_all()
  writeTENxMatrix(x = x, filepath = filename, group = name, verbose = FALSE)
  hfile <- hdf5r::H5File$new(filename = filename, mode = "r+")
  xgroup <- hfile[[name]]
  xgroup$create_attr(
    attr_name = "dims",
    robj = dim(x = x),
    dtype = GuessDType(dim(x = x))
  )
  xgroup$close()
  assign("hgroup", hgroup, envir = .BaseNamespaceEnv)
  return(invisible(x = NULL))
}

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# WriteH5Group generic
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Write data to an HDF5 group
#'
#' Writing data to HDF5 files can be done simply with usually sensible defaults.
#' However, when wanting any semblance of control over how an R object is
#' written out, the code constructs get complicated quickly. \code{WriteH5Group}
#' provides a wrapper with sensible defaults over some of these complex code
#' constructs to provide greater control over how data are written to disk.
#' These defaults were chosen to fit best with \link{h5Seurat} files, see
#' \code{\href{../doc/h5Seurat-spec.html}{vignette("h5Seurat-spec")}} for more
#' details
#'
#' @param x An object
#' @param name Name to save data as
#' @param hgroup An HDF5 file or group (\code{H5File} or \code{H5Group} objects
#' from hdf5r)
#' @param verbose Show progress updates
#'
#' @return Invisibly returns \code{NULL}
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Setup an HDF5 file
#' hfile <- hdf5r::H5File$new(filename = tempfile(fileext = ".h5"), mode = "a")
#' }
#'
setGeneric(
  name = "WriteH5Group",
  def = function(x, name, group, hfile, verbose = TRUE) {
    if (!inherits(x = hfile, what = "H5File")) {
      stop(
        "'hfile' must be an HDF5 file object from hdf5r",
        call. = FALSE
      )
    }
    standardGeneric(f = "WriteH5Group")
  },
  signature = c("x")
)

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# WriteH5Group definitions
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom SeuratObject S4ToList
#'
#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c(x = "ANY"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    if (group != "/") {
      tryCatch(
        expr = {
          hgroup <- hfile[[group]]
        },
        error = function(...) {
          stop(
            "Unable to find the group [", group, "] in the HDF5 file",
            call. = FALSE
          )
        }
      )
    }
    if (inherits(x = x, what = "SpatialImage")) {
      ImageWrite(x = x, name = name, group = group, hfile = hfile, verbose = verbose)
    } else if (isS4(x)) {
      if (group != "/") {
        hgroup <- hfile[[group]]
      } else {
        hgroup <- hfile
      }
      xgroup <- hgroup$create_group(name = name)
      classdef <- attr(x = S4ToList(object = x), which = "classDef")
      xgroup$create_attr(
        attr_name = "s4class",
        robj = classdef,
        dtype = GuessDType(x = classdef)
      )
      # class <- GetClass(class = class(x = x))
      # xgroup$create_attr(
      #   attr_name = 's4class',
      #   robj = class,
      #   dtype = GuessDType(x = class)
      # )
      for (i in slotNames(x = x)) {
        WriteH5Group(
          x = slot(object = x, name = i),
          name = i,
          group = paste0(group, "/", name),
          hfile = hfile,
          verbose = verbose
        )
      }
    } else if (!is.null(x = x)) {
      BasicWrite(x = x, name = name, group = group, hfile = hfile, verbose = verbose)
    }
    return(invisible(x = NULL))
  }
)

#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "array"),
  definition = BasicWrite
)

#' @importClassesFrom Seurat Assay
#'
#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "Assay"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    if (group != "/") {
      hgroup <- hfile[[group]]
    } else {
      hgroup <- hfile
    }
    xgroup <- hgroup$create_group(name = name)
    # Write out expression data
    # TODO: determine if empty matrices should be present
    for (i in c("counts", "data", "scale.data")) {
      dat <- GetAssayData(object = x, slot = i)
      if (!IsMatrixEmpty(x = dat)) {
        if (verbose) {
          message("Adding ", i, " for ", name)
        }
        WriteH5Group(x = dat, name = i, group = paste0(group, "/", name), hfile = hfile, verbose = verbose)
      }
      # For scale.data, ensure we have the features used
      if (i == "scale.data") {
        WriteH5Group(
          x = rownames(x = dat),
          name = "scaled.features",
          group = paste0(group, "/", name),
          hfile = hfile,
          verbose = verbose
        )
      }
    }
    # Write out feature names
    WriteH5Group(
      x = rownames(x = x),
      name = "features",
      group = paste0(group, "/", name),
      hfile = hfile,
      verbose = verbose
    )
    # Write out the key
    xgroup$create_attr(
      attr_name = "key",
      robj = Key(object = x),
      dtype = GuessDType(x = Key(object = x))
    )
    # Write out variable features
    if (length(x = VariableFeatures(object = x))) {
      if (verbose) {
        message("Adding variable features for ", name)
      }
      WriteH5Group(
        x = VariableFeatures(object = x),
        name = "variable.features",
        group = paste0(group, "/", name),
        hfile = hfile,
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
        name = "meta.features",
        group = paste0(group, "/", name),
        hfile = hfile,
        verbose = verbose
      )
    } else if (verbose) {
      message("No feature-level metadata found for ", name)
    }
    # Write out miscellaneous data
    WriteH5Group(
      x = Misc(object = x),
      name = "misc",
      group = paste0(group, "/", name),
      hfile = hfile,
      verbose = verbose
    )
    # Write out other slots for extended assay objects
    if (class(x = x)[1] != "Assay") {
      # extclass <- GetClass(class = class(x = x))
      extclass <- attr(x = SeuratObject::S4ToList(object = x), which = "classDef")
      xgroup$create_attr(
        attr_name = "s4class",
        robj = extclass,
        dtype = GuessDType(x = extclass)
      )
      slots.extended <- setdiff(
        x = slotNames(x = x),
        # y = slotNames(x = tryNew(Class = 'Assay'))
        y = slotNames(x = methods::getClassDef(Class = "Assay"))
      )
      for (slot in slots.extended) {
        if (verbose) {
          message("Writing out ", slot, " for ", name)
        }
        WriteH5Group(
          x = slot(object = x, name = slot),
          name = slot,
          group = paste0(group, "/", name),
          hfile = hfile,
          verbose = verbose
        )
      }
    }
    return(invisible(x = NULL))
  }
)

#' @rdname WriteH5Group
#'
#' @examples
#' \donttest{
#' # Data frames are stored as either datasets or groups, depending on the
#' # presence of factor columns
#' df <- data.frame(
#'   x = c("g1", "g1", "g2", "g1", "g2"),
#'   y = 1:5,
#'   stringsAsFactors = FALSE
#' )
#'
#' # When no factor columns are present, the data frame is written as a single
#' # HDF5 compound dataset
#' WriteH5Group(x = df, name = "df", hgroup = hfile)
#' hfile[["df"]]
#'
#' # When factors are present, the data frame is written as a group
#' # This is because h5py does not implement HDF5 Enums, so factor level
#' # information would be lost
#' df$x <- factor(x = df$x)
#' WriteH5Group(x = df, name = "df.factor", hgroup = hfile)
#' hfile[["df.factor"]]
#' }
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "data.frame"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    factor.cols <- vapply(
      X = colnames(x = x),
      FUN = function(i) {
        return(is.factor(x = x[, i, drop = TRUE]))
      },
      FUN.VALUE = logical(length = 1L),
      USE.NAMES = FALSE
    )
    bool.cols <- vapply(
      X = colnames(x = x),
      FUN = function(i) {
        return(is.logical(x = x[, i, drop = TRUE]))
      },
      FUN.VALUE = logical(length = 1L),
      USE.NAMES = FALSE
    )
    if (any(factor.cols) || getOption(x = "SeuratDisk.dtypes.dataframe_as_group", default = FALSE)) {
      if (group != "/") {
        hgroup <- hfile[[group]]
      } else {
        hgroup <- hfile
      }
      xgroup <- hgroup$create_group(name = name)
      for (i in colnames(x = x)) {
        WriteH5Group(
          x = x[, i, drop = TRUE],
          name = i,
          group = paste0(group, "/", name),
          hfile = hfile,
          verbose = verbose
        )
      }
      xgroup$create_attr(
        attr_name = "colnames",
        robj = intersect(x = colnames(x = x), y = names(x = xgroup)),
        dtype = GuessDType(x = colnames(x = x))
      )
      if (any(bool.cols)) {
        xgroup$create_attr(
          attr_name = "logicals",
          robj = intersect(x = colnames(x = x)[bool.cols], y = names(x = xgroup)),
          dtype = GuessDType(x = colnames(x = x))
        )
      }
    } else {
      for (i in colnames(x = x)[bool.cols]) {
        x[[i]] <- BoolToInt(x = x[[i]])
      }
      hgroup$create_dataset(name = name, robj = x, dtype = GuessDType(x = x))
      if (any(bool.cols)) {
        hgroup[[name]]$create_attr(
          attr_name = "logicals",
          robj = intersect(
            x = colnames(x = x)[bool.cols],
            y = hgroup[[name]]$get_type()$get_cpd_labels()
          ),
          dtype = GuessDType(x = colnames(x = x))
        )
      }
    }
    return(invisible(x = NULL))
  }
)

#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "dgCMatrix"),
  definition = SparseWrite
)

#' @importClassesFrom Seurat DimReduc
#'
#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "DimReduc"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    if (group != "/") {
      hgroup <- hfile[[group]]
    } else {
      hgroup <- hfile
    }
    xgroup <- hgroup$create_group(name = name)
    # Add cell embeddings
    if (verbose) {
      message("Adding cell embeddings for ", name)
    }
    WriteH5Group(
      x = Embeddings(object = x),
      name = "cell.embeddings",
      group = paste0(group, "/", name),
      hfile = hfile,
      verbose = verbose
    )
    # Add feature loadings
    for (i in c("feature.loadings", "feature.loadings.projected")) {
      projected <- grepl(pattern = "projected", x = i)
      type <- ifelse(test = projected, yes = "projected loadings", no = "loadings")
      if (!IsMatrixEmpty(x = Loadings(object = x, projected = projected))) {
        if (verbose) {
          message("Adding ", type, " for ", name)
        }
        loadings <- Loadings(object = x, projected = projected)
        WriteH5Group(
          x = loadings, name = i, group = paste0(group, "/", name),
          hfile = hfile, verbose = verbose
        )
        WriteH5Group(
          x = rownames(x = loadings),
          name = ifelse(test = projected, yes = "projected.features", no = "features"),
          group = paste0(group, "/", name),
          hfile = hfile,
          verbose = verbose
        )
      } else if (verbose) {
        message("No ", type, " for ", name)
      }
    }
    # Add assay, key, and global status
    xgroup$create_attr(
      attr_name = "active.assay",
      robj = DefaultAssay(object = x),
      dtype = GuessDType(x = DefaultAssay(object = x))
    )
    xgroup$create_attr(
      attr_name = "key",
      robj = Key(object = x),
      dtype = GuessDType(x = Key(object = x))
    )
    xgroup$create_attr(
      attr_name = "global",
      robj = BoolToInt(x = IsGlobal(object = x)),
      dtype = GuessDType(x = IsGlobal(object = x))
    )
    # Add standard deviations
    if (length(x = Stdev(object = x)) > 0) {
      if (verbose) {
        message("Adding standard deviations for ", name)
      }
      WriteH5Group(
        x = Stdev(object = x),
        name = "stdev",
        group = paste0(group, "/", name),
        hfile = hfile,
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
        name = "jackstraw",
        group = paste0(group, "/", name),
        hfile = hfile,
        verbose = verbose
      )
    } else if (verbose) {
      message("No JackStraw data for ", name)
    }
    # Add misc
    WriteH5Group(
      x = slot(object = x, name = "misc"),
      name = "misc",
      group = paste0(group, "/", name),
      hfile = hfile,
      verbose = verbose
    )
    return(invisible(x = NULL))
  }
)

#' @rdname WriteH5Group
#'
#' @examples
#' \donttest{
#' # Factors turn into a group with two components: values and levels
#' # This is to preserve level information for HDF5 APIs that don't implement
#' # the HDF5 Enum type (eg. h5py)
#' # values corresponds to the integer values of each member of a factor
#' # levels is a string dataset with one entry per level
#' fctr <- factor(x = c("g1", "g1", "g2", "g1", "g2"))
#' WriteH5Group(x = fctr, name = "factor", hgroup = hfile)
#' hfile[["factor"]]
#' }
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "factor"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    if (group != "/") {
      hgroup <- hfile[[group]]
    } else {
      hgroup <- hfile
    }
    xgroup <- hgroup$create_group(name = name)
    # Write the integer values out
    WriteH5Group(
      x = as.integer(x = x),
      name = "values",
      group = paste0(group, "/", name),
      hfile = hfile,
      verbose = verbose
    )
    # Write the levels out
    WriteH5Group(
      x = levels(x = x),
      name = "levels",
      group = paste0(group, "/", name),
      hfile = hfile,
      verbose = verbose
    )
  }
)

#' @importClassesFrom Seurat Graph
#'
#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "Graph"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    hfile <- SparseWrite(x = x, name = name, group = group, hfile, verbose = verbose)
    if (group != "/") {
      hgroup <- hfile[[group]]
    } else {
      hgroup <- hfile
    }
    if (!is.null(x = DefaultAssay(object = x))) {
      hgroup[[name]]$create_attr(
        attr_name = "assay.used",
        robj = DefaultAssay(object = x),
        dtype = GuessDType(x = DefaultAssay(object = x))
      )
    }
    return(invisible(x = NULL))
  }
)

#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "list"),
  definition = BasicWrite
)

#' @rdname WriteH5Group
#'
#' @examples
#' \donttest{
#' # Logicals get encoded as integers with the following mapping
#' # FALSE becomes 0L
#' # TRUE becomes 1L
#' # NA becomes 2L
#' # These are stored as H5T_INTEGERS instead of H5T_LOGICALS
#' # Additionally, an attribute called "s3class" is written with the value of "logical"
#' WriteH5Group(c(TRUE, FALSE, NA), name = "logicals", hgroup = hfile)
#' hfile[["logicals"]]
#' hfile[["logicals"]]$attr_open("s3class")$read()
#' }
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "logical"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    WriteH5Group(
      x = BoolToInt(x = x),
      name = name,
      group = group,
      hfile = hfile,
      verbose = verbose
    )
    if (group != "/") {
      hgroup <- hfile[[group]]
    } else {
      hgroup <- hfile
    }
    hgroup[[name]]$create_attr(
      attr_name = "s3class",
      robj = "logical",
      dtype = GuessDType(x = "logical")
    )
    return(invisible(x = NULL))
  }
)

#' @importClassesFrom Seurat Neighbor
#'
#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "Neighbor"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    if (group != "/") {
      hgroup <- hfile[[group]]
    } else {
      hgroup <- hfile
    }
    xgroup <- hgroup$create_group(name = name)
    for (i in slotNames(x = x)) {
      if (i == "alg.idx" && !is.null(x = slot(object = x, name = i))) {
        warning(
          "We cannot save neighbor indexes at this time; ",
          "please save the index separately",
          call. = FALSE,
          immediate. = TRUE
        )
        next
      }
      WriteH5Group(
        x = slot(object = x, name = i),
        name = i,
        group = paste0(group, "/", name),
        hfile = hfile,
        verbose = verbose
      )
    }
  }
)

#' @importClassesFrom Seurat SeuratCommand
#'
#' @rdname WriteH5Group
#'
setMethod(
  f = "WriteH5Group",
  signature = c("x" = "SeuratCommand"),
  definition = function(x, name, group, hfile, verbose = TRUE) {
    # Write out params
    WriteH5Group(
      x = Filter(
        f = Negate(f = is.function),
        x = as.list(x = x)
      ),
      name = name,
      group = group,
      hfile = hfile,
      verbose = verbose
    )
    # Add extra information as HDF5 attributes
    for (slot in slotNames(x = x)) {
      if (slot == "params") {
        next
      }
      slot.val <- slot(object = x, name = slot)
      if (!is.null(x = slot.val)) {
        slot.val <- as.character(x = slot.val)
        if (group != "/") {
          hgroup <- hfile[[group]]
        } else {
          hgroup <- hfile
        }
        hgroup[[name]]$create_attr(
          attr_name = slot,
          robj = slot.val,
          dtype = GuessDType(x = slot.val)
        )
      }
    }
    return(invisible(x = NULL))
  }
)

#' @name WriteH5Group
#' @rdname WriteH5Group
#'
#' @examples
#' \donttest{
#' # Close and remove the HDF5 file
#' hfile$close_all()
#' file.remove(hfile$filename)
#' }
#'
NULL
