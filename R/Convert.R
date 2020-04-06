#' @include zzz.R
#' @include Connect.R
#' @include TestObject.R
#' @importFrom utils setTxtProgressBar
#' @importFrom hdf5r H5File h5attr H5S
#' @importFrom tools file_path_sans_ext
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Generics
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert an on-disk single-cell dataset to another format
#'
#' HDF5-based single-cell datasets can be converted from one format to another
#' using minimal memory. Details about conversion formats implemented are
#' provided below
#'
# @inheritSection H5ADToH5Seurat AnnData/H5AD to h5Seurat
#' @inheritSection H5SeuratToH5AD h5Seurat to AnnData/H5AD
#'
#' @param source Source dataset
#' @param dest Name of destination dataset
#' @param assay Converting from \code{\link{h5Seurat}}: name of assay to write
#' out; converting to \code{\link{h5Seurat}}: name to store assay data as
#' @param overwrite Overwrite existing \code{dest}
#' @param verbose Show progress updates
#' @param ... Arguments passed to other methods
#'
#' @return If \code{source} is a \code{character}, invisibly returns
#' \code{dest}; otherwise, returns an \code{\link[hdf5r]{H5File}}, or
#' filetype-specific subclass of \code{H5File} (eg. \code{\link{h5Seurat}}),
#' connection to \code{dest}
#'
#' @name Convert
#' @rdname Convert
#'
#' @export
#'
Convert <- function(source, dest, assay, overwrite = FALSE, verbose = TRUE, ...) {
  if (!missing(x = dest) && !is.character(x = dest)) {
    stop("'dest' must be a filename or type", call. = FALSE)
  }
  UseMethod(generic = 'Convert', object = source)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Methods
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' @importFrom hdf5r H5File
#'
#' @rdname Convert
#' @method Convert character
#' @export
#'
Convert.character <- function(
  source,
  dest,
  assay,
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  hfile <- Connect(filename = source, force = TRUE)
  on.exit(expr = hfile$close_all())
  dfile <- Convert(
    source = hfile,
    dest = dest,
    assay = assay,
    overwrite = overwrite,
    verbose = verbose,
    ...
  )
  dfile$close_all()
  return(invisible(x = dfile$filename))
}

#' @rdname Convert
#' @method Convert H5File
#' @export
#'
Convert.H5File <- function(
  source,
  dest,
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  .NotYetImplemented()
  stype <- FileType(file = source$filename)
  dtype <- FileType(file = dest)
  if (tolower(x = dest) == dtype) {
    dest <- paste(file_path_sans_ext(x = source$filename), dtype, sep = '.')
  }
  dfile <- switch(
    EXPR = stype,
    'h5ad' = switch(
      EXPR = dtype,
      'h5seurat' = '',
      stop("Unable to convert H5AD files to ", dtype, " files", call. = FALSE)
    ),
    stop("Unknown file type: ", stype, call. = FALSE)
  )
  return(dfile)
}

#' @rdname Convert
#' @method Convert h5Seurat
#' @export
#'
Convert.h5Seurat <- function(
  source,
  dest = 'h5ad',
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  type <- FileType(file = dest)
  if (tolower(x = dest) == type) {
    dest <- paste(file_path_sans_ext(x = source$filename), type, sep = '.')
  }
  dfile <- switch(
    EXPR = type,
    'h5ad' = H5SeuratToH5AD(
      source = source,
      dest = dest,
      assay = assay,
      overwrite = overwrite,
      verbose = verbose
    ),
    stop("Unable to convert h5Seurat files to ", type, " files", call. = FALSE)
  )
  return(dfile)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Implementations
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Convert AnnData/H5AD files to h5Seurat files
#'
#' @inheritParams Convert
#'
#' @return Returns a handle to \code{dest} as an \code{\link{h5Seurat}} object
#'
#' @section AnnData/H5AD to h5Seurat:
#' The AnnData/H5AD to h5Seurat conversion will try to automatically fill in
#' datasets based on data presence. It works in the following manner:
#' \subsection{Expression data:}{
#'  The expression matrices \code{counts}, \code{data}, and \code{scale.data}
#'  are filled by \code{/X} and \code{/raw/X} in the following manner:
#'  \itemize{
#'   \item \code{counts} will be filled with \code{/raw/X} if present;
#'   otherwise, it will be filled with \code{/X}
#'   \item \code{data} will be filled with \code{/raw/X} if \code{/raw/X} is
#'   present and \code{/X} is dense; otherwise, it will be filled with \code{/X}
#'   \item \code{scale.data} will be filled with \code{/X} if it dense;
#'   otherwise, it will be empty
#'  }
#'  TODO: note about meta.features
#' }
#' \subsection{Cell-level metadata:}{
#' Cell-level metadata is added to \code{meta.data}; the row names of the
#' metadata (as determined by the value of the "_index" attribute or the "index"
#' dataset, in that order) are added to the "cell.names" dataset instead. If the
#' \code{__categories} dataset is present, each dataset within
#' \code{__categories} will be stored as a factor group. Cell-level metadata
#' will be added as an HDF5 group unless factors are \strong{not} present and
#' \code{\link[SeuratDisk]{SeuratDisk.dtype.dataframe_as_group}} is \code{FALSE}
#' }
#' \subsection{Dimensional reduction information:}{
#'  TODO: this
#' }
#' \subsection{Nearest-neighbor graph:}{
#' If a nearest neighbor graph is present in \code{/uns/neighbors/distances}, it
#' will be added as a graph dataset in the h5Seurat file and associated with
#' \code{assay}; if a value is present in \code{/uns/neighbors/params/metric},
#' the name of the graph will be \code{assay_metric}, otherwise, it will be
#' \code{assay_anndata}}
#'
#' @keywords internal
#'
H5ADToH5Seurat <- function(
  source,
  dest,
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE
) {
  .NotYetImplemented()
  if (file.exists(dest)) {
    if (overwrite) {
      file.remove(dest)
    } else {
      stop("Destination h5Seurat file exists", call. = FALSE)
    }
  }
  dfile <- h5Seurat$new(filename = dest, mode = WriteMode(overwrite = FALSE))
  GetRownames <- function(dset) {
    if (!IsDataFrame(x = source[[dset]])) {
      stop("'", dset, "' is not a data frame", call. = FALSE)
    }
    if (inherits(x = source[[dset]], what = 'H5Group')) {
      rownames <- if (source[[dset]]$attr_exists(attr_name = '_index')) {
        h5attr(x = source[[dset]], which = '_index')
      } else if (source[[dset]]$exists(name = '_index')) {
        '_index'
      } else if (source[[dset]]$exists(name = 'index')) {
        'index'
      } else {
        stop("Cannot find rownames in ", dset, call. = FALSE)
      }
    } else {
      # TODO: fix this
      stop("Don't know how to handle datasets", call. = FALSE)
    }
    return(rownames)
  }
  if (source$exists(name = 'raw')) {
    shape.var <- 'raw/var'
    shape.x <- 'raw/X'
  } else {
    shape.var <- 'var'
    shape.x <- 'X'
  }
  ds.map <- c(
    scale.data = if (inherits(x = source[['X']], what = 'H5D')) {
      'X'
    } else {
      NULL
    },
    data = if (inherits(x = source[['X']], what = 'H5D') && source$exists(name = 'raw')) {
      'raw/X'
    } else {
      'X'
    },
    counts = if (source$exists(name = 'raw')) {
      'raw/X'
    } else {
      'X'
    }
  )
  # Add assay data
  assay.group <- dfile[['assays']]$create_group(name = assay)
  for (i in seq_along(along.with = ds.map)) {
    assay.group$obj_copy_from(
      src_loc = source,
      src_name = ds.map[[i]],
      dst_name = names(x = ds.map)[i]
    )
  }
  if (source$exists(name = 'raw/var')) {
    ''
  }
  # TODO: add meta.features
  # Add cell-level metadata
  # Add dimensional reduction information
  # Add nearest-neighbor graph
  return(dfile)
}

#' Convert h5Seurat files to H5AD files
#'
#' @inheritParams Convert
#'
#' @return Returns a handle to \code{dest} as an \code{\link[hdf5r]{H5File}}
#' object
#'
#' @section h5Seurat to AnnData/H5AD:
#' The h5Seurat to AnnData/H5AD conversion will try to automatically fill in
#' datasets based on data presence. Data presense is determined by the h5Seurat
#' index (\code{source$index()}). It works in the following manner:
#' \subsection{Assay data:}{
#'  \itemize{
#'   \item \code{X} will be filled with \code{scale.data} if \code{scale.data}
#'   is present; otherwise, it will be filled with \code{data}
#'   \item \code{var} will be filled with \code{meta.features} \strong{only} for
#'   the features present in \code{X}; for example, if \code{X} is filled with
#'   \code{scale.data}, then \code{var} will contain only features that have
#'   been scaled
#'   \item \code{raw.X} will be filled with \code{data} if \code{X} is filled
#'   with \code{scale.data}; otherwise, it will be filled with \code{counts}. If
#'   \code{counts} is not present, then \code{raw} will not be filled
#'   \item \code{raw.var} will be filled with \code{meta.features} with the
#'   features present in \code{raw.X}; if \code{raw.X} is not filled, then
#'   \code{raw.var} will not be filled
#'  }
#' }
#' \subsection{Cell-level metadata:}{
#'  Cell-level metadata is added to \code{obs}
#' }
#' \subsection{Dimensional reduction information}{
#'  Only dimensional reductions associated with \code{assay} or marked as
#'  \link[Seurat:IsGlobal]{global} will be transfered to the H5AD file. For
#'  every reduction \code{reduc}:
#'  \itemize{
#'   \item cell embeddings are placed in \code{obsm} and renamed to
#'   \code{X_reduc}
#'   \item feature loadings, if present, are placed in \code{varm} and renamed
#'   to either \code{PCs} if \code{reduc} is "pca" otherwise \code{reduc} in all
#'    caps
#'  }
#'  For example, if \code{reduc} is "ica", then cell embeddings will be "X_ica"
#'  in \code{obsm} and feature loaodings, if present, will be "ICA" in
#'  \code{varm}
#' }
#' \subsection{Nearest-neighbor graphs}{
#'  If a nearest-neighbor graph is associated with \code{assay}, it will be
#'  added to \code{uns/neighbors/distances}; if more than one graph is present,
#'  then \strong{only} the last graph according to the index will be added.
#' }
#'
#' @keywords internal
#'
H5SeuratToH5AD <- function(
  source,
  dest,
  assay = DefaultAssay(object = source),
  overwrite = FALSE,
  verbose = TRUE
) {
  if (file.exists(dest)) {
    if (overwrite) {
      file.remove(dest)
    } else {
      stop("Destination H5AD file exists", call. = FALSE)
    }
  }
  rownames <- '_index'
  dfile <- H5File$new(filename = dest, mode = WriteMode(overwrite = FALSE))
  # Transfer data frames from h5Seurat files to H5AD files
  #
  # @param src Source dataset
  # @param dname Name of destination
  # @param index Integer values of rows to take
  #
  # @return Invisibly returns \code{NULL}
  #
  TransferDF <- function(src, dname, index) {
    if (verbose) {
      message("Transfering ", basename(path = src$get_obj_name()), " to ", dname)
    }
    if (inherits(x = src, what = 'H5D')) {
      CompoundToGroup(
        src = src,
        dest = dfile,
        dst.name = dname,
        order = 'column-order',
        index = index
      )
    } else {
      dfile$create_group(name = dname)
      for (i in src$names) {
        if (IsFactor(x = src[[i]])) {
          dfile[[dname]]$create_dataset(
            name = i,
            robj = src[[i]][['values']][index],
            dtype = GuessDType(x = 1L)
          )
          if (!dfile[[dname]]$exists(name = '__categories')) {
            dfile[[dname]]$create_group(name = '__categories')
          }
          dfile[[dname]][['__categories']]$create_dataset(
            name = i,
            robj = src[[i]][['levels']][],
            GuessDType(x = src[[i]][['levels']][])
          )
        } else {
          dfile[[dname]]$create_dataset(
            name = i,
            robj = src[[i]][index],
            dtype = GuessDType(x = src[[i]][1])
          )
        }
      }
      if (src$attr_exists(attr_name = 'colnames')) {
        dfile[[dname]]$create_attr(
          attr_name = 'column-order',
          robj = h5attr(x = src, which = 'colnames'),
          dtype = GuessDType(x = h5attr(x = src, which = 'colnames'))
        )
      }
    }
    return(invisible(x = NULL))
  }
  # @rdname TransferDF
  #
  TransferMatrix <- function(src, dname) {
    if (verbose) {
      pb <- PB()
    }
    mtx.dims <- rev(x = src$dims)
    dset <- dfile$create_dataset(
      name = dname,
      dtype = src$get_type(),
      space = H5S$new(dims = mtx.dims, maxdims = mtx.dims)
    )
    chunk.points <- source$chunk.points(
      dataset = src$get_obj_name(),
      MARGIN = 1L
    )
    for (i in 1:nrow(x = chunk.points)) {
      select <- seq.default(
        from = chunk.points[i, 'start'],
        to = chunk.points[i, 'end']
      )
      dset[, select] <- src[select, ]
      if (verbose) {
        setTxtProgressBar(pb = pb, value = i / nrow(x = chunk.points))
      }
    }
    if (verbose) {
      close(con = pb)
    }
    return(invisible(x = NULL))
  }
  # Add assay data
  assay.group <- source[['assays']][[assay]]
  if (source$index()[[assay]]$slots[['scale.data']]) {
    x.data <- 'scale.data'
    raw.data <- 'data'
  } else {
    x.data <- 'data'
    raw.data <- if (source$index()[[assay]]$slots[['counts']]) {
      'counts'
    } else {
      NULL
    }
  }
  if (verbose) {
    message("Adding ", x.data, " from ", assay, " as X")
  }
  assay.group$obj_copy_to(dst_loc = dfile, dst_name = 'X', src_name = x.data)
  if (dfile[['X']]$attr_exists(attr_name = 'dims')) {
    dims <- h5attr(x = dfile[['X']], which = 'dims')
    dfile[['X']]$create_attr(
      attr_name = 'shape',
      robj = rev(x = dims),
      dtype = GuessDType(x = dims)
    )
    dfile$attr_delete(attr_name = 'dims')
  }
  x.features <- switch(
    EXPR = x.data,
    'scale.data' = which(x = assay.group[['features']][] %in% assay.group[['scaled.features']][]),
    seq.default(from = 1, to = assay.group[['features']]$dims)
  )
  # Add meta.features
  if (assay.group$exists(name = 'meta.features')) {
    TransferDF(
      src = assay.group[['meta.features']],
      dname = 'var',
      index = x.features
    )
  } else {
    dfile$create_group(name = 'var')
  }
  # Add feature names
  dfile[['var']]$create_dataset(
    name = rownames,
    robj = assay.group[['features']][x.features],
    dtype = GuessDType(x = assay.group[['features']][1])
  )
  dfile[['var']]$create_attr(
    attr_name = rownames,
    robj = rownames,
    dtype = GuessDType(x = rownames)
  )
  # Add raw
  if (!is.null(x = raw.data)) {
    if (verbose) {
      message("Adding ", raw.data, " from ", assay, " as raw")
    }
    dfile$create_group(name = 'raw')
    assay.group$obj_copy_to(
      dst_loc = dfile[['raw']],
      dst_name = 'X',
      src_name = raw.data
    )
    if (dfile[['raw/X']]$attr_exists(attr_name = 'dims')) {
      dims <- h5attr(x = dfile[['raw/X']], which = 'dims')
      dfile[['raw/X']]$create_attr(
        attr_name = 'shape',
        robj = rev(x = dims),
        dtype = GuessDType(x = dims)
      )
      dfile[['raw/X']]$attr_delete(attr_name = 'dims')
    }
    # Add meta.features
    if (assay.group$exists(name = 'meta.features')) {
      TransferDF(
        src = assay.group[['meta.features']],
        dname = 'raw/var',
        index = seq.default(from = 1, to = assay.group[['features']]$dims)
      )
    } else {
      dfile[['raw']]$create_group(name = 'var')
    }
    # Add feature names
    dfile[['raw/var']]$create_dataset(
      name = rownames,
      robj = assay.group[['features']][],
      dtype = GuessDType(x = assay.group[['features']][1])
    )
    dfile[['raw/var']]$create_attr(
      attr_name = rownames,
      robj = rownames,
      dtype = GuessDType(x = rownames)
    )
  }
  # Add cell-level metadata
  TransferDF(
    src = source[['meta.data']],
    dname = 'obs',
    index = seq.default(from = 1, to = length(x = Cells(x = source)))
  )
  dfile[['obs']]$create_dataset(
    name = rownames,
    robj = Cells(x = source),
    dtype = GuessDType(x = Cells(x = source))
  )
  dfile[['obs']]$create_attr(
    attr_name = rownames,
    robj = rownames,
    dtype = GuessDType(x = rownames)
  )
  # Add dimensional reduction information
  dfile$create_group(name = 'obsm')
  dfile$create_group(name = 'varm')
  reductions <- source$index()[[assay]]$reductions
  for (reduc in names(x = reductions)) {
    if (verbose) {
      message("Adding dimensional reduction information for ", reduc)
    }
    TransferMatrix(
      src = source[['reductions']][[reduc]][['cell.embeddings']],
      dname = paste0('obsm/X_', reduc)
    )
    if (reductions[[reduc]]['feature.loadings']) {
      if (verbose) {
        message("Adding feature loadings for ", reduc)
      }
      TransferMatrix(
        src = source[['reductions']][[reduc]][['feature.loadings']],
        dname = paste0(
          'varm/',
          switch(EXPR = reduc, 'pca' = 'PCs', toupper(x = reduc))
        )
      )
    }
  }
  # Add global dimensional reduction information
  global.reduc <- source$index()[['globals']][['reductions']]
  for (reduc in global.reduc) {
    if (reduc %in% names(x = reductions)) {
      next
    } else if (verbose) {
      message("Adding dimensional reduction information for ", reduc, " (global)")
    }
    TransferMatrix(
      src = source[['reductions']][[reduc]][['cell.embeddings']],
      dname = paste0('obsm/X_', reduc)
    )
  }
  # Create uns
  dfile$create_group(name = 'uns')
  # Add graph
  graph <- source$index()[[assay]]$graphs
  graph <- graph[length(x = graph)]
  if (!is.null(x = graph)) {
    if (verbose) {
      message("Adding ", graph, " as neighbors")
    }
    dgraph <- dfile[['uns']]$create_group(name = 'neighbors')
    source[['graphs']]$obj_copy_to(
      dst_loc = dgraph,
      dst_name = 'distances',
      src_name = graph
    )
    if (source[['graphs']][[graph]]$attr_exists(attr_name = 'dims')) {
      dims <- h5attr(x = source[['graphs']][[graph]], which = 'dims')
      dgraph[['distances']]$create_attr(
        attr_name = 'shape',
        robj = rev(x = dims),
        dtype = GuessDType(x = dims)
      )
    }
    dgraph$create_group(name = 'params')
    dgraph[['params']]$create_dataset(
      name = 'metric',
      robj = gsub(pattern = paste0('^', assay, '_'), replacement = '', x = graph),
      dtype = GuessDType(x = graph)
    )
  }
  return(dfile)
}
