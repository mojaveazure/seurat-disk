#' @include zzz.R
#' @include Connect.R
#' @include TestObject.R
#' @include Transpose.R
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
#' @inheritSection H5ADToH5Seurat AnnData/H5AD to h5Seurat
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
  dest = 'h5seurat',
  assay = 'RNA',
  overwrite = FALSE,
  verbose = TRUE,
  ...
) {
  stype <- FileType(file = source$filename)
  dtype <- FileType(file = dest)
  if (tolower(x = dest) == dtype) {
    dest <- paste(file_path_sans_ext(x = source$filename), dtype, sep = '.')
  }
  dfile <- switch(
    EXPR = stype,
    'h5ad' = switch(
      EXPR = dtype,
      'h5seurat' = H5ADToH5Seurat(
        source = source,
        dest = dest,
        assay = assay,
        overwrite = overwrite,
        verbose = verbose
      ),
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
#'  Feature names are taken from the feature-level metadata
#' }
#' \subsection{Feature-level metadata:}{
#' Feature-level metadata is added to the \code{meta.features} datasets in each
#' assay. Feature names are taken from the dataset specified by the
#' \dQuote{_index} attribute, the \dQuote{_index} dataset, or the \dQuote{index}
#' dataset, in that order. Metadata is populated with \code{/raw/var} if
#' present, otherwise with \code{/var}; if both \code{/raw/var} and \code{/var}
#' are present, then \code{meta.features} will be populated with \code{/raw/var}
#' first, then \code{/var} will be added to it. For columns present in both
#' \code{/raw/var} and \code{/var}, the values in \code{/var} will be used
#' instead. \strong{Note}: it is possible for \code{/var} to have fewer features
#' than \code{/raw/var}; if this is the case, then only the features present in
#' \code{/var} will be overwritten, with the metadata for features \emph{not}
#' present in \code{/var} remaining as they were in \code{/raw/var} or empty
#' }
#' \subsection{Cell-level metadata:}{
#' Cell-level metadata is added to \code{meta.data}; the row names of the
#' metadata (as determined by the value of the \dQuote{_index} attribute, the
#' \dQuote{_index} dataset, or the \dQuote{index} dataset, in that order) are
#' added to the \dQuote{cell.names} dataset instead. If the
#' \dQuote{__categories} dataset is present, each dataset within
#' \dQuote{__categories} will be stored as a factor group. Cell-level metadata
#' will be added as an HDF5 group unless factors are \strong{not} present and
#' \code{\link[SeuratDisk]{SeuratDisk.dtype.dataframe_as_group}} is \code{FALSE}
#' }
#' \subsection{Dimensional reduction information:}{
#'  Cell embeddings are taken from \code{/obsm}; dimensional reductions are
#'  named based on their names from \code{obsm} by removing the preceding
#'  \dQuote{X_}.For example, if a dimensional reduction is named \dQuote{X_pca}
#'  in \code{/obsm}, the resulting dimensional reduction information will be
#'  named \dQuote{pca}. The key will be set to one of the following:
#'  \itemize{
#'   \item \dQuote{PC_} if \dQuote{pca} is present in the dimensional reduction
#'   name (\code{grepl("pca", reduction.name, ignore.case = TRUE)})
#'   \item \dQuote{tSNE_} if \dQuote{tsne} is present in the dimensional
#'   reduction name (\code{grepl("tsne", reduction.name, ignore.case = TRUE)})
#'   \item \code{reduction.name_} for all other reductions
#'  }
#'  Remember that the preceding \dQuote{X_} will be removed from the reduction
#'  name before converting to a key. Feature loadings are taken from
#'  \code{/varm} and placed in the associated dimensional reduction. The
#'  dimensional reduction is determine from the loadings name in \code{/varm}:
#'  \itemize{
#'   \item \dQuote{PCs} will be added to a dimensional reduction named
#'   \dQuote{pca}
#'   \item All other loadings in \code{/varm} will be added to a dimensional
#'   reduction named \code{tolower(loading)} (eg. a loading named \dQuote{ICA}
#'   will be added to a dimensional reduction named \dQuote{ica})
#'  }
#'  If a dimensional reduction cannot be found according to the rules above, the
#'  loading will not be taken from the AnnData/H5AD file. Miscellaneous
#'  information will be taken from \code{/uns/reduction} where \code{reduction}
#'  is the name of the reduction in \code{/obsm} without the preceding
#'  \dQuote{X_}; if no dimensional reduction information present, then
#'  miscellaneous information will not be taken from the AnnData/H5AD file.
#'  Standard deviations are taken from a dataset \code{/uns/reduction/variance};
#'  the variances will be converted to standard deviations and added to the
#'  \code{stdev} dataset of a dimensional reduction
#' }
#' \subsection{Nearest-neighbor graph:}{
#' If a nearest neighbor graph is present in \code{/uns/neighbors/distances}, it
#' will be added as a graph dataset in the h5Seurat file and associated with
#' \code{assay}; if a value is present in \code{/uns/neighbors/params/method},
#' the name of the graph will be \code{assay_method}, otherwise, it will be
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
  if (file.exists(dest)) {
    if (overwrite) {
      file.remove(dest)
    } else {
      stop("Destination h5Seurat file exists", call. = FALSE)
    }
  }
  dfile <- h5Seurat$new(filename = dest, mode = WriteMode(overwrite = FALSE))
  # Get rownames from an H5AD data frame
  #
  # @param dset Name of data frame
  #
  # @return Returns the name of the dataset that contains the rownames
  #
  GetRownames <- function(dset) {
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
  ColToFactor <- function(dfgroup) {
    if (dfgroup$exists(name = '__categories')) {
      for (i in names(x = dfgroup[['__categories']])) {
        tname <- basename(path = tempfile(tmpdir = ''))
        dfgroup$obj_copy_to(dst_loc = dfgroup, dst_name = tname, src_name = i)
        dfgroup$link_delete(name = i)
        dfgroup$create_group(name = i)
        dfgroup$obj_copy_to(
          dst_loc = dfgroup,
          dst_name = paste0(i, '/values'),
          src_name = tname
        )
        dfgroup$obj_copy_to(
          dst_loc = dfgroup,
          dst_name = paste0(i, '/levels'),
          src_name = paste0('__categories/', i)
        )
        dfgroup$link_delete(name = tname)
      }
      dfgroup$link_delete(name = '__categories')
    }
    return(invisible(x = NULL))
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
    if (verbose) {
      message("Adding ", ds.map[[i]], " as ", names(x = ds.map)[i])
    }
    dst <- names(x = ds.map)[i]
    assay.group$obj_copy_from(
      src_loc = source,
      src_name = ds.map[[i]],
      dst_name = dst
    )
    if (assay.group[[dst]]$attr_exists(attr_name = 'shape')) {
      dims <- rev(x = h5attr(x = assay.group[[dst]], which = 'shape'))
      assay.group[[dst]]$create_attr(
        attr_name = 'dims',
        robj = dims,
        dtype = GuessDType(x = dims)
      )
    }
  }
  features.source <- ifelse(
    test = source$exists(name = 'raw/var'),
    yes = 'raw/var',
    no = 'var'
  )
  features.dset <- GetRownames(dset = features.source)
  assay.group$obj_copy_from(
    src_loc = source,
    src_name = paste(features.source, features.dset, sep = '/'),
    dst_name = 'features'
  )
  scaled <- !is.null(x = ds.map['scale.data']) && !is.null(x = ds.map['scale.data'])
  if (scaled) {
    scaled.dset <- GetRownames(dset = 'var')
    assay.group$obj_copy_from(
      src_loc = source,
      src_name = paste0('var/', scaled.dset),
      dst_name = 'scaled.features'
    )
  }
  # Add feature-level metadata
  if (!getOption(x = "Seuratdisk.dtypes.dataframe_as_group", default = FALSE)) {
    warning(
      "Adding feature-level metadata as a compound is not yet supported",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (source$exists(name = 'raw/var')) {
    if (verbose) {
      message("Adding meta.features from raw/var")
    }
    assay.group$obj_copy_from(
      src_loc = source,
      src_name = 'raw/var',
      dst_name = 'meta.features'
    )
    if (scaled) {
      features.use <- assay.group[['features']][] %in% assay.group[['scaled.features']][]
      features.use <- which(x = features.use)
      meta.scaled <- names(x = source[['var']])
      meta.scaled <- meta.scaled[!meta.scaled %in% c('__categories', scaled.dset)]
      for (mf in meta.scaled) {
        if (!mf %in% names(x = assay.group[['meta.features']])) {
          if (verbose) {
            message("Adding ", mf, " from scaled feature-level metadata")
          }
          assay.group[['meta.features']]$create_dataset(
            name = mf,
            dtype = source[['var']][[mf]]$get_type(),
            space = H5S$new(dims = assay.group[['features']]$dims)
          )
        } else if (verbose) {
          message("Merging ", mf, " from scaled feature-level metadata")
        }
        assay.group[['meta.features']][[mf]][features.use] <- source[['var']][[mf]]$read()
      }
    }
  } else {
    if (verbose) {
      message("Adding meta.features from var")
    }
    assay.group$obj_copy_from(
      src_loc = source,
      src_name = 'var',
      dst_name = 'meta.features'
    )
  }
  ColToFactor(dfgroup = assay.group[['meta.features']])
  if (assay.group[['meta.features']]$attr_exists(attr_name = 'column-order')) {
    colnames <- h5attr(
      x = assay.group[['meta.features']],
      which = 'column-order'
    )
    assay.group[['meta.features']]$create_attr(
      attr_name = 'colnames',
      robj = colnames,
      dtype = GuessDType(x = colnames)
    )
  }
  assay.group[['meta.features']]$link_delete(name = GetRownames(dset = 'var'))
  # Add cell-level metadata
  if (source$exists(name = 'obs')) {
    if (!source[['obs']]$exists(name = '__categories') && !getOption("SeuratDisk.dtypes.dataframe_as_group", x = TRUE)) {
      warning(
        "Conversion from H5AD to h5Seurat allowing compound datasets is not yet implemented",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    dfile$obj_copy_from(
      src_loc = source,
      src_name = 'obs',
      dst_name = 'meta.data'
    )
    ColToFactor(dfgroup = dfile[['meta.data']])
    if (dfile[['meta.data']]$attr_exists(attr_name = 'column-order')) {
      colnames <- h5attr(x = dfile[['meta.data']], which = 'column-order')
      dfile[['meta.data']]$create_attr(
        attr_name = 'colnames',
        robj = colnames,
        dtype = GuessDType(x = colnames)
      )
    }
    rownames <- GetRownames(dset = 'obs')
    dfile$obj_copy_from(
      src_loc = dfile,
      src_name = paste0('meta.data/', rownames),
      dst_name = 'cell.names'
    )
    dfile[['meta.data']]$link_delete(name = rownames)
  } else {
    warning(
      "No cell-level metadata present, creating fake cell names",
      call. = FALSE,
      immediate. = TRUE
    )
    ncells <- if (inherits(x = assay.group[['data']], what = 'H5Group')) {
      assay.group[['data/indptr']]$dims - 1
    } else {
      assay.group[['data']]$dims[2]
    }
    dfile$create_group(name = 'meta.data')
    dfile$create_dataset(
      name = 'cell.names',
      robj = paste0('Cell', seq.default(from = 1, to = ncells)),
      dtype = GuessDType(x = 'Cell1')
    )
  }
  # Add dimensional reduction information
  if (source$exists(name = 'obsm')) {
    # Add cell embeddings
    if (inherits(x = source[['obsm']], what = 'H5Group')) {
      for (reduc in names(x = source[['obsm']])) {
        sreduc <- gsub(pattern = '^X_', replacement = '', x = reduc)
        reduc.group <- dfile[['reductions']]$create_group(name = sreduc)
        message("Adding ", reduc, " as cell embeddings for ", sreduc)
        Transpose(
          x = source[['obsm']][[reduc]],
          dest = reduc.group,
          dname = 'cell.embeddings',
          verbose = verbose
        )
        reduc.group$create_group(name = 'misc')
        reduc.group$create_attr(
          attr_name = 'active.assay',
          robj = assay,
          dtype = GuessDType(x = assay)
        )
        key <- paste0(
          if (grepl(pattern = 'pca', x = sreduc, ignore.case = TRUE)) {
            'PC'
          } else if (grepl(pattern = 'tsne', x = sreduc, ignore.case = TRUE)) {
            'tSNE'
          } else {
            sreduc
          },
          '_'
        )
        reduc.group$create_attr(
          attr_name = 'key',
          robj = key,
          dtype = GuessDType(x = reduc)
        )
        global <- BoolToInt(x = grepl(
          pattern = 'tsne|umap',
          x = sreduc,
          ignore.case = TRUE
        ))
        reduc.group$create_attr(
          attr_name = 'global',
          robj = global,
          dtype = GuessDType(x = global)
        )
      }
    } else {
      warning(
        "Reading compound dimensional reductions not yet supported, please update your H5AD file",
        call. = FALSE,
        immediate. = TRUE
      )
    }
    # Add feature loadings
    if (source$exists(name = 'varm')) {
      if (inherits(x = source[['varm']], what = 'H5Group')) {
        for (reduc in names(x = source[['varm']])) {
          sreduc <- switch(EXPR = reduc, 'PCs' = 'pca', tolower(x = reduc))
          if (!isTRUE(x = sreduc %in% names(x = dfile[['reductions']]))) {
            warning(
              "Cannot find a reduction named ",
              sreduc,
              " (",
              reduc,
              " in varm)",
              call. = FALSE,
              immediate. = TRUE
            )
            next
          }
          if (isTRUE(x = verbose)) {
            message("Adding ", reduc, " as feature loadings fpr ", sreduc)
          }
          Transpose(
            x = source[['varm']][[reduc]],
            dest = dfile[['reductions']][[sreduc]],
            dname = 'feature.loadings',
            verbose = verbose
          )
        }
      } else {
        warning(
          "Reading compound dimensional reductions not yet supported",
          call. = FALSE,
          immediate. = TRUE
        )
      }
    }
    # Add miscellaneous information
    if (source$exists(name = 'uns')) {
      for (reduc in names(x = source[['uns']])) {
        if (!isTRUE(x = reduc %in% names(x = dfile[['reductions']]))) {
          next
        }
        if (verbose) {
          message("Adding miscellaneous information for ", reduc)
        }
        dfile[['reductions']][[reduc]]$link_delete(name = 'misc')
        dfile[['reductions']][[reduc]]$obj_copy_from(
          src_loc = source[['uns']],
          src_name = reduc,
          dst_name = 'misc'
        )
        if ('variance' %in% names(x = dfile[['reductions']][[reduc]][['misc']])) {
          if (verbose) {
            message("Adding standard deviations for ", reduc)
          }
          dfile[['reductions']][[reduc]]$create_dataset(
            name = 'stdev',
            robj = sqrt(x = dfile[['reductions']][[reduc]][['misc']][['variance']][]),
            dtype = GuessDType(x = 1.0)
          )
        }
      }
    }
  }
  # Add nearest-neighbor graph
  if (source$exists('uns/neighbors/distances')) {
    graph.name <- paste(
      assay,
      ifelse(
        test = source$exists(name = 'uns/neighbors/params/method'),
        yes = source[['uns/neighbors/params/method']][1],
        no = 'anndata'
      ),
      sep = '_'
    )
    if (verbose) {
      message("Saving nearest-neighbor graph as ", graph.name)
    }
    dfile[['graphs']]$obj_copy_from(
      src_loc = source,
      src_name = 'uns/neighbors/distances',
      dst_name = graph.name
    )
    if (dfile[['graphs']][[graph.name]]$attr_exists(attr_name = 'shape')) {
      dfile[['graphs']][[graph.name]]$create_attr(
        attr_name = 'dims',
        robj = h5attr(x = dfile[['graphs']][[graph.name]], which = 'shape'),
        dtype = GuessDType(x = h5attr(
          x = dfile[['graphs']][[graph.name]],
          which = 'shape'
        ))
      )
    }
    dfile[['graphs']][[graph.name]]$create_attr(
      attr_name = 'assay.used',
      robj = assay,
      dtype = GuessDType(x = assay)
    )
  }
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
#'   to either \dQuote{PCs} if \code{reduc} is \dQuote{pca} otherwise
#'   \code{reduc} in all caps
#'  }
#'  For example, if \code{reduc} is \dQuote{ica}, then cell embeddings will be
#'  \dQuote{X_ica} in \code{obsm} and feature loaodings, if present, will be
#'  \dQuote{ICA} in \code{varm}
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
      name = 'method',
      robj = gsub(pattern = paste0('^', assay, '_'), replacement = '', x = graph),
      dtype = GuessDType(x = graph)
    )
  }
  return(dfile)
}
