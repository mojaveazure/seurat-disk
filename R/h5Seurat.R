#' @include scdisk.R
#'
NULL

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Class definition
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' A class for connections to h5Seurat files
#'
#' @docType class
#' @name h5Seurat-class
#' @rdname h5Seurat-class
#' @aliases h5Seurat
#' @format An \code{\link[R6]{R6Class}} object
#' @seealso \code{\link[hdf5r]{H5File}}
#'
#' @importFrom R6 R6Class
#' @importFrom hdf5r H5File h5attr
#' @importFrom utils packageVersion
#'
#' @export
#'
h5Seurat <- R6Class(
  classname = 'h5Seurat',
  inherit = scdisk,
  cloneable = FALSE,
  portable = TRUE,
  lock_class = TRUE,
  public = list(
    # Methods
    #' @description Get the index for this h5Seurat file
    index = function() {
      if (!length(x = private$index.internal)) {
        private$build.index(
          version = ClosestVersion(
            query = self$version(),
            targets = private$versions
          )
        )
      }
      return(private$index.internal)
    },
    #' @description Set the version attribute
    #' @param version A version number matching the regex
    #' \code{^\\d+(\\.\\d+){2}(\\.9\\d{3})?$}
    set.version = function(version) {
      version <- as.character(x = version)
      if (!grepl(pattern = version.regex, x = version)) {
        stop("Invalid version specification: ", version, call. = FALSE)
      }
      if (self$attr_exists(attr_name = 'version')) {
        self$attr_delete(attr_name = 'version')
      }
      self$create_attr(
        attr_name = 'version',
        robj = version,
        dtype = GuessDType(x = version)
      )
      if (!self$mode %in% modes$new) {
        private$validate()
      }
      return(invisible(x = self))
    },
    #' @description Get the version attribute
    version = function() {
      return(h5attr(x = self, which = 'version'))
    }
  ),
  private = list(
    # Fields
    index.internal = list(),
    versions = c('3.1.2', '3.1.4.9900'),
    # Methods
    build.index = function(version) {
      version <- match.arg(arg = version, choices = private$versions)
      version <- numeric_version(x = version)
      # Get Assay information
      index <- sapply(
        X = names(x = self[['assays']]),
        FUN = function(x) {
          slots <- c('counts', 'data', 'scale.data')
          check <- slots %in% names(x = self[['assays']][[x]])
          names(x = check) <- slots
          check[['scale.data']] <- check[['scale.data']] && self[['assays']][[x]]$exists(name = 'scaled.features')
          check <- list(check)
          names(x = check) <- 'slots'
          return(check)
        },
        simplify = FALSE,
        USE.NAMES = TRUE
      )
      # Get DimReduc information
      reduc.slots <- c(
        'cell.embeddings',
        'feature.loadings',
        'feature.loadings.projected',
        'jackstraw'
      )
      for (reduc in names(x = self[['reductions']])) {
        reduc.assay <- h5attr(
          x = self[['reductions']][[reduc]],
          which = 'active.assay'
        )
        if (!reduc.assay %in% names(x = index)) {
          warning(
            "Cannot find assay ",
            reduc.assay,
            " in the H5Seurat file",
            call. = FALSE,
            immediate. = TRUE
          )
          next
        }
        check <- reduc.slots %in% names(x = self[['reductions']][[reduc]])
        names(x = check) <- reduc.slots
        if (check[['feature.loadings']]) {
          check[['feature.loadings']] <- self[['reductions']][[reduc]]$exists(name = 'features')
        }
        if (check[['feature.loadings.projected']]) {
          check[['feature.loadings.projected']] <- self[['reductions']][[reduc]]$exists(name = 'projected.features')
        }
        index[[reduc.assay]][['reductions']][[reduc]] <- check
        if (IsGlobal(object = self[['reductions']][[reduc]])) {
          index$global$reductions <- c(index$global$reductions, reduc)
        }
      }
      # Get graph information
      for (graph in names(x = self[['graphs']])) {
        if (self[['graphs']][[graph]]$attr_exists(attr_name = 'assay.used')) {
          graph.assay <- h5attr(x = self[['graphs']][[graph]], which = 'assay.used')
          if (graph.assay %in% names(x = index)) {
            index[[graph.assay]]$graphs <- c(index[[graph.assay]]$graphs, graph)
          } else {
            warning(
              "Cannot find assay ",
              graph.assay,
              " in the h5Seurat file",
              call. = FALSE,
              immediate. = TRUE
            )
          }
        } else {
          index$no.assay$graphs <- c(index$no.assay$graphs, graph)
        }
      }
      # Get images
      if (version >= numeric_version(x = '3.1.4.9900')) {
        for (image in names(x = self[['images']])) {
          img.group <- self[['images']][[image]]
          if (!img.group$attr_exists(attr_name = 'assay') || !img.group$attr_exists(attr_name = 's4class')) {
            next
          }
          img.assay <- h5attr(x = img.group, which = 'assay')
          if (!img.assay %in% names(x = index)) {
            warning(
              "Cannot find assay ",
              img.assay,
              " in the H5Seurat file",
              call. = FALSE,
              immediate. = TRUE
            )
            index$no.assay$images <- c(index$no.assay$images, image)
          } else {
            index[[img.assay]]$images <- c(index[[img.assay]]$images, image)
          }
          if (IsGlobal(object = img.group)) {
            index$global$images <- c(index$global$images, image)
          }
        }
      }
      # Get commands
      for (cmd in names(x = self[['commands']])) {
        assay <- ifelse(
          test = self[['commands']][[cmd]]$attr_exists(attr_name = 'assay.used'),
          yes = h5attr(x = self[['commands']][[cmd]], which = 'assay.used'),
          no = NA_character_
        )
        if (assay %in% setdiff(x = names(x = index), y = c('global', 'no.assay'))) {
          index[[assay]]$commands <- c(index[[assay]]$commands, cmd)
        } else if (!is.na(x = assay)) {
          warning(
            "Cannot find assay",
            assay,
            " in the h5Seurat file",
            call. = FALSE,
            immediate. = TRUE
          )
        } else {
          index$no.assay$commands <- c(index$no.assay$commands, cmd)
        }
      }
      # TODO: Get metadata
      # TODO: Get miscellaneous data
      # TODO: Get tool-specific results
      # Finalize the index
      private$index.internal <- structure(
        .Data = index,
        class = c('h5SI', 'list'),
        # active.assay = DefaultAssay(object = self)
        active.assay = 'SCT'
      )
      return(invisible(x = NULL))
    },
    create = function(version, verbose = TRUE) {
      if (self$mode == 'r') {
        stop(private$errors(type = 'mode'), call. = FALSE)
      }
      version <- ClosestVersion(query = version, targets = private$versions)
      if (verbose) {
        message("Creating h5Seurat file for version ", version)
      }
      self$set.version(version = version)
      if (numeric_version(x = version) >= numeric_version(x = '3.1.2')) {
        for (group in c('assays', 'commands', 'graphs', 'misc', 'reductions', 'tools')) {
          if (!private$is.data(name = group, type = 'H5Group')) {
            self$create_group(name = group)
          }
        }
        attrs <- c(
          'active.assay' = '',
          'project' = 'SeuratDiskProject'
        )
        for (i in seq_along(along.with = attrs)) {
          if (!self$attr_exists(attr_name = names(x = attrs)[i])) {
            self$create_attr(
              attr_name = names(x = attrs)[i],
              robj = attrs[i],
              dtype = GuessDType(x = attrs[i])
            )
          }
        }
      }
      if (numeric_version(x = version) >= numeric_version(x = '3.1.4.9900')) {
        self$create_group(name = 'images')
      }
      return(invisible(x = self))
    },
    validate = function(verbose = TRUE, ...) {
      if (self$mode %in% modes$new) {
        private$create(
          version = packageVersion(pkg = 'Seurat'),
          verbose = verbose
        )
        return(invisible(x = NULL))
      }
      if (verbose) {
        message("Validating h5Seurat file")
      }
      if (!self$attr_exists(attr_name = 'version')) {
        stop(
          "Invalid h5Seurat file: cannot find attribute 'version'",
          call. = FALSE
        )
      }
      version <- h5attr(x = self, which = 'version')
      version <- ClosestVersion(query = version, targets = private$versions)
      version <- numeric_version(x = version)
      if (version >= numeric_version(x = '3.1.2')) {
        private$v3.1.2()
      }
      if (version >= numeric_version(x = '3.1.3.9900')) {
        private$v3.2.0()
      }
      private$build.index(version = as.character(x = version))
      return(invisible(x = NULL))
    },
    v3.1.2 = function() {
      # TODO: Check top-level attributes
      attrs <- c('project', 'active.assay', 'version')
      for (attr in attrs) {
        if (!self$attr_exists(attr_name = attr)) {
          stop("Missing attribute ", attr, call. = FALSE)
        }
      }
      # TODO: Check cell.names and meta.data
      if (!private$is.data(name = 'cell.names')) {
        stop("Cannot find dataset with cell names", call. = FALSE)
      }
      ncells <- self[['cell.names']]$dims
      if (length(x = ncells) != 1) {
        stop("Cell names must be one-dimensional", call. = FALSE)
      }
      if (private$is.data(name = 'meta.data')) {
        if (length(x = self[['meta.data']]$dims) != 1) {
          stop("Cell-level metadata must be one-dimensional")
        } else if (self[['meta.data']]$dims != ncells) {
          stop(
            "Cell number mismatch between cell names and cell-level metadata",
            call. = FALSE
          )
        }
        if (!inherits(x = self[['meta.data']]$get_type(), what = 'H5T_COMPOUND')) {
          stop("Cell-level metadata must be a data frame", call. = FALSE)
        }
      } else if (private$is.data(name = 'meta.data', type = 'H5Group')) {
        warning("Validation for group meta data not yet implemented", call. = FALSE, immediate. = TRUE)
      } else {
        stop("Cannot find cell-level metadata")
      }
      # TODO: Check Assays
      if (!private$is.data(name = 'assays', type = 'H5Group')) {
        stop("Cannot find assay expression data", call. = FALSE)
      }
      if (!DefaultAssay(object = self) %in% names(x = self[['assays']])) {
        stop("Default assay not present", call. = FALSE)
      }
      for (assay in names(x = self[['assays']])) {
        if (!private$is.data(name = file.path('assays', assay), type = 'H5Group')) {
          stop(
            "Assay representations must be HDF5 groups, offending entry: ",
            assay,
            call. = FALSE
          )
        }
      }
      # TODO: Check DimReducs
      # TODO: Check Graphs
      # TODO: Check SeuratCommands
      # TODO: Check miscellaneous data
      # TODO: Check tool-specific results
      return(invisible(x = NULL))
    },
    v3.2.0 = function() {
      return(invisible(x = NULL))
      .NotYetImplemented()
    }
  )
)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Tools for handling h5Seurat indexes (h5SI objects)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#' Tools for handling h5Seurat indexes
#'
#' @param x,object An h5Seurat index (\code{h5SI})
#'
#' @name h5SI
#' @rdname h5SI
#'
#' @seealso \code{\link[Seurat]{DefaultAssay}} \code{\link[base]{print}}
#'
#' @keywords internal
#'
NULL

#' @importFrom Seurat DefaultAssay
#'
#' @inheritParams Seurat::DefaultAssay
#'
#' @return \code{DefaultAssay}: the assay set as the default assay
#'
#' @rdname h5SI
#' @method DefaultAssay h5SI
#' @export
#'
DefaultAssay.h5SI <- function(object) {
  return(attr(x = object, which = 'active.assay'))
}

#' @inheritParams base::print
#'
#' @return \code{print}: Invisibly returns \code{x}
#'
#' @importFrom cli symbol
#' @importFrom tools toTitleCase
#' @importFrom crayon red green yellow col_align
#'
#' @rdname h5SI
#' @method print h5SI
#' @export
#'
print.h5SI <- function(x, ...) {
  # Some constants
  catn <- function(...) {
    cat(..., '\n', sep = '')
  }
  reduc.header <- c(
    'Embeddings' = 'cell.embeddings',
    'Loadings' = 'feature.loadings',
    'Projected' = 'feature.loadings.projected',
    'JackStraw' = 'jackstraw'
  )
  symbols <- c(red(symbol$cross), green(symbol$tick))
  # Get the assays
  assays <- setdiff(x = names(x = x), y = c('global', 'no.assay'))
  assays <- assays[order(assays == DefaultAssay(object = x), decreasing = TRUE)]
  for (assay in assays) {
    header <- paste("Data for assay", assay)
    if (assay == DefaultAssay(object = x)) {
      header <- paste0(header, yellow(symbol$star), ' (default assay)')
    }
    catn(header)
    # Show slot information
    catn(col_align(
      text = c('counts', 'data', 'scale.data'),
      width = nchar(x = 'scale.data') + 1,
      align = 'center'
    ))
    catn(col_align(
      text = symbols[x[[assay]]$slots + 1],
      width = nchar(x = 'scale.data') + 1,
      align = 'center'
    ))
    # Show dimensional reduction information
    if (!is.null(x = x[[assay]]$reductions)) {
      catn("Dimensional reductions:")
      reductions <- names(x = x[[assay]]$reductions)
      reductions <- paste0(' ', reductions, ': ')
      reductions <- col_align(text = reductions, width = max(nchar(x = reductions)))
      catn(
        MakeSpace(n = max(nchar(x = reductions))),
        col_align(
          text = names(x = reduc.header),
          width = max(nchar(x = names(x = reduc.header))) + 1,
          align = 'center'
        )
      )
      for (i in seq_along(along.with = reductions)) {
        reduc <- names(x = x[[assay]]$reductions)[i]
        catn(
          reductions[i],
          col_align(
            text = symbols[x[[assay]]$reductions[[reduc]] + 1],
            width = max(nchar(x = names(x = reduc.header))) + 1,
            align = 'center'
          )
        )
      }
    }
    # Show graph information
    if (!is.null(x = x[[assay]]$graphs)) {
      catn("Graphs:")
      catn(paste0(' ', symbol$line, ' ', x[[assay]]$graphs, collapse = '\n'))
    }
    # Show image information
    if (!is.null(x = x[[assay]]$images)) {
      catn("Images:")
      catn(paste0(' ', symbol$line, ' ', x[[assay]]$images, collapse = '\n'))
    }
    # TODO: Show command information
    # TODO: Show globals
    # if (!is.null(x = x$global)) {
    #   catn("Globally available information:")
    # }
    # Show no assay
  }
  return(invisible(x = x))
}
