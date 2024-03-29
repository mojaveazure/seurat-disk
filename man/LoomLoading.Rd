% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoadLoom.R
\name{LoomLoading}
\alias{LoomLoading}
\alias{LoadLoom0.1}
\alias{LoadLoom3.0}
\title{Loom-file Loading}
\usage{
LoadLoom0.1(
  file,
  assay = NULL,
  cells = "col_atts/CellID",
  features = "row_attrs/Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE
)

LoadLoom3.0(
  file,
  assay = NULL,
  cells = "col_attrs/CellID",
  features = "row_attrs/Gene",
  normalized = NULL,
  scaled = NULL,
  filter = c("cells", "features", "all", "none"),
  verbose = TRUE
)
}
\arguments{
\item{assay}{Name of assay to store expression data as; if \code{NULL}, will
search for an HDF5 attribute named \code{SEURAT_ASSAY} or an attribute
dataset named \code{/attrs/SEURAT_ASSAY} for assay name. If not found,
defaults to \dQuote{RNA}}

\item{cells}{Name of dataset in \code{/col_attrs} with cell names}

\item{features}{Name of dataset in \code{/row_attrs} with feature names}

\item{normalized}{Name of matrix in \code{/layers} to store normalized data
as; pass \dQuote{/matrix} to store \code{/matrix} as normalized data instead
of raw counts}

\item{scaled}{Name of dataset in \code{/layers} to store scaled data as}

\item{filter}{Keep only selected cells and/or features as specified by
\code{/col_attrs/Valid} and \code{/row_attrs/Valid}, respectively}

\item{verbose}{Show progress updates}
}
\value{
A \code{\link[Seurat]{Seurat}} object
}
\description{
Version-specific loom-file loading functions
}
\details{
\code{LoadLoom} will try to automatically fill slots of a \code{Seurat}
object based on data presence or absence in a given loom file. This method
varies by loom specification version. For version-specific details, see
sections below
}
\section{Loom 0.1 Loading}{

Loading data from loom files less than version 3.0.0 is not
currently supported
}

\section{Loom 3.0.0 Loading}{

blah
}

\keyword{internal}
