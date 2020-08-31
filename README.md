
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SeuratDisk v0.0.0.9013

<!-- badges: start -->

[![CRAN/METACRAN](https://img.shields.io/cran/v/SeuratDisk)](https://cran.r-project.org/package=SeuratDisk)
[![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://github.com/mojaveazure/seurat-disk)
<!-- badges: end -->

<!-- Interfaces for HDF5-based Single Cell File Formats -->

The h5Seurat file format is specifically designed for the storage and
analysis of multi-modal single-cell and spatially-resolved expression
experiments, for example, from CITE-seq or 10X Visium technologies. It
holds all molecular information and associated metadata, including (for
example) nearest-neighbor graphs, dimensional reduction information,
spatial coordinates and image data, and cluster labels. We also support
rapid and on-disk conversion between h5Seurat and AnnData objects, with
the goal of enhancing interoperability between Seurat and Scanpy.

## Installation

SeuratDisk is not currently available on CRAN. You can install it from
[GitHub](https://github.com/mojaveazure/seurat-disk) with:

``` r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
```

## Dependencies

SeuratDisk depends on the following non-standard packages:

| Package |                   CRAN Webpage                    |                      Source                      |                     Website                     |
| :-----: | :-----------------------------------------------: | :----------------------------------------------: | :---------------------------------------------: |
|   cli   |  [CRAN](https://cran.r-project.org/package=cli)   |  [GitHub](https://github.com/r-lib/cli#readme)   |                        –                        |
| crayon  | [CRAN](https://cran.r-project.org/package=crayon) | [GitHub](https://github.com/r-lib/crayon#readme) |                        –                        |
|  hdf5r  | [CRAN](https://cran.r-project.org/package=hdf5r)  |   [GitHub](https://github.com/hhoeflin/hdf5r)    |  [Website](https://hhoeflin.github.io/hdf5r/)   |
| Matrix  | [CRAN](https://cran.r-project.org/package=Matrix) |                        –                         | [Website](http://Matrix.R-forge.R-project.org/) |
|   R6    |   [CRAN](https://cran.r-project.org/package=R6)   |      [GitHub](https://github.com/r-lib/R6/)      |         [Website](https://r6.r-lib.org)         |
|  rlang  | [CRAN](https://cran.r-project.org/package=rlang)  |     [GitHub](https://github.com/r-lib/rlang)     |        [Website](http://rlang.r-lib.org)        |
| Seurat  | [CRAN](https://cran.r-project.org/package=Seurat) |  [GitHub](https://github.com/satijalab/seurat)   |   [Website](http://www.satijalab.org/seurat)    |
|  withr  | [CRAN](https://cran.r-project.org/package=withr)  |                        –                         |        [Website](http://withr.r-lib.org)        |
