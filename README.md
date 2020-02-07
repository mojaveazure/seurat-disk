
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SeuratDisk v0.0.0.9000

<!-- badges: start -->

![CRAN/METACRAN](https://img.shields.io/cran/v/SeuratDisk)
<!-- badges: end -->

Interfaces for HDF5-based Single Cell File Formats

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

SeuratDisk depends on the following packages:

  - hdf5r [:link:](https://cran.r-project.org/package=hdf5r)
    [:octocat:](https://github.com/hhoeflin/hdf5r)
    [:page\_facing\_up:](https://hhoeflin.github.io/hdf5r)
  - methods
  - R6 [:link:](https://cran.r-project.org/package=R6)
    [:octocat:](https://github.com/r-lib/R6/)
    [:page\_facing\_up:](https://r6.r-lib.org)
  - rlang [:link:](https://cran.r-project.org/package=rlang)
    [:octocat:](https://github.com/r-lib/rlang)
    [:page\_facing\_up:](http://rlang.r-lib.org)
  - Seurat [:link:](https://cran.r-project.org/package=Seurat)
    [:octocat:](https://github.com/satijalab/seurat)
    [:page\_facing\_up:](http://www.satijalab.org/seurat)
  - tools
