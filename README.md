# survcomp #

Overview
--------

R package providing functions to assess and to compare the performance of risk prediction (survival) models.

Author: *Benjamin Haibe-Kains*, Markus Schroeder, Catharina Olsen, Christos Sotiriou, Gianluca Bontempi, John Quackenbush

*survcomp* package is also available on [Bioconductor](https://bioconductor.org/packages/release/bioc/html/survcomp.html#since).

Installation
------------

``` r
# Installing the development version from GitHub:
# install.packages("devtools")
devtools::install_github("bhklab/survcomp")

# Installing from Bioconductor
## try http:// if https:// URLs are not supported
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("survcomp")
```

Usage
-----

Extensive usage examples are provided in the survcomp vignette on Bioconductor: [survcomp.pdf](https://bioconductor.org/packages/release/bioc/vignettes/survcomp/inst/doc/survcomp.pdf)

Getting help
------------

To view documentation for the version of this package installed in your system, start R and enter:

``` r
browseVignettes("survcomp")
```

Contact us by filing an issue in the survcomp [issues](https://github.com/bhklab/survcomp/issues) page.
