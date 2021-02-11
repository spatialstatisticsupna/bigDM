# bigDM
Scalable Bayesian disease mapping models for high-dimensional data.

## Table of contents

- [The package](#the-package)
- [Installation](#installation)
- [Basic Use](#basic-use)
- [Copyright and license](#copyright-and-license)


# The package
This package implements several (scalable) spatial generalised linear mixed models for high-dimensional areal count
data, with inference in a fully Bayesian setting using the integrated nested Laplace approximation
(INLA) technique.

Below, there is a list with a brief overview of all package functions:

* ```add_neighbour``` Adds isolated areas (polygons) to its nearest neighbour.
* ```CAR_INLA``` Fits several spatial CAR models for high-dimensional count data.
* ```connect_subgraphs``` Merges disjoint connected subgraphs.
* ```divide_carto``` Divides the spatial domain into subregions.
* ```mergeINLA``` Merges inla objects for partition models.
* ```random_partition``` Defines a random partition of the spatial domain based on a regular grid.


# Installation

[Installing Rtools40 for Windows](https://cran.r-project.org/bin/windows/Rtools/)

R version 4.0.0 and newer for Windows requires the new Rtools40 to build R packages with C/C++/Fortran code from source.


## Install from CRAN
Not in CRAN

## Install from GitHub
```
# Install devtools package from CRAN repository
install.packages("devtools")

# Load devtools library
library(devtools)

# Install the R-INLA package
install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# Install bigDM from GitHub repositoy
install_github("spatialstatisticsupna/bigDM")
```

_*Important note:_ The __bigDM__ package uses the `parallel:makeCluster()` function to speed up computations when fitting scalable models. It seems that some unresolved issues have been detected when running this function in RStudio on macOS with R 4.0.x (https://github.com/rstudio/rstudio/issues/6692)

## Basic Use
See the [vignette](https://emi-sstcdapp.unavarra.es/bigDM/bigDM.html) for further details and examples using this package.

When using this package, please cite the following paper:

[Orozco-Acosta, E., Adin, A., and Ugarte, M.D. (2021). Scalable Bayesian modeling for smoothing disease risks in large spatial data sets using INLA. _Spatial Statistics_, __41__, 100496.] (https://doi.org/10.1016/j.spasta.2021.100496)


## Copyright and license
Licensed under the GPL-3 License. [Full license here](/LICENSE.md).
