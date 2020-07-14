# bigDM
Scalable Bayesian disease mapping models for high-dimensional data.

## Table of contents

- [The package](#the-package)
- [Installation](#installation)
- [Basic Use](#basic-use)
- [Copyright and license](#copyright-and-license)


# The package
This package implements several (scalable) spatial generalised linear mixed models to areal count
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
## Install from CRAN
Not in CRAN

## Install from GitHub
```
# Install devtools package from cran repository
install.packages("devtools")

# load devtools library
library(devtools)

# Install rsat from GitHub repositoy
install_github("spatialstatisticsupna/bigDM")
```

## Basic Use
See the [vignette](/) for further details and examples using this package.


## Copyright and license
Licensed under the GPL-3 License. [Full license here](/LICENSE.md).
