#' Scalable Bayesian Disease Mapping Models for High-Dimensional Data
#'
#' @description This package implements several (scalable) spatial Poisson mixed models for high-dimensional areal count data
#' in a fully Bayesian setting using the integrated nested Laplace approximation (INLA) technique.
#'
#' @details Below, there is a list with a brief overview of all package functions:
#' \tabular{ll}{
#'   \code{ \link{add_neighbour}}\tab Adds isolated areas (polygons) to its nearest neighbour \cr
#'   \code{\link{CAR_INLA}} \tab  Fits several spatial CAR models for high-dimensional count data \cr
#'   \code{\link{connect_subgraphs}} \tab Merges disjoint connected subgraphs \cr
#'   \code{\link{divide_carto}} \tab Divides the spatial domain into subregions \cr
#'   \code{\link{mergeINLA}} \tab Merges \code{inla} objects for partition models \cr
#'   \code{\link{random_partition}} \tab Defines a random partition of the spatial domain based on a regular grid \cr
#'   ----------------------\tab ---------------------------------------------------------------------------------- \cr
#'   }
#'
#' @author
#' Maintainer: Aritz Adin <aritz.adin@unavarra.es>
#' \cr \cr
#' This work has been supported by the Spanish Ministry of Science and Innovation (project MTM 2017-82553-R (AEI/FEDER, UE)).
#' It has also been partially funded by la Caixa Foundation (ID 1000010434), Caja Navarra Foundation, and UNED Pamplona, under agreement LCF/PR/PR15/51100007 (project REF P/13/20).
#'
#' @references \insertRef{orozco2020}{bigDM}
#'
#' @seealso See the following vignettes for further details and examples using this package:
#' \enumerate{
#'    \item{\href{https://emi-sstcdapp.unavarra.es/bigDM/bigDM-1-fitting-spatial-models.html}{\code{bigDM: fitting spatial models}}}
#'    \item{\href{https://emi-sstcdapp.unavarra.es/bigDM/bigDM-2-parallel-and-distributed-modelling.html}{\code{bigDM: parallel and distributed modelling}}}
#' }
#'
#' @examples
#' ## See the examples for CAR_INLA
#'
#' @docType package
#' @name bigDM-package
#' @aliases bigDM
NULL
