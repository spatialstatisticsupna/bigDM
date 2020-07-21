#' Spanish colorectal cancer mortality data
#'
#' @description \code{sf} object containing the polygons of the municipalities of continental Spain
#' and simulated colorectal cancer mortality data.
#'
#' @usage Carto_SpainMUN
#'
#' @format Formal class \code{sf}; the data contains a data.frame with 7907 rows and 10 variables.
#' \itemize{
#'   \item ID: character vector of geographic identifiers
#'   \item name: character vector of municipality names
#'   \item lat: numeric vector of longitude values
#'   \item long: numeric vector of latitude values
#'   \item area: municipality polygon areas in square meters
#'   \item perimeter: municipality polygon perimeters in degree units
#'   \item obs: observed number of cases
#'   \item exp: expected number of cases
#'   \item region: character vector of autonomous region
#'   \item geometry: sfc_MULTIPOLYGON
#' }
#' @name Carto_SpainMUN
#' @docType data
#' @keywords data
NULL

#' @docType package
#' @bibliography system.file("REFERENCES.bib", package = "bigDM")
