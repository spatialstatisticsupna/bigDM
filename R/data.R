#' Spanish colorectal cancer mortality data
#'
#' @description \code{sf} object containing the polygons of the municipalities of continental Spain
#' and simulated colorectal cancer mortality data.
#'
#' @usage Carto_SpainMUN
#'
#' @format Formal class \code{sf}; the data contains a data.frame with 7907 rows and 11 variables.
#' \itemize{
#'   \item ID: character vector of geographic identifiers
#'   \item name: character vector of municipality names
#'   \item area: municipality polygon areas (in square meters)
#'   \item perimeter: municipality polygon perimeters (in meters)
#'   \item obs: observed number of cases
#'   \item exp: expected number of cases
#'   \item SMR: standardized mortality ratios
#'   \item region: character vector of autonomous regions
#'   \item geometry: sfc_MULTIPOLYGON
#' }
#' @name Carto_SpainMUN
#' @docType data
#' @keywords data
NULL

#' Spanish lung cancer mortality data
#'
#' @description \code{data.frame} object containing simulated lung cancer mortality data in the 7907 municipalities of continental Spain during the period 1991-2015.
#'
#' @usage Data_LungCancer
#'
#' @format Formal class \code{data.frame} with 197.675 rows and 5 colunmns.
#' \itemize{
#'   \item ID: character vector of geographic identifiers
#'   \item year: numeric vector of year's identifiers
#'   \item obs: observed number of cases
#'   \item exp: expected number of cases
#'   \item SMR: standardized mortality ratios
#'   \item pop: population at risk
#' }
#' @name Data_LungCancer
#' @docType data
#' @keywords data
NULL

#' Spanish cancer mortality data for the joint analysis of multiple diseases
#'
#' \code{data.frame} object containing simulated cancer mortality data for three diseases in the 7907 municipalities of continental Spain.
#'
#' @usage Data_MultiCancer
#'
#' @format Formal class \code{data.frame} with 237.271 rows and 5 colunmns.
#' \itemize{
#'   \item ID: character vector of geographic identifiers
#'   \item disease: numeric vector of disease identifiers
#'   \item obs: observed number of cases
#'   \item exp: expected number of cases
#'   \item SMR: standardized mortality ratios
#' }
#' @name Data_MultiCancer
#' @docType data
#' @keywords data
NULL

#' @docType package
#' @bibliography system.file("REFERENCES.bib", package = "bigDM")
