#' Add isolated areas (polygons) to its nearest neighbour
#'
#' @description The function returns a neighbour list of class \code{nb} and its associated spatial adjacency matrix
#' computed by adding isolated areas to its nearest neighbour (in terms of Euclidean distance between centroids) using the \code{knearneigh} function of 'spdep' package.
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
#' @param nb optional argument with the neighbour list of class \code{nb}.
#' If \code{NULL} (default), this object is computed from the \code{carto} argument.
#' @param plot logical value (default \code{FALSE}), if \code{TRUE} then the computed neighbourhood graph is plotted.
#'
#' @return This function returns a list with the following two elements:
#' \itemize{
#'   \item \code{nb}: the modified neighbours's list
#'   \item \code{W}: associated spatial adjacency matrix of class \code{dgCMatrix}
#' }
#'
#' @importFrom methods as
#' @importFrom sf st_as_sf st_centroid st_geometry
#' @importFrom spdep card knearneigh knn2nb nb2listw
#'
#' @examples
#' library(spdep)
#'
#' ## Load the Spanish colorectal cancer mortality data ##
#' data(Carto_SpainMUN)
#'
#' ## Compute the neighbour list from spatial polygons ##
#' nb_SpainMUN <- poly2nb(Carto_SpainMUN)
#' summary(nb_SpainMUN) # 1 region with no links
#'
#' ## Add isolated area to its nearest neighbour ####
#' carto.mod <- add_neighbour(carto=Carto_SpainMUN, nb=nb_SpainMUN)
#' summary(carto.mod$nb) # 0 region with no links
#'
#' @export
add_neighbour <- function(carto, nb=NULL, plot=FALSE){

  ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class
  carto <- sf::st_as_sf(carto)

  ## Compute the neighbours list of class 'nb'
  if(is.null(nb)) nb <- spdep::poly2nb(carto)

  ## Search for regions with no links
  card.nb <- spdep::card(nb)
  cat(" ",sum(card.nb==0),"region(s) with no links\n")

  if(any(card.nb==0)){

    ## Add to these regions its nearest neighbour
    knn.nb <- spdep::knn2nb(spdep::knearneigh(sf::st_centroid(sf::st_geometry(carto), of_largest_polygon=TRUE), k=1))

    for(i in which(card(nb)==0)){
      j <- knn.nb[[i]]
      nb[[i]] <- as.integer(j)
      nb[[j]] <- as.integer(unique(sort(c(i,nb[[j]]))))
    }
  }

  ## Plot the spatial polygons and the computed neighbourhood graph
  if(plot){
    plot(carto$geometry)
    plot(nb, sf::st_centroid(sf::st_geometry(carto), of_largest_polygon=TRUE), pch=19, cex=0.5, col="red", add=TRUE)
  }

  ## Compute the spatial adjacency matrix
  W <- spdep::nb2listw(nb, style="B")
  W <- methods::as(W,"CsparseMatrix")

  return(list(nb=nb,W=W))
}
