#' Merge disjoint connected subgraphs
#'
#' @description The function returns a neighbour list of class \code{nb} and its associated spatial adjacency matrix
#' computed by merging disjoint connected subgraphs through its nearest polygon centroids.
#'
#' @details This function first calls the \code{\link{add_neighbour}} function to search for isolated areas.
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
#' @param ID.area character vector of geographic identifiers.
#' @param nb optional argument with the neighbours list of class \code{nb}.
#' If \code{NULL} (default), this object is computed from the \code{carto} argument.
#' @param plot logical value (default \code{FALSE}), if \code{TRUE} then the computed neighbourhood graph is ploted.
#'
#' @return This function returns a list with the following two elements:
#' \itemize{
#'   \item \code{nb}: the modified neighbours list
#'   \item \code{W}: associated spatial adjacency matrix of class \code{CsparseMatrix}
#' }
#'
#' @import spatialreg
#' @importFrom methods as
#' @importFrom sf st_as_sf st_centroid st_distance st_drop_geometry st_geometry
#' @importFrom spdep card n.comp.nb poly2nb
#'
#' @examples
#' library(spdep)
#'
#' ## Load the Spain colorectal cancer mortality data ##
#' data(Carto_SpainMUN)
#'
#' ## Select the polygons (municipalities) of the 'Comunidad Valenciana' region ##
#' carto <- Carto_SpainMUN[Carto_SpainMUN$region=="Comunidad Valenciana",]
#'
#' carto.nb <- poly2nb(carto)
#' n.comp.nb(carto.nb)$nc # 2 disjoint connected subgraphs
#'
#' ## Plot the spatial polygons and its neighbourhood graph
#' op <- par(mfrow=c(1,2), pty="s")
#'
#' plot(carto$geometry, main="Original neighbourhood graph")
#' plot(carto.nb, st_centroid(st_geometry(carto), of_largest_polygon=TRUE),
#'      pch=19, cex=0.5, col="red", add=TRUE)
#'
#' ## Use the 'connect_subgraphs' function ##
#' carto.mod <- connect_subgraphs(carto=carto, ID.area="ID", nb=carto.nb, plot=TRUE)
#' title(main="Modified neighbourhood graph")
#'
#' n.comp.nb(carto.mod$nb)$nc==1
#'
#' par(op)
#'
#' @export
connect_subgraphs <- function(carto, ID.area=NULL, nb=NULL, plot=FALSE){

  ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class
  carto <- sf::st_as_sf(carto)

  ## Compute the neighbours list of class 'nb'
  if(is.null(nb)) nb <- spdep::poly2nb(carto)

  ## Search for regions with no links ##
  cat("Searching for isolated areas:\n")
  nb <- add_neighbour(carto=carto, nb=nb, plot=FALSE)$nb

  ## Search for disjoint connected subgraphs ##
  nc <- spdep::n.comp.nb(nb)$nc

  cat("\nSearching for disjoint connected subgraphs:\n")
  if(nc==1){
    cat(" No disjoint connected subgraphs\n")
  }else{
    cat(" ",nc,"disjoint connected subgraphs\n")
  }

  while(nc>1){

    ## Compute distance matrix between centroids ##
    dist.matrix <- sf::st_distance(sf::st_centroid(sf::st_geometry(carto), of_largest_polygon=TRUE))

    rownames(dist.matrix) <- st_drop_geometry(carto)[,ID.area]
    colnames(dist.matrix) <- st_drop_geometry(carto)[,ID.area]

    smallest.subgraph <- which.min(table(spdep::n.comp.nb(nb)$comp.id))
    smallest.loc <- which(spdep::n.comp.nb(nb)$comp.id==smallest.subgraph)

    ## Distance matrix between centroids of smallest subgraph and the rest of polygons ##
    smallest.dist <- dist.matrix[smallest.loc,-smallest.loc]
    min.pos <- which.min(smallest.dist)

    pos <- c((min.pos-1) %% nrow(smallest.dist), (min.pos-1) %/% nrow(smallest.dist))+1

    i.ID <- rownames(smallest.dist)[pos[1]]
    j.ID <- colnames(smallest.dist)[pos[2]]

    i <- which(rownames(dist.matrix)==i.ID)
    j <- which(colnames(dist.matrix)==j.ID)

    nb[[i]] <- as.integer(sort(c(j,nb[[i]])))
    nb[[j]] <- as.integer(sort(c(nb[[j]],i)))

    nc <- spdep::n.comp.nb(nb)$nc
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
