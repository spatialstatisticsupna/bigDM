#' Define a random partition of the spatial domain based on a regular grid
#'
#' @description The function takes an object of class \code{SpatialPolygonsDataFrame} or \code{sf} and
#' defines a random partition of the spatial polygons based on a regular grid over the whole domain
#' using the \code{st_make_grid} function of 'sf' package.
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
#' @param rows numeric; number of rows to define the regular grid. Defaults to 3.
#' @param columns numeric; number of columns to define the regular grid. Defaults to 3.
#'
#' @return \code{sf} object with the original data and a grouping variable named 'ID.group'
#'
#' @importFrom sf st_as_sf st_centroid st_geometry st_intersects st_make_grid
#'
#' @examples
#' library(tmap)
#'
#' ## load the Spain colorectal cancer mortality data ##
#' data(Carto_SpainMUN)
#'
#' ## define a random partition based on a 3x3 regular grid ##
#' carto.new <- random_partition(carto=Carto_SpainMUN, rows=3, columns=3)
#'
#' ## plot of the grouping variable 'ID.group' ##
#' tm_shape(carto.new) +
#'   tm_polygons(col="ID.group") +
#'   tm_layout(legend.outside=TRUE)
#'
#' @export
random_partition <- function(carto, rows=3, columns=3){

        ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class
        carto <- sf::st_as_sf(carto)

        ## Create the regular grid over the whole spatial domain ##
        carto.grid <- sf::st_make_grid(carto, n=c(columns,rows))

        ## Compute the intersection between polygon centroids and the regular grid ##
        aux <- sf::st_centroid(sf::st_geometry(carto), of_largest_polygon=TRUE)
        aux <- sf::st_intersects(carto.grid,aux)

        ## Check for errors ##
        if(nrow(carto)!=sum(unlist(lapply(aux, function(x) length(x))))){
                stop("WARNING: The spatial partition is not well defined.")
        }

        ## Add the grouping variable to the input data ##
        ID.group <- numeric()
        for(i in 1:length(aux)){
                ID.group[aux[[i]]] <- i
        }
        carto$ID.group <- factor(ID.group)

        return(carto)
}
