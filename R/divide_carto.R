#' Divide the spatial domain into subregions
#'
#' @description The function takes an object of class \code{SpatialPolygonsDataFrame} or \code{sf}
#' and divides it into subregions according to some grouping variable.
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
#' @param ID.group character vector of grouping identifiers.
#' @param k numeric value with the neighbourhood order to add polygons at the border of the spatial subdomains.
#' If k=0 (default) a disjoint partition is defined.
#' @param plot logical value (default \code{FALSE}), if \code{TRUE} then the spatial polygons within each subdomain are ploted.
#'
#' @return List of \code{sf} objects with the spatial polygons of each subdomain.
#'
#' @importFrom geos geos_prepared_intersects
#' @importFrom RColorBrewer brewer.pal
#' @importFrom sf st_as_sf st_set_geometry st_union
#' @importFrom stats aggregate
#'
#' @examples
#' \dontrun{
#' library(tmap)
#'
#' ## Load the Spain colorectal cancer mortality data ##
#' data(Carto_SpainMUN)
#'
#' ## Plot of the grouping variable 'region' ##
#' tmap4 <- packageVersion("tmap") >= "3.99"
#'
#' if(tmap4){
#'         tm_shape(Carto_SpainMUN) +
#'                 tm_polygons(fill="region",
#'                             fill.scale=tm_scale(values="brewer.set3"),
#'                             fill.legend=tm_legend(frame=FALSE))
#' }else{
#'         tm_shape(Carto_SpainMUN) +
#'                 tm_polygons(col="region") +
#'                 tm_layout(legend.outside=TRUE)
#' }
#'
#' ## Disjoint partition ##
#' carto.k0 <- divide_carto(carto=Carto_SpainMUN, ID.group="region", k=0)
#'
#' ## Partition + 1st order neighbours ##
#' carto.k1 <- divide_carto(carto=Carto_SpainMUN, ID.group="region", k=1)
#'
#' ## Partition + 2nd order neighbours ##
#' carto.k2 <- divide_carto(carto=Carto_SpainMUN, ID.group="region", k=2)
#'
#' ## Plot the spatial polygons for the autonomous region of Castilla y Leon ##
#' plot(carto.k2$`Castilla y Leon`$geometry, col="dodgerblue4", main="Castilla y Leon")
#' plot(carto.k1$`Castilla y Leon`$geometry, col="dodgerblue", add=TRUE)
#' plot(carto.k0$`Castilla y Leon`$geometry, col="lightgrey", add=TRUE)
#' }
#'
#' @export
divide_carto <- function(carto, ID.group=NULL, k=0, plot=FALSE){

        ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class
        carto <- sf::st_as_sf(carto)

        ## Construct the grouped 'sf' object by ID.group variable ##
        Data <- sf::st_set_geometry(carto, NULL)
        carto.group <- stats::aggregate(carto[,"geometry"], list(ID.group=Data[,ID.group]), utils::head)
        D <- nrow(carto.group)

        ########################
        ## Disjoint partition ##
        ########################
        group.names <- carto.group$ID.group
        names(group.names) <- group.names

        carto.k0 <- lapply(group.names, function(x) {
                aux <- carto[st_set_geometry(carto[,ID.group],NULL)==x,]
                rownames(aux) <- NULL
                aux
        })

        if(k==0){
                if(plot) lapply(carto.k0, function(x) plot(x$geometry, main=unique(st_set_geometry(x,NULL)[,ID.group])))
                return(carto.k0)
        }

        ############################################
        ## Partition including k-order neighbours ##
        ############################################
        if(k>0){

                carto.k <- vector("list",D)
                names(carto.k) <- names(carto.k0)

                if(plot) color <- RColorBrewer::brewer.pal(k+2,"Blues")

                for(i in 1:D){
                        aux.carto <- carto.group$geometry[i]

                        for(j in 1:k){
                                # loc <- sf::st_intersects(aux.carto, carto[,"geometry"])
                                # carto.k[[i]] <- merge(carto, data.frame(loc=unlist(loc)), by.x="temp", by.y="loc")
                                loc <- geos::geos_prepared_intersects(aux.carto, carto[,"geometry"])
                                carto.k[[i]] <- carto[loc,]
                                rownames(carto.k[[i]]) <- NULL

                                if(plot){
                                        plot(carto.k[[i]]$geometry, col=color[j+2], main=sort(unique(Data[,ID.group]))[i],
                                             xlim=sf::st_bbox(carto.group$geometry[i])[c(1,3)]*c(0.99,1.01),
                                             ylim=sf::st_bbox(carto.group$geometry[i])[c(2,4)]*c(0.99,1.01))

                                        plot(carto.k0[[i]]$geometry, col=color[1], add=TRUE)
                                }

                                aux.carto <- sf::st_union(carto.k[[i]])
                        }
                }

                return(carto.k)
        }
}
