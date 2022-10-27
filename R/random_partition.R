#' Define a random partition of the spatial domain based on a regular grid
#'
#' @description The function takes an object of class \code{SpatialPolygonsDataFrame} or \code{sf} and
#' defines a random partition of the spatial polygons based on a regular grid over the whole domain
#' using the \code{st_make_grid} function of the \code{sf} package.
#'
#' @details After defining a random partition of the spatial polygons based on a regular grid, the subregions with number of areas smaller than the value given by the \code{min.size} are merged to its nearest neighbour.
#' Then, the subregions with number of areas greater than the value given by the \code{max.size} argument are divided.
#' Finally, if \code{prop.zero} argument is set, the subregions with proportion of areas with zero cases below that threshold are merged to its smallest neighbour.
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
#' @param rows integer; number of rows to define the regular grid. Default to 3.
#' @param columns integer; number of columns to define the regular grid. Default to 3.
#' @param min.size numeric; value to fix the minimum number of areas in each spatial partition (if \code{NULL}, this step is skipped). Default to 50.
#' @param max.size numeric; value to fix the maximum number of areas in each spatial partition (if \code{NULL}, this step is skipped). Default to 600.
#' @param prop.zero numeric; value between 0 and 1 that indicates the maximum proportion of areas with no cases for each spatial partition.
#' @param O character; name of the variable that contains the observed number of disease cases for each areal units. Only required if \code{prop.zero} argument is set.
#'
#' @return \code{sf} object with the original data and a grouping variable named 'ID.group'
#'
#' @importFrom sf st_as_sf st_bbox st_centroid st_geometry st_intersects st_make_grid
#' @importFrom spdep knearneigh knn2nb
#' @importFrom stats aggregate
#'
#' @examples
#' \dontrun{
#' library(tmap)
#'
#' ## Load the Spain colorectal cancer mortality data ##
#' data(Carto_SpainMUN)
#'
#' ## Random partition based on a 3x3 regular grid (with no size restrictions) ##
#' carto.r1 <- random_partition(carto=Carto_SpainMUN, rows=3, columns=3,
#'                              min.size=NULL, max.size=NULL)
#' table(carto.r1$ID.group)
#'
#' part1 <- aggregate(carto.r1[,"geometry"], by=list(ID.group=carto.r1$ID.group), head)
#'
#' tm_shape(carto.r1) +
#'   tm_polygons(col="ID.group") +
#'   tm_shape(part1) + tm_borders(col="black", lwd=2) +
#'   tm_layout(main.title="3x3 regular grid (with no size restrictions)",
#'             main.title.position="center", main.title.size=1,
#'             legend.outside=TRUE)
#'
#'
#' ## Random partition based on a 6x4 regular grid (with size restrictions) ##
#' carto.r2 <- random_partition(carto=Carto_SpainMUN, rows=6, columns=4,
#'                              min.size=50, max.size=600)
#' table(carto.r2$ID.group)
#'
#' part2 <- aggregate(carto.r2[,"geometry"], by=list(ID.group=carto.r2$ID.group), head)
#'
#' tm_shape(carto.r2) +
#'   tm_polygons(col="ID.group") +
#'   tm_shape(part2) + tm_borders(col="black", lwd=2) +
#'   tm_layout(main.title="6x4 regular grid (min.size=50, max.size=600)",
#'             main.title.position="center", main.title.size=1,
#'             legend.outside=TRUE)
#'
#'
#' ## Random partition based on a 6x4 regular grid (with size and proportion of zero restrictions) ##
#' carto.r3 <- random_partition(carto=Carto_SpainMUN, rows=6, columns=4,
#'                              min.size=50, max.size=600, prop.zero=0.5, O="obs")
#' table(carto.r3$ID.group)
#'
#' part3 <- aggregate(carto.r3[,"geometry"], by=list(ID.group=carto.r3$ID.group), head)
#'
#' tm_shape(carto.r3) +
#'   tm_polygons(col="ID.group") +
#'   tm_shape(part3) + tm_borders(col="black", lwd=2) +
#'   tm_layout(main.title="6x4 regular grid (min.size=50, max.size=600, prop.zero=0.5)",
#'             main.title.position="center", main.title.size=1,
#'             legend.outside=TRUE)
#' }
#'
#' @export
random_partition <- function(carto, rows=3, columns=3, min.size=50, max.size=1000, prop.zero=NULL, O=NULL){

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
        carto$ID.group <- factor(as.numeric(factor(ID.group)))

        ## Merge the subregions with lower number of areas than min.size ##
        if(!is.null(min.size)){
                partition.size <- table(carto$ID.group)

                while(any(partition.size<min.size)){
                        cat(sprintf("+ Merging small subregions (min.size=%d)\n",min.size))

                        data <- st_set_geometry(carto,NULL)
                        partition <- stats::aggregate(carto[,"geometry"], list(ID.group=data$ID.group), head)

                        pos <- which(partition.size<min.size)
                        knn.nb <- spdep::knn2nb(spdep::knearneigh(sf::st_centroid(sf::st_geometry(partition), of_largest_polygon=TRUE), k=1))

                        for(i in pos){
                                carto$ID.group[carto$ID.group==i] <- knn.nb[[i]]
                        }

                        carto$ID.group <- factor(as.numeric(factor(carto$ID.group)))
                        partition.size <- table(carto$ID.group)
                }
        }

        ## Divide the subregions with greater number of areas than max.size ##
        if(!is.null(max.size)){
                while(any(partition.size>max.size)){
                        cat(sprintf("+ Dividing big subregions (max.size=%d)\n",max.size))

                        carto$ID.group <- as.numeric(carto$ID.group)
                        pos <- which(partition.size>max.size)

                        for(i in pos){
                                carto.aux <- carto[carto$ID.group==i,]
                                bbox <- sf::st_bbox(carto.aux)
                                largest.dim <- as.numeric(which.max(c(bbox["xmax"]-bbox["xmin"],bbox["ymax"]-bbox["ymin"])))

                                if(largest.dim==1) carto.grid <- sf::st_make_grid(carto.aux, n=c(2,1))
                                if(largest.dim==2) carto.grid <- sf::st_make_grid(carto.aux, n=c(1,2))

                                aux <- sf::st_centroid(sf::st_geometry(carto.aux), of_largest_polygon=TRUE)
                                aux <- sf::st_intersects(carto.grid,carto.aux)

                                ID.aux <- numeric()
                                for(j in 1:length(aux)){
                                        ID.aux[aux[[j]]] <- j
                                }
                                ID.aux <- ID.aux+max(as.numeric(carto$ID.group))

                                carto$ID.group[carto$ID.group==i] <- ID.aux
                                summary(as.factor(carto$ID.group))
                        }

                        carto$ID.group <- factor(as.numeric(factor(carto$ID.group)))
                        partition.size <- table(carto$ID.group)
                }
        }

        ## Check if proportion of areas with no cases for each subregion is below the maximum value given by prop.zero ##
        if(!is.null(prop.zero)){
                cat(sprintf("+ Checking if the proportion of areas with no cases for each subregion is below prop.zero=%g\n",prop.zero))
                if(is.null(O)) stop("WARNING: the 'O' argument is missing")

                partition <- aggregate(carto[,O], by=list(ID.group=carto$ID.group), function(x) mean(x==0))
                pos <- which(sf::st_set_geometry(partition, NULL)[,O]>prop.zero)

                it <- 1
                while(length(pos)>0){
                        cat(sprintf("  -> Iteration %d: %d subregion(s) are been merged\n",it,length(pos)))
                        knn.nb <- spdep::knn2nb(spdep::knearneigh(sf::st_centroid(sf::st_geometry(partition), of_largest_polygon=TRUE), k=4))

                        for(i in pos){
                                carto$ID.group[carto$ID.group==i] <- as.numeric(names(which.min(partition.size[knn.nb[[i]]])))
                        }

                        carto$ID.group <- factor(as.numeric(factor(carto$ID.group)))

                        partition <- aggregate(carto[,O], by=list(ID.group=carto$ID.group), function(x) mean(x==0))
                        pos <- which(sf::st_set_geometry(partition, NULL)[,O]>prop.zero)

                        it <- it+1
                }
        }

        if(any(table(carto$ID.group)>max.size)) warning(sprintf("%d subregion(s) have more than %d areas",sum(table(carto$ID.group)>max.size),max.size), call.=FALSE)

        return(carto)
}
