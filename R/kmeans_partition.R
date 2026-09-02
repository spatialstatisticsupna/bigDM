#' Define a random spatial partition of the domain based on k-means clustering of polygon centroids
#'
#' @description The function takes an object of class \code{SpatialPolygonsDataFrame} or \code{sf} and
#' defines a random spatial partition of the polygons using the k-means algorithm applied to
#' the longitude and latitude coordinates of their centroids.
#'
#' @details The k-means algorithm is applied to the longitude and latitude coordinates of the
#' polygon centroids to obtain spatially compact clusters. The resulting partition may
#' depend on the CRS of \code{carto}, particularly for large geographic domains.
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
#' @param centers integer; number of clusters into which the spatial domain is partitioned. Default to 10.
#' @param min.size numeric; value to fix the minimum number of areas in each spatial partition (if \code{NULL}, this step is skipped). Default to 150.
#' @param max.size numeric; value to fix the maximum number of areas in each spatial partition (if \code{NULL}, this step is skipped). Default to 600.
#' @param prop.zero numeric; value between 0 and 1 that indicates the maximum proportion of areas with no cases for each spatial partition.
#' @param O character; name of the variable that contains the observed number of disease cases for each areal units. Only required if \code{prop.zero} argument is set.
#' @param ... additional arguments passed to \code{\link[stats]{kmeans}}.
#'
#' @return \code{sf} object with the original data and a grouping variable named 'ID.group'
#'
#' @importFrom sf st_as_sf st_centroid st_coordinates st_drop_geometry
#' @importFrom spdep knearneigh knn2nb
#' @importFrom stats kmeans
#'
#' @seealso
#' \code{\link{grid_partition}} for a spatial partition based on a regular grid.
#'
#' @examples
#' \dontrun{
#' library(tmap)
#'
#' ## Load the Spain colorectal cancer mortality data ##
#' data(Carto_SpainMUN)
#'
#  ##' Random partition based on 10 clusters (with no size restrictions) ##
#' set.seed(1234)
#' carto.r1 <- kmeans_partition(carto=Carto_SpainMUN, centers=10,
#'                              min.size=NULL, max.size=NULL)
#' table(carto.r1$ID.group)
#'
#' part1 <- aggregate(carto.r1[,"geometry"], by=list(ID.group=carto.r1$ID.group), head)
#'
#' tm_shape(carto.r1) +
#'         tm_polygons(fill="ID.group",
#'                     fill.scale=tm_scale(values="brewer.set3"),
#'                     fill.legend=tm_legend(frame=FALSE)) +
#'         tm_shape(part1) + tm_borders(col="black", lwd=2) +
#'         tm_title(text="Random partition with 10 clusters (with no size restrictions)")
#'
#'
#' ## Random partition based on 20 clusters (with size restrictions) ##
#' set.seed(1234)
#' carto.r2 <- kmeans_partition(carto=Carto_SpainMUN, centers=20,
#'                              min.size=300, max.size=600)
#'
#' table(carto.r2$ID.group)
#'
#' part2 <- aggregate(carto.r2[,"geometry"], by=list(ID.group=carto.r2$ID.group), head)
#'
#' tm_shape(carto.r2) +
#'         tm_polygons(fill="ID.group",
#'                     fill.scale=tm_scale(values="brewer.set3"),
#'                     fill.legend=tm_legend(frame=FALSE)) +
#'         tm_shape(part2) + tm_borders(col="black", lwd=2) +
#'         tm_title(text="Random partition with 20 clusters (min.size=300, max.size=600)")
#'
#'
#' ## Random partition based on 30 clusters (with size and proportion of zero restrictions) ##
#' carto.r3 <- kmeans_partition(carto=Carto_SpainMUN, centers=30,
#'                              min.size=150, max.size=600, prop.zero=0.5, O="obs")
#'
#' table(carto.r3$ID.group)
#'
#' part3 <- aggregate(carto.r3[,"geometry"], by=list(ID.group=carto.r3$ID.group), head)
#'
#' tm_shape(carto.r3) +
#'         tm_polygons(fill="ID.group",
#'                     fill.scale=tm_scale(values="brewer.set3"),
#'                     fill.legend=tm_legend(frame=FALSE)) +
#'         tm_shape(part3) + tm_borders(col="black", lwd=2) +
#'         tm_title(text="Random partition with 30 clusters
#'                        (min.size=150, max.size=600, prop.zero=0.5)")
#' }
#'
#' @export
kmeans_partition <- function(carto, centers=10, min.size=150, max.size=600, prop.zero=NULL, O=NULL, ...){

        ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class
        carto <- sf::st_as_sf(carto)

        ## Extract the coordinates of the polygon centroids
        coords <- suppressWarnings(sf::st_centroid(carto, of_largest_polygon=TRUE)) |>
                sf::st_coordinates()

        ## Apply the k-means algorithm
        cl <- stats::kmeans(coords, centers=centers, ...)$cluster
        partition.size <- table(cl)

        carto$ID.group <- factor(cl)

        ## Merge the subregions with lower number of areas than min.size ##
        if(!is.null(min.size)){
                while(any(partition.size<min.size)){
                        cat(sprintf("+ Merging small subregions (min.size=%d)\n",min.size))

                        data <- sf::st_drop_geometry(carto, NULL)
                        partition <- stats::aggregate(carto[,"geometry"], list(ID.group=data$ID.group), head)

                        pos <- which(partition.size<min.size)
                        knn.nb <- suppressWarnings(
                                spdep::knn2nb(spdep::knearneigh(sf::st_centroid(sf::st_geometry(partition), of_largest_polygon=TRUE), k=1))
                        )

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
                pos <- which(sf::st_drop_geometry(partition)[,O]>prop.zero)

                it <- 1
                while(length(pos)>0){
                        cat(sprintf("  -> Iteration %d: %d subregion(s) are been merged\n",it,length(pos)))
                        knn.nb <- suppressWarnings(
                                spdep::knn2nb(spdep::knearneigh(sf::st_centroid(sf::st_geometry(partition), of_largest_polygon=TRUE), k=4))
                        )

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
