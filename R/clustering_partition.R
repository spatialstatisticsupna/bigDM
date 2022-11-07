#' Obtain a partition of the spatial domain using the density-based spatial clustering (DBSC) algorithm described in Santafé et al. (2021)
#'
#' @description The function takes an object of class \code{SpatialPolygonsDataFrame} or \code{sf} and defines a spatial partition using the DBSC algorithm described in \insertCite{santafe2021;textual}{bigDM}.
#'
#' @details The DBSC algorithm implemented in this function is a new spatial clustering algorithm based on the density clustering algorithm introduced by \insertCite{rodriguez2014clustering;textual}{bigDM} and the posterior modification presented by \insertCite{wang2016automatic;textual}{bigDM}.
#' This algorithm is able to obtain a single clustering partition of the data by automatically detecting clustering centers and assigning each area to its nearest cluster centroid.
#' The algorithm has its basis in the assumption that cluster centers are points with high local density and relatively large distance to other points with higher local densities.
#' See \insertCite{santafe2021;textual}{bigDM} for more details.
#'
#' @references
#' \insertRef{rodriguez2014clustering}{bigDM}
#'
#' \insertRef{santafe2021}{bigDM}
#'
#' \insertRef{wang2016automatic}{bigDM}
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
#' @param ID.area character; name of the variable that contains the IDs of spatial areal units.
#' @param var character; name of the variable that contains the data of interest to compute spatial clusters, usually the vector of log-SMR.
#' @param n.cluster numeric; value to fix the number of cluster centers in the DBSC algorithm. Default to 10.
#' @param min.size numeric (default \code{NULL}); value to fix the minimum size of areas in each spatial partition.
#' @param W optional argument with the binary adjacency matrix of the spatial areal units. If \code{NULL} (default), this object is computed from the \code{carto} argument (two areas are considered as neighbours if they share a common border).
#' @param l numeric value with the neighbourhood order used to assign areas to each cluster. If \code{k=1} (default), only areas that share a common border are considered.
#' @param Wk previously computed binary adjacency matrix of l-order neighbours. If this argument is included (default \code{NULL}), the parameter \code{l} is ignored.
#' @param distance the distance measure to be used (default \code{"euclidean"}). See the \code{method} argument of \code{\link[stats]{dist}} function for other options.
#' @param verbose logical value (default \code{TRUE}); indicates if the function runs in verbose mode.
#'
#' @return \code{sf} object with the original data and a grouping variable named 'ID.group'.
#'
#' @importFrom grDevices boxplot.stats
#' @importFrom Matrix colSums
#' @importFrom stats dist
#'
#' @examples
#'\dontrun{
#' library(foreign)
#' library(maptools)
#' library(rgdal)
#' library(tmap)
#'
#' ## Load the Spain colorectal cancer mortality data ##
#' data(Carto_SpainMUN)
#'
#' ## Define a spatial partition using the DBSC algorithm ##
#' Carto_SpainMUN$logSMR <- log(Carto_SpainMUN$obs/Carto_SpainMUN$exp+0.0001)
#'
#' carto.new <- clustering_partition(carto=Carto_SpainMUN, ID.area="ID", var="logSMR",
#'                                   n.cluster=20, l=2, min.size=100, verbose=TRUE)
#' table(carto.new$ID.group)
#'
#' ## Plot of the grouping variable 'ID.group' ##
#' carto.partition <- unionSpatialPolygons(as(carto.new,"Spatial"),carto.new$ID.group)
#'
#' tm_shape(carto.new) +
#'         tm_polygons(col="ID.group") +
#'         tm_shape(carto.partition) +
#'         tm_borders(col="black", lwd=2) +
#'         tm_layout(legend.outside=TRUE)
#'}
#'
#' @export
clustering_partition <- function(carto, ID.area=NULL, var=NULL, n.cluster=10, min.size=NULL, W=NULL, l=1, Wk=NULL, distance="euclidean", verbose=TRUE){

        ## Check for errors ##
        if(is.null(carto))
                stop("the carto argument is missing")
        if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
                stop("the carto argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
        if(!is.null(n.cluster)){
                if(n.cluster<=0 | n.cluster>nrow(carto))
                        stop("invalid value for 'n.cluster' argument: 1<=n.cluster<=nrow(carto)")
        }else{
                n.cluster <- -1
        }

        ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class ##
        carto <- sf::st_as_sf(carto)
        data <- sf::st_set_geometry(carto, NULL)
        n <- nrow(data)

        if(n.cluster==1){
                ## Return only 1 cluster ##
                carto$ID.group <- rep(1,n)
        }else{
                ## Compute the binary adjacency matrix of l-order neighbours ##
                if(verbose) cat("STEP 1: computing adjacency matrix...\n")

                if(is.null(W)){
                        invisible(utils::capture.output(aux <- connect_subgraphs(carto, ID.area)))
                        W <- aux$W
                }

                if(is.null(Wk)){
                        if(l>1){
                                Wk <- W
                                for(i in 2:l){
                                        Wk <- Wk%*%Wk
                                        diag(Wk) <- 0
                                        Wk[Wk>0] <- 1
                                }
                        }else{
                                Wk <- W
                        }
                }

                ## DBSC algorithm (Santafé et al., 2021) ##
                if(!var %in% colnames(data)){
                        stop(paste("no",var,"variable found in carto object"))
                }else{
                        if(verbose) cat("STEP 2: computing distance matrix...\n")
                        distanceMatrix <- as.matrix(stats::dist(data[,var], method=distance, diag=TRUE, upper=TRUE))
                        rho <- Matrix::colSums(Wk)/(Matrix::colSums(Wk*distanceMatrix) + 1e-10)
                        rho.order <- sort(rho,index.return=TRUE)

                        delta <- numeric(n)
                        for(index in seq(n-1)){
                                i <- rho.order$ix[index]
                                delta[i] <- min(distanceMatrix[i,rho.order$ix[(index+1):n]])
                        }
                        delta[rho.order$ix[n]] <- max(delta)

                        gamma <- rho*delta
                        outlierValues <- grDevices::boxplot.stats(gamma, coef=2)$out
                        if(n.cluster!=-1) outlierValues <- sort(outlierValues, decreasing=TRUE)[1:n.cluster]
                        centers <- which(gamma %in% outlierValues)
                        nClusters <- length(centers)

                        if(verbose) cat("STEP 3: Computing spatial partition...\n")
                        clusteringPartition <- numeric(n)
                        clusteringPartition[centers] <- seq(nClusters)

                        connected <- as.matrix(W[centers,])
                        connected[connected==0] <- Inf

                        clusterDistances <- distanceMatrix[centers,]
                        clusterDistances[,centers] <- Inf

                        used <- logical(n)
                        used[centers] <- TRUE

                        it <- 1
                        while(any(!used) & it<=n){
                                if(verbose) cat(sprintf(" + Iteration %d \n",it))
                                distanceToNeighbors <- clusterDistances*connected
                                index <- which.min(distanceToNeighbors)
                                clusterIndex <- index %% nClusters
                                objectIndex <- index %/% nClusters + 1
                                if(clusterIndex == 0){
                                        clusterIndex <- nClusters
                                        objectIndex <- objectIndex - 1
                                }
                                clusteringPartition[objectIndex] <- clusterIndex
                                connected[clusterIndex,W[objectIndex,]==1] <- 1
                                used[objectIndex] <- TRUE
                                clusterDistances[,objectIndex] <- Inf

                                it <- it+1
                        }

                        if(!is.null(min.size)){
                                if(verbose) cat(sprintf("STEP 4: merging small clusters (min.size=%d)...\n",min.size))

                                cluster.size <- sort(table(clusteringPartition))
                                D <- distanceMatrix[centers,centers]

                                while(any(cluster.size<min.size)){
                                        clusterIndex <- as.numeric(names(cluster.size)[1])
                                        objectIndex <- which(clusteringPartition==clusterIndex)

                                        aux <- unique(clusteringPartition[setdiff(which(connected[clusterIndex,]==1),objectIndex)])
                                        mergeCluster <- as.numeric(which(D[clusterIndex,]==min(D[clusterIndex,aux])))

                                        clusteringPartition[objectIndex] <- mergeCluster
                                        D[clusterIndex,] <- Inf
                                        D[,clusterIndex] <- Inf

                                        cluster.size <- sort(table(clusteringPartition))
                                        it <- it+1
                                }
                        }
                }

                carto$ID.group <- factor(as.numeric(factor(clusteringPartition)))
        }

        return(carto)
}

utils::globalVariables(c("ID"))
