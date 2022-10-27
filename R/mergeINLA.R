#' Merge \code{inla} objects for partition models
#'
#' @description The function takes local models fitted for each subregion of the whole spatial domain and unifies them into a single \code{inla} object.
#' This function is valid for both disjoint and \emph{k}-order neighbourhood models.
#'
#' @details If the disjoint model is fitted (\code{k=0} argument), the log-risk surface is just the union of the posterior estimates of each submodel.
#' However, a single estimate of the overall log-risk can be computed by sampling from the joint posterior distribution of the linear predictors using the \code{inla.posterior.sample} function of R-INLA.
#' After joining the \eqn{S} samples from each submodel, the following global intercepts could be defined:
#' \itemize{
#' \item for spatial models fitted with \code{\link{CAR_INLA}} function \deqn{\tilde{\alpha}^s=\frac{1}{n}\sum_{i=1}^n (\log{r_{it}}-\mathbf{x_i}^{'}\beta), \quad \mbox{for} \quad s=1,\ldots,S}
#' \item for spatio-temporal models fitted with \code{\link{STCAR_INLA}} function \deqn{\tilde{\alpha}^s=\frac{1}{nT}\sum_{i=1}^n\sum_{t=1}^T \log{r_{it}}, \quad \mbox{for} \quad s=1,\ldots,S}
#' \item for multivariate spatial models fitted with \code{\link{MCAR_INLA}} function \deqn{\tilde{\alpha}_j^s=\frac{1}{n}\sum_{i=1}^n \log{r_{ij}}, \quad \mbox{for} \quad s=1,\ldots,S}
#' }
#' and then it is possible to compute the kernel density estimate of \eqn{\tilde{\alpha}}.
#' \cr \cr
#' If the \emph{k}-order neighbourhood model is fitted (\code{k>0} argument), note that the final risk surface \eqn{{\bf r}=(r_1,\ldots,r_{nT})^{'}} is no longer the union of the posterior estimates obtained from each submodel.
#' Since multiple log-risk estimates can be obtained for some areal-time units from the different local submodel, their posterior estimates must be properly combined to obtain a single posterior distribution for each \eqn{r_{it}}.
#' Two different merging strategies could be considered. If the \code{merge.strategy="mixture"} argument is specified, mixture distributions of the estimated posterior probability density functions with weights proportional to the conditional predictive ordinates (CPO) are computed.
#' If the \code{merge.strategy="original"} argument is specified (default option), the posterior marginal estimate ot the areal-unit corresponding to the original submodel is selected.
#' \cr \cr
#' See \insertCite{orozco2020;textual}{bigDM} and \insertCite{orozco2022;textual}{bigDM} for more details.
#'
#' @param inla.models list of multiple objects of class \code{inla}.
#' @param k numeric value with the neighbourhood order used for the partition model. If k=0 the \emph{Disjoint model} is considered.
#' @param ID.area character; name of the variable that contains the IDs of spatial areal units. Default to \code{"Area"}.
#' @param ID.year character; name of the variable that contains the IDs of time points. Default to \code{"NULL"} (for spatial models).
#' @param ID.disease character; name of the variable that contains the IDs of the diseases. Default to \code{"NULL"} (only required for multivariate models).
#' @param O character; name of the variable that contains the observed number of disease cases for each areal units. Default to \code{"O"}.
#' @param E character; name of the variable that contains either the expected number of disease cases or the population at risk for each areal unit. Default to \code{"E"}.
#' @param seed numeric; control the RNG of \code{inla.qsample} (see \code{help(inla.qsample)} for further information). Defaults to \code{NULL}.
#' @param n.sample numeric; number of samples to generate from the posterior marginal distribution of the risks. Default to 1000.
#' @param compute.intercept logical value (default \code{FALSE}); if \code{TRUE} then the overall log-risk \eqn{\alpha} (for spatial or spatio-temporal models fitted with either \code{\link{CAR_INLA}} or \code{\link{STCAR_INLA}} functions), or disease-specific overall log-risk \eqn{\alpha_j} (for multivariate models fitted with \code{\link{MCAR_INLA}} function) are computed.
#' Only works if \code{k=0} argument (\emph{Disjoint model}) is specified. CAUTION: This method might be very time consuming.
#' @param compute.DIC logical value; if \code{TRUE} (default) then approximate values of the Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) are computed.
#' @param merge.strategy one of either \code{"mixture"} or \code{"original"} (default), which specifies the merging strategy to compute posterior marginal estimates of relative risks.
#'
#' @return This function returns an object of class \code{inla} containing the following elements:
#' \item{\code{summary.fixed}}{A data.frame containing the mean, standard deviation, quantiles and mode of the model's fixed effects.}
#' \item{\code{marginals.fixed}}{A list containing the posterior marginal density of the model's fixed effects.}
#' \item{\code{summary.fixed.partition}}{A data.frame containing the mean, standard deviation, quantiles and mode of the model's fixed effects in each partition.}
#' \item{\code{marginals.fixed.partition}}{A list containing the posterior marginal density of the model's fixed effects in each partition.}
#' \item{\code{summary.random}}{If \code{k=0} a list with a data.frame containing the mean, standard deviation, quantiles and mode of the model's random effects.}
#' \item{\code{marginals.random}}{If \code{k=0} a list containing the posterior marginal densities of the model's random effects.}
#' \item{\code{summary.linear.predictor}}{If \code{k=0} a data.frame containing the mean, standard deviation, quantiles and mode of the log-risks (or log-rates) in the model.}
#' \item{\code{marginals.linear.predictor}}{If \code{k=0} a list containing the posterior marginal densities of the log-risks (or log-rates) in the model.}
#' \item{\code{summary.fitted.values}}{A data.frame containing the mean, standard deviation, quantiles, mode and cdf of the risks (or rates) in the model.}
#' \item{\code{marginals.fitted.values}}{A list containing the posterior marginal densities of the risks (or rates) in the model.}
#' \item{\code{summary.cor}}{A data.frame containing the mean, standard deviation, quantiles and mode of the between-disease correlation coefficients. Only for the multivariate spatial models fitted using the \code{\link{MCAR_INLA}} function.}
#' \item{\code{marginals.cor}}{A list containing the posterior marginal densities of the between-disease correlation coefficients. Only for the multivariate spatial models fitted using the \code{\link{MCAR_INLA}} function.}
#' \item{\code{summary.cor.partition}}{A data.frame containing the mean, standard deviation, quantiles and mode of the between-disease correlation coefficients in each partition. Only for the multivariate spatial models fitted using the \code{\link{MCAR_INLA}} function.}
#' \item{\code{marginals.cor.partition}}{A list containing the posterior marginal densities of the between-disease correlation coefficients in each partition. Only for the multivariate spatial models fitted using the \code{\link{MCAR_INLA}} function.}
#' \item{\code{summary.var}}{A data.frame containing the mean, standard deviation, quantiles and mode of the within-disease variances for each disease. Only for the multivariate spatial models fitted using the \code{\link{MCAR_INLA}} function.}
#' \item{\code{marginals.var}}{A list containing the posterior marginal densities of the within-disease variances for each disease. Only for the multivariate spatial models fitted using the \code{\link{MCAR_INLA}} function.}
#' \item{\code{summary.var.partition}}{A data.frame containing the mean, standard deviation, quantiles and mode of the within-disease variances in each partition. Only for the multivariate spatial models fitted using the \code{\link{MCAR_INLA}} function.}
#' \item{\code{marginals.var.partition}}{A list containing the posterior marginal densities of the within-disease variances in each partition. Only for the multivariate spatial models fitted using the \code{\link{MCAR_INLA}} function.}
#' \item{\code{logfile}}{A list of the log files of each submodel.}
#' \item{\code{version}}{A list containing information about the R-INLA version.}
#' \item{\code{cpu.used}}{The sum of cpu times used by the \code{inla} function for each submodel (\code{Pre}, \code{Running} and \code{Post}), and the cpu time of the merging process \code{Merging}.}
#'
#' @import parallel
#' @importFrom methods is
#' @importFrom stats runif density sd quantile dpois var
#' @importFrom rlist list.flatten
#'
#' @examples
#' ## See the vignettes accompanying this package for an example of its use.
#'
#' @export
mergeINLA <- function(inla.models=list(), k=NULL, ID.area="Area", ID.year=NULL, ID.disease=NULL, O="O", E="E",
                      seed=NULL, n.sample=1000, compute.intercept=FALSE, compute.DIC=TRUE, merge.strategy="original"){

  if(requireNamespace("INLA", quietly=TRUE)){

          D <- length(inla.models)

          ## Check for errors ##
          if(!is(inla.models,"list") | D==0)
                  stop("the 'inla.models' argument must be a non-empty list")
          if(unlist(unique(lapply(inla.models, function(x) class(x))))!="inla")
                  stop("the 'inla.models' argument must contain 'inla' class objects")
          if(is.null(k))
                  stop("the 'k' argument is missing")
          if(is.null(ID.area))
                  stop("the 'ID.area' argument is missing")
          if(is.null(O))
                  stop("the 'O' argument is missing")
          if(is.null(E))
                  stop("the 'E' argument is missing")
          if(!is.null(ID.year) & !is.null(ID.disease))
                  stop("both 'ID.year' and 'ID.disease' arguments cannot be non-null")
          if(!(merge.strategy %in% c("mixture","original")))
                  stop("invalid 'merge.strategy' argument")

          tt <- system.time({

                  if(is.null(seed)){
                          set.seed(format(Sys.Date(), "%Y%m%d"))
                  }else{
                          set.seed(seed)
                  }

                  result <- vector("list",51)
                  attr(result,"class") <- "inla"

                  names(result) <- c("names.fixed","summary.fixed","marginals.fixed","summary.lincomb","marginals.lincomb","size.lincomb",
                                     "summary.lincomb.derived","marginals.lincomb.derived","size.lincomb.derived","mlik","cpo","po","waic",
                                     "model.random","summary.random","marginals.random","size.random","summary.linear.predictor",
                                     "marginals.linear.predictor","summary.fitted.values","marginals.fitted.values","size.linear.predictor",
                                     "summary.hyperpar","marginals.hyperpar","internal.summary.hyperpar","internal.marginals.hyperpar",
                                     "offset.linear.predictor","model.spde2.blc","summary.spde2.blc","marginals.spde2.blc","size.spde2.blc",
                                     "model.spde3.blc","summary.spde3.blc","marginals.spde3.blc","size.spde3.blc","logfile","misc","dic",
                                     "mode","neffp","joint.hyper","nhyper","version","Q","graph","ok","cpu.used","all.hyper",".args",
                                     "call","model.matrix")

                  if(is.null(ID.year)){
                    if(is.null(ID.disease)){
                      ID.list <- lapply(inla.models, function(x) x$.args$data[,ID.area])
                      ID <- sort(as.character(unique(unlist(ID.list))))
                      
                      ID.group <- do.call(rbind,lapply(inla.models, function(x) x$.args$data[,c(ID.area,"ID.group")]))
                      ID.group <- ID.group[!duplicated(ID.group[,ID.area]),]
                      ID.group <- ID.group[order(ID.group[,ID.area]),]
                      rownames(ID.group) <- ID
                    }else{
                      ID.list <- lapply(inla.models, function(x) paste(x$.args$data[,ID.disease],x$.args$data[,ID.area],sep="."))
                      ID <- sort(as.character(unique(unlist(ID.list))))
                      
                      ID.group <- do.call(rbind,lapply(inla.models, function(x) x$.args$data[,c(ID.area,ID.disease,"ID.group")]))
                      ID.group <- ID.group[!duplicated(ID.group[,c(ID.area,ID.disease)]),]
                      ID.group <- ID.group[order(ID.group[,ID.disease],ID.group[,ID.area]),]
                      rownames(ID.group) <- ID
                    }
                          
                  }else{
                    ID.list <- lapply(inla.models, function(x) paste(x$.args$data[,ID.year],x$.args$data[,ID.area],sep="."))
                    ID <- sort(as.character(unique(unlist(ID.list))))
                    
                    ID.group <- do.call(rbind,lapply(inla.models, function(x) x$.args$data[,c(ID.area,ID.year,"ID.group")]))
                    ID.group <- ID.group[!duplicated(ID.group[,c(ID.area,ID.year)]),]
                    ID.group <- ID.group[order(ID.group[,ID.year],ID.group[,ID.area]),]
                    rownames(ID.group) <- ID
                  }

                  ## Fixed effects ##
                  names.fixed <- unique(lapply(inla.models, function(x) x$names.fixed))

                  if(length(names.fixed)>1){
                          stop("Different 'names.fixed' arguments for INLA models")
                  }else{
                    
                          result$names.fixed <- unlist(names.fixed)
                          
                          aux <- do.call(rbind,lapply(inla.models, function(x) x$summary.fixed))
                          rownames(aux) <- paste(rep(unlist(names.fixed),D),rep(formatC(1:D, width=ceiling(log(D+1,10)), flag='0'),each=length(unlist(names.fixed))),sep=".")
                          result$summary.fixed.partition <- aux[order(rownames(aux)),]

                          aux <- rlist::list.flatten(lapply(inla.models, function(x) x$marginals.fixed))
                          names(aux) <- paste(rep(unlist(names.fixed),D),rep(formatC(1:D, width=ceiling(log(D+1,10)), flag='0'),each=length(unlist(names.fixed))),sep=".")
                          result$marginals.fixed.partition <- aux[order(names(aux))]
                          
                          ## CMC algorithm for the fixed effects ##
                          fixed.CMC <- compute.CMC(marginals=result$marginals.fixed.partition, names=unlist(names.fixed))
                          result$summary.fixed <- fixed.CMC$summary.CMC
                          result$marginals.fixed <- fixed.CMC$marginals.CMC
                          
                          ## Compute model's intercept(s) by extracting samples from the joint posterior distribution of the linear predictors ##
                          if(k==0 & compute.intercept==TRUE){
                            
                            inla.seed <- as.integer(runif(1)*.Machine$integer.max)
                            
                            suppressWarnings({
                              cl <- makeCluster(detectCores())
                              clusterExport(cl, varlist=c("n.sample","inla.seed"), envir=environment())
                              clusterEvalQ(cl, INLA::inla.posterior.sample)
                              model.sample <- parLapply(cl, inla.models, modelSampling)
                              stopCluster(cl)
                            })
                           
                            if(is.null(ID.disease)){ # For univariate models
                              
                              # linear.predictor.sample <- do.call(rbind,lapply(model.sample, function(x) inla.posterior.sample.eval("Predictor",x)))
                              
                              X.names <- unlist(names.fixed)[-1]
                              X.beta <- lapply(inla.models, function(x) as.matrix(x$.args$data[,X.names])%*%x$summary.fixed[X.names,"mean"])
                              linear.predictor.sample <- lapply(model.sample, function(x) inla.posterior.sample.eval("Predictor",x))
                              linear.predictor.sample <- mapply(function(x,y) apply(x,2,function(xx) xx-y), x=linear.predictor.sample, y=X.beta, SIMPLIFY=FALSE)
                              linear.predictor.sample <- do.call(rbind,linear.predictor.sample)
                              
                              alpha.sample <- apply(linear.predictor.sample,2,function(x) mean(x))
                              alpha.density <- density(alpha.sample, n=75, bw="SJ")
                              
                              marginals.alpha <- cbind(x=alpha.density$x, y=alpha.density$y)
                              dimnames(marginals.alpha) <- list(NULL, c("x","y"))
                              
                              result$summary.intercept <- compute.summary(alpha.density)
                              rownames(result$summary.intercept) <- unlist(names.fixed)[1]
                              result$marginals.intercept <- list(marginals.alpha)
                              names(result$marginals.intercept) <- rownames(result$summary.intercept)

                            }else{ # For multivariate models

                              linear.predictor.sample <- lapply(model.sample, function(x) inla.posterior.sample.eval("Predictor",x))
                              
                              J <- length(unique(inla.models[[1]]$.args$data[,"ID.disease"]))
                              summary.intercept <- vector("list",J)
                              marginals.intercept <- vector("list",J)
                              names(summary.intercept) <- unlist(names.fixed)
                              names(marginals.intercept) <- unlist(names.fixed)
                              
                              for(j in seq(J)){
                                aux <- do.call(rbind,lapply(linear.predictor.sample, function(x) x[seq(nrow(x)/J*(j-1)+1,nrow(x)/J*j),]))
                                
                                alpha.sample <- apply(aux,2,function(x) mean(x))
                                alpha.density <- density(alpha.sample, n=75, bw="SJ")
                                
                                marginals.alpha <- cbind(x=alpha.density$x, y=alpha.density$y)
                                dimnames(marginals.alpha) <- list(NULL, c("x","y"))
                                
                                summary.intercept[[j]] <- compute.summary(alpha.density)
                                marginals.intercept[[j]] <- marginals.alpha
                              }
                              
                              result$summary.intercept <- do.call(rbind,summary.intercept)
                              result$marginals.intercept <- marginals.intercept

                            }
                          }
                  }

                  
                  ## lincomb / lincomb.derived ##
                  result$summary.lincomb <- data.frame()
                  # result$marginals.lincomb <- NULL
                  # result$size.lincomb <- NULL

                  result$summary.lincomb.derived <- data.frame()
                  # result$marginals.lincomb.derived <- NULL
                  # result$size.lincomb.derived <- NULL


                  ## mlik / cpo / po / waic ##
                  result$cpo <- list(cpo=logical(), pit=logical(), failure=logical())
                  result$po <- list(po=logical())
                  # result$mlik <- NULL
                  # result$waic <- NULL


                  ## Random effects ##
                  if(k==0){
                          model.random <- unique(lapply(inla.models, function(x) x$model.random))

                          if(length(model.random)>1){
                                  stop("Different 'model.random' arguments for INLA models")
                          }else{
                                  if(is.null(unlist(model.random))){
                                          # result$model.random <- NULL
                                          result$summary.random <- list()
                                          # result$marginals.random <- NULL
                                          # result$size.random <- NULL
                                  }else{

                                          result$model.random <- unlist(model.random)
                                          result$summary.random <- vector("list",length(result$model.random))
                                          names(result$summary.random) <- names(inla.models[[1]]$summary.random)
                                          
                                          if(is.null(ID.disease)){ # For univariate models
                                            
                                            result$summary.random$ID.area <- do.call(rbind,lapply(inla.models, function(x) x$summary.random$ID.area))
                                            result$summary.random$ID.area$ID <- as.character(unlist(lapply(inla.models, function(x) unique(x$.args$data[,ID.area]))))
                                            result$summary.random$ID.area <- result$summary.random$ID.area[order(result$summary.random$ID.area$ID),]
                                            result$summary.random$ID.area$ID <- seq(1,nrow(result$summary.random$ID.area))
                                            rownames(result$summary.random$ID.area) <- seq(1,nrow(result$summary.random$ID.area))
                                            
                                            if("ID.area.year" %in% names(result$summary.random)){
                                              result$summary.random$ID.area.year <- do.call(rbind,lapply(inla.models, function(x) x$summary.random$ID.area.year))
                                              result$summary.random$ID.area.year$ID <- as.character(unlist(lapply(inla.models, function(x) paste(x$.args$data[,ID.year],x$.args$data[,ID.area],sep="."))))
                                              result$summary.random$ID.area.year <- result$summary.random$ID.area.year[order(result$summary.random$ID.area.year$ID),]
                                              result$summary.random$ID.area.year$ID <- seq(1,nrow(result$summary.random$ID.area.year))
                                              rownames(result$summary.random$ID.area.year) <- seq(1,nrow(result$summary.random$ID.area.year))
                                            }
                                            
                                            result$marginals.random <- vector("list",length(unique(result$model.random)))
                                            names(result$marginals.random) <- names(inla.models[[1]]$summary.random)
                                            
                                            result$marginals.random$ID.area <- rlist::list.flatten(lapply(inla.models, function(x) x$marginals.random$ID.area))
                                            names(result$marginals.random$ID.area) <- as.character(unlist(lapply(inla.models, function(x) unique(x$.args$data[,ID.area]))))
                                            result$marginals.random$ID.area <- result$marginals.random$ID.area[order(names(result$marginals.random$ID.area))]
                                            names(result$marginals.random$ID.area) <- paste("index",seq(1:length(result$marginals.random$ID.area)),sep=".")
                                            
                                            if("ID.area.year" %in% names(result$summary.random)){
                                              result$marginals.random$ID.area.year <- rlist::list.flatten(lapply(inla.models, function(x) x$marginals.random$ID.area.year))
                                              names(result$marginals.random$ID.area.year) <- as.character(unlist(lapply(inla.models, function(x) paste(x$.args$data[,ID.year],x$.args$data[,ID.area],sep="."))))
                                              result$marginals.random$ID.area.year <- result$marginals.random$ID.area.year[order(names(result$marginals.random$ID.area.year))]
                                              names(result$marginals.random$ID.area.year) <- paste("index",seq(1:length(result$marginals.random$ID.area.year)),sep=".")
                                            }
                                            
                                            result$size.random <- vector("list",length(result$model.random))
                                            for(i in 1:length(result$size.random)){
                                              result$size.random[[i]] <- list(n=nrow(result$summary.random[[i]]), N=nrow(result$summary.random[[i]]), Ntotal=nrow(result$summary.random[[i]]), ngroup=1, nrep=1)
                                            }
                                            
                                          }else{ # For multivariate models
                                            result$summary.random$idx <- do.call(rbind,lapply(inla.models, function(x) x$summary.random$idx))
                                            result$summary.random$idx$ID <- as.character(unlist(lapply(inla.models, function(x) paste(x$.args$data[,ID.disease],x$.args$data[,ID.area],sep="."))))
                                            result$summary.random$idx <- result$summary.random$idx[order(result$summary.random$idx$ID),]
                                            result$summary.random$idx$ID <- seq(1,nrow(result$summary.random$idx))
                                            rownames(result$summary.random$idx) <- seq(1,nrow(result$summary.random$idx))
                                            
                                            result$marginals.random <- vector("list",length(unique(result$model.random)))
                                            names(result$marginals.random) <- names(inla.models[[1]]$summary.random)
                                            
                                            result$marginals.random$idx <- rlist::list.flatten(lapply(inla.models, function(x) x$marginals.random$idx))
                                            names(result$marginals.random$idx) <- as.character(unlist(lapply(inla.models, function(x) paste(x$.args$data[,ID.disease],x$.args$data[,ID.area],sep="."))))
                                            result$marginals.random$idx <- result$marginals.random$idx[order(names(result$marginals.random$idx))]
                                            names(result$marginals.random$idx) <- paste("index",seq(1:length(result$marginals.random$idx)),sep=".")
                                            
                                            result$size.random <- vector("list",length(result$model.random))
                                            for(i in 1:length(result$size.random)){
                                              result$size.random[[i]] <- list(n=nrow(result$summary.random[[i]]), N=nrow(result$summary.random[[i]]), Ntotal=nrow(result$summary.random[[i]]), ngroup=1, nrep=1)
                                            }
                                          }
                                  }
                          }
                  }else{
                          # result$model.random <- NULL
                          result$summary.random <- list()
                          # result$marginals.random <- NULL
                          # result$size.random <- NULL
                  }


                  ## Linear predictor ##
                  if(k==0){
                          summary.linear.predictor <- do.call(rbind,lapply(inla.models, function(x) x$summary.linear.predictor))
                          rownames(summary.linear.predictor) <- unlist(ID.list)
                          result$summary.linear.predictor <- summary.linear.predictor[order(rownames(summary.linear.predictor)),]
                          aux <- as.character(1:nrow(result$summary.linear.predictor))
                          l <- max(nchar(aux))
                          while(min(nchar(aux))<l) aux[nchar(aux)==min(nchar(aux))] <- paste("0",aux[nchar(aux)==min(nchar(aux))],sep="")
                          rownames(result$summary.linear.predictor) <- paste("Predictor",aux,sep=".")

                          marginals.linear.predictor <- rlist::list.flatten(lapply(inla.models, function(x) x$marginals.linear.predictor))
                          names(marginals.linear.predictor) <- unlist(ID.list)
                          result$marginals.linear.predictor <- marginals.linear.predictor[order(names(marginals.linear.predictor))]
                          names(result$marginals.linear.predictor) <- paste("Predictor",1:length(result$marginals.linear.predictor),sep=".")

                          result$size.linear.predictor <- list(n=nrow(result$summary.linear.predictor), N=nrow(result$summary.linear.predictor), Ntotal=nrow(result$summary.linear.predictor), ngroup=1, nrep=1)
                          result$offset.linear.predictor <- rep(0,nrow(result$summary.linear.predictor))
                  }else{
                          result$summary.linear.predictor <- data.frame()
                          result$marginals.linear.predictor <- vector("list",0)
                  }


                  ## Fitted values ##
                  result$summary.fitted.values <- data.frame()
                  result$marginals.fitted.values <- vector("list",0)

                  if(k==0){
                          summary.fitted.values <- do.call(rbind,lapply(inla.models, function(x) x$summary.fitted.values))
                          rownames(summary.fitted.values) <- unlist(ID.list)
                          result$summary.fitted.values <- summary.fitted.values[order(rownames(summary.fitted.values)),]
                          aux <- as.character(1:nrow(result$summary.fitted.values))
                          l <- max(nchar(aux))
                          while(min(nchar(aux))<l) aux[nchar(aux)==min(nchar(aux))] <- paste("0",aux[nchar(aux)==min(nchar(aux))],sep="")
                          rownames(result$summary.fitted.values) <- paste("fitted.Predictor",aux,sep=".")

                          marginals.fitted.values <- rlist::list.flatten(lapply(inla.models, function(x) x$marginals.fitted.values))
                          names(marginals.fitted.values) <- unlist(ID.list)
                          result$marginals.fitted.values <- marginals.fitted.values[order(names(marginals.fitted.values))]
                          names(result$marginals.fitted.values) <- paste("fitted.Predictor",aux,sep=".")
                  }else{
                          models.summary.fitted.values <- lapply(inla.models, function(x) x$summary.fitted.values)
                          models.marginals.fitted.values <- lapply(inla.models, function(x) x$marginals.fitted.values)
                          models.cpo <- lapply(inla.models, function(x) x$cpo)

                          suppressWarnings({
                                  cl <- makeCluster(detectCores())
                                  clusterExport(cl, varlist=c("models.summary.fitted.values","models.marginals.fitted.values","models.cpo","ID","ID.list","merge.strategy","ID.group"), envir=environment())
                                  clusterEvalQ(cl,{
                                          INLA::inla.dmarginal
                                          INLA::inla.qmarginal
                                          INLA::inla.rmarginal
                                          INLA::inla.mmarginal
                                  }
                                  )
                                  fitted.values.summary.marginal <- parLapply(cl,ID,computeFittedValues)
                                  stopCluster(cl)
                          })

                          result$summary.fitted.values <- do.call(rbind,lapply(fitted.values.summary.marginal, function(x) x[[1]]))
                          result$marginals.fitted.values <- lapply(fitted.values.summary.marginal, function(x) matrix(unlist(x[[2]]),ncol = 2))

                          aux <- as.character(1:nrow(result$summary.fitted.values))
                          l <- max(nchar(aux))
                          while(min(nchar(aux))<l) aux[nchar(aux)==min(nchar(aux))] <- paste("0",aux[nchar(aux)==min(nchar(aux))],sep="")
                          rownames(result$summary.fitted.values) <- paste("fitted.Predictor",aux,sep=".")

                          names(result$marginals.fitted.values) <-  rownames(result$summary.fitted.values)
                          result$marginals.fitted.values <- result$marginals.fitted.values[order(names(result$marginals.fitted.values))]
                  }


                  ## Hyperparameters ##
                  result$summary.hyperpar <- data.frame()
                  result$internal.summary.hyperpar <- data.frame()
                  # result$marginals.hyperpar <- NULL
                  # result$internal.marginals.hyperpar <- NULL


                  ## SPDE ##
                  result$summary.spde2.blc <- list()
                  result$summary.spde3.blc <- list()
                  # result$model.spde2.blc <- NULL
                  # result$marginals.spde2.blc <- NULL
                  # result$size.spde2.blc <- NULL
                  # result$model.spde3.blc <- NULL
                  # result$marginals.spde3.blc <- NULL
                  # result$size.spde3.blc <- NULL


                  ## logfile ##
                  result$logfile <- lapply(inla.models,function(x) x$logfile)
                  names(result$logfile) <- paste("Submodel",1:length(result$logfile),sep=".")


                  ## misc ##
                  result$misc <- list()

                  ## .args ##
                  result$.args <- vector("list",25)
                  names(result$.args) <- c("formula","family","data","quantiles","E","offset","verbose","control.compute","control.predictor",
                                           "control.family","control.inla","control.results","control.fixed","control.mode","control.expert",
                                           "control.lincomb","control.update","only.hyperparam","inla.call","num.threads","blas.num.threads",
                                           "keep","silent","debug",".parent.frame")

                  formula <- unique(lapply(inla.models, function(x) x$.args$formula))
                  if(length(formula)==1) result$.args$formula <- formula[[1]]

                  family <- unique(lapply(inla.models, function(x) x$.args$family))
                  if(length(family)==1) result$.args$family <- family[[1]]

                  result$.args$data <- do.call(rbind,lapply(inla.models, function(x) x$.args$data))
                  if(is.null(ID.year)){
                    if(is.null(ID.disease)){
                          result$.args$data$ID <- as.character(result$.args$data[,ID.area])
                          result$.args$data <- result$.args$data[!duplicated(result$.args$data$ID),]
                          result$.args$data <- result$.args$data[order(result$.args$data$ID),]
                    }else{
                      result$.args$data$ID <- paste(result$.args$data[,ID.disease],result$.args$data[,ID.area],sep=".")
                      result$.args$data <- result$.args$data[!duplicated(result$.args$data$ID),]
                      result$.args$data <- result$.args$data[order(result$.args$data$ID),]
                      result$.args$data$ID.area <- rep(1:length(unique(result$.args$data[,ID.area])), length(unique(result$.args$data[,ID.disease])))
                      result$.args$data$ID.disease <- rep(1:length(unique(result$.args$data[,ID.disease])), each=length(unique(result$.args$data[,ID.area])))
                      result$.args$data$idx <- seq(1:nrow(result$.args$data))
                    }
                  }else{
                          result$.args$data$ID <- paste(result$.args$data[,ID.year],result$.args$data[,ID.area],sep=".")
                          result$.args$data <- result$.args$data[!duplicated(result$.args$data$ID),]
                          result$.args$data <- result$.args$data[order(result$.args$data$ID),]
                          result$.args$data$ID.area <- rep(1:length(unique(result$.args$data[,ID.area])), length(unique(result$.args$data[,ID.year])))
                          result$.args$data$ID.year <- rep(1:length(unique(result$.args$data[,ID.year])), each=length(unique(result$.args$data[,ID.area])))
                          result$.args$data$ID.area.year <- seq(1:nrow(result$.args$data))
                  }
                  result$.args$data$ID <- NULL
                  rownames(result$.args$data) <- seq(1:nrow(result$.args$data))

                  quantiles <- unique(lapply(inla.models, function(x) x$.args$quantiles))
                  if(length(quantiles)==1) result$.args$quantiles <- quantiles[[1]]

                  result$.args$E <- result$.args$data[,E]

                  if(all(unlist(lapply(inla.models, function(x) x$.args$offset))==0)) result$.args$offset <- numeric(nrow(result$.args$data))

                  result$.args$verbose <- all(unlist(lapply(inla.models, function(x) x$.args$verbose)))

                  control.compute <- unique(lapply(inla.models, function(x) x$.args$control.compute))
                  if(length(control.compute)==1) result$.args$control.compute <- control.compute[[1]]

                  control.predictor <- unique(lapply(inla.models, function(x) x$.args$control.predictor))
                  if(length(control.predictor)==1) result$.args$control.predictor <- control.predictor[[1]]

                  control.family <- unique(lapply(inla.models, function(x) x$.args$control.family))
                  if(length(control.family)==1) result$.args$control.family <- control.family[[1]]

                  control.inla <- unique(lapply(inla.models, function(x) x$.args$control.inla))
                  if(length(control.inla)==1) result$.args$control.inla <- control.inla[[1]]

                  control.results <- unique(lapply(inla.models, function(x) x$.args$control.results))
                  if(length(control.results)==1) result$.args$control.results <- control.results[[1]]

                  control.fixed <- unique(lapply(inla.models, function(x) x$.args$control.fixed))
                  if(length(control.fixed)==1) result$.args$control.fixed <- control.fixed[[1]]

                  control.mode <- unique(lapply(inla.models, function(x) x$.args$control.mode))
                  if(length(control.mode)==1) result$.args$control.mode <- control.mode[[1]]

                  control.expert <- unique(lapply(inla.models, function(x) x$.args$control.expert))
                  if(length(control.expert)==1) result$.args$control.expert <- control.expert[[1]]

                  control.lincomb <- unique(lapply(inla.models, function(x) x$.args$control.lincomb))
                  if(length(control.lincomb)==1) result$.args$control.lincomb <- control.lincomb[[1]]

                  control.update <- unique(lapply(inla.models, function(x) x$.args$control.update))
                  if(length(control.update)==1) result$.args$control.update <- control.update[[1]]

                  result$.args$only.hyperparam <- all(unlist(lapply(inla.models, function(x) x$.args$only.hyperparam)))

                  inla.call <- unique(lapply(inla.models, function(x) x$.args$inla.call))
                  if(length(inla.call)==1) result$.args$inla.call <- inla.call[[1]]

                  num.threads <- unique(lapply(inla.models, function(x) x$.args$num.threads))
                  if(length(num.threads)==1) result$.args$num.threads <- num.threads[[1]]

                  blas.num.threads <- unique(lapply(inla.models, function(x) x$.args$blas.num.threads))
                  if(length(blas.num.threads)==1) result$.args$blas.num.threads <- blas.num.threads[[1]]

                  result$.args$keep <- all(unlist(lapply(inla.models, function(x) x$.args$keep)))

                  result$.args$silent <- all(unlist(lapply(inla.models, function(x) x$.args$silent)))

                  result$.args$debug <- all(unlist(lapply(inla.models, function(x) x$.args$debug)))

                  # parent.frame <- unique(lapply(inla.models, function(x) x$.args$.parent.frame))
                  # if(length(parent.frame)==1) result$.args$.parent.frame <- parent.frame[[1]]


                  ## Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) ##
                  if(compute.DIC){
                          marginals.fitted.values <- result$marginals.fitted.values

                          suppressWarnings({
                                  cl <- makeCluster(detectCores())
                                  clusterExport(cl, varlist=c("n.sample","marginals.fitted.values"), envir=environment())
                                  clusterEvalQ(cl, INLA::inla.rmarginal)
                                  risk.sample <- do.call(rbind,parLapply(cl,marginals.fitted.values,riskSampleDeviance))
                                  stopCluster(cl)
                          })

                          mu.sample <- apply(risk.sample, 2, function(x) result$.args$data[,E]*x)
                          result$dic$mean.deviance <- mean(apply(mu.sample, 2, function(x) -2*sum(log(dpois(result$.args$data[,O],x)))))
                          result$dic$deviance.mean <- -2*sum(log(dpois(result$.args$data[,O],apply(mu.sample,1,mean))))
                          result$dic$dic <- 2*result$dic$mean.deviance-result$dic$deviance.mean
                          result$dic$p.eff <- result$dic$mean.deviance-result$dic$deviance.mean

                          satured.deviance <- -2*sum(log(dpois(result$.args$data[,O],result$.args$data[,O])))
                          result$dic$mean.deviance.sat <- result$dic$mean.deviance-satured.deviance
                          result$dic$deviance.mean.sat <- result$dic$deviance.mean-satured.deviance
                          result$dic$dic.sat <- 2*result$dic$mean.deviance.sat-result$dic$deviance.mean.sat

                          lppd <- sum(log(apply(apply(mu.sample, 2, function(x) dpois(result$.args$data[,O],x)),1,mean)))
                          p.eff <- sum(apply(apply(mu.sample, 2, function(x) log(dpois(result$.args$data[,O],x))),1,var))
                          result$waic$waic <- -2*lppd+2*p.eff
                          result$waic$p.eff <- p.eff
                  }


                  ## mode / neffp / joint.hyper ##
                  result$mode <- list()
                  result$neffp <- data.frame()
                  # result$joint.hyper <- NULL


                  ## nhyper ##
                  nhyper <- unique(lapply(inla.models, function(x) x$nhyper))
                  if(length(nhyper)>1){
                          stop("Different 'nhyper' arguments for INLA models")
                  }else{
                          result$nhyper <- unlist(nhyper)
                  }


                  ## version ##
                  version <- unique(lapply(inla.models, function(x) x$version))
                  if(length(version)>1){
                          stop("Different 'version' arguments for INLA models")
                  }else{
                          result$version <- version[[1]]
                  }

                  ## Q / graph / ok ##
                  # result$Q <- NULL
                  # result$graph <- NULL
                  result$ok <- all(unlist(lapply(inla.models, function(x) x$ok)))


                  ## cpu.used ##
                  cpu.used <- lapply(inla.models, function(x) x$cpu.used)
                  result$cpu.used <- apply(do.call(rbind,cpu.used),2,sum)


                  ## all.hyper ##
                  all.hyper <- unique(lapply(inla.models, function(x) x$all.hyper))
                  if(length(all.hyper)==1) result$all.hyper <- aux[[1]]


                  ## call ##
                  call <- unique(lapply(inla.models, function(x) x$call))
                  if(length(call)==1) result$call <- call[[1]]


                  ## model.matrix ##
                  # result$model.matrix <- NULL
                  
                  ## Correlation coefficients between diseases and variances (only for multivariate models) ##
                  if(!is.null(ID.disease)){
                    names.cor <- unique(lapply(inla.models, function(x) rownames(x$summary.cor)))
                    names.var <- unique(lapply(inla.models, function(x) rownames(x$summary.var)))
                    
                    if(length(names.cor)>1){
                      stop("Different correlation coefficients for INLA models")
                    }else{
                      aux <- do.call(rbind,lapply(inla.models, function(x) x$summary.cor))
                      rownames(aux) <- paste(rep(unlist(names.cor),D),rep(formatC(1:D, width=ceiling(log(D+1,10)), flag='0'),each=length(unlist(names.cor))),sep=".")
                      result$summary.cor.partition <- aux[order(rownames(aux)),]
                      
                      aux <- rlist::list.flatten(lapply(inla.models, function(x) x$marginals.cor))
                      names(aux) <- paste(rep(unlist(names.cor),D),rep(formatC(1:D, width=ceiling(log(D+1,10)), flag='0'),each=length(unlist(names.cor))),sep=".")
                      result$marginals.cor.partition <- aux[order(names(aux))]
                      
                      ## CMC algorithm for the correlation coefficients ##
                      cor.CMC <- compute.CMC(marginals=result$marginals.cor.partition, names=unlist(names.cor))
                      result$summary.cor <- cor.CMC$summary.CMC
                      result$marginals.cor <- cor.CMC$marginals.CMC
                    }
                    
                    if(length(names.var)>1){
                      stop("Different variances for INLA models")
                    }else{
                      aux <- do.call(rbind,lapply(inla.models, function(x) x$summary.var))
                      rownames(aux) <- paste(rep(unlist(names.var),D),rep(formatC(1:D, width=ceiling(log(D+1,10)), flag='0'),each=length(unlist(names.var))),sep=".")
                      result$summary.var.partition <- aux[order(rownames(aux)),]
                      
                      aux <- rlist::list.flatten(lapply(inla.models, function(x) x$marginals.var))
                      names(aux) <- paste(rep(unlist(names.var),D),rep(formatC(1:D, width=ceiling(log(D+1,10)), flag='0'),each=length(unlist(names.var))),sep=".")
                      result$marginals.var.partition <- aux[order(names(aux))]
                      
                      ## CMC algorithm for the correlation coefficients ##
                      var.CMC <- compute.CMC(marginals=result$marginals.var.partition, names=unlist(names.var))
                      result$summary.var <- var.CMC$summary.CMC
                      result$marginals.var <- var.CMC$marginals.CMC
                    }
                  }
          })

          result$cpu.used <- c(result$cpu.used[1:3], Merging=as.numeric(tt[3]), Total=as.numeric(result$cpu.used[4]+tt[3]))

          return(result)
  }else{
          stop("INLA library is not installed! Please use following command to install the stable version of the R-INLA package:\n install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }
}


#########################
## Auxiliary functions ##
#########################
modelSampling <- function(x){
  INLA::inla.posterior.sample(n.sample, x, seed=inla.seed)
}

computeFittedValues <- function(q){

  pos <- lapply(ID.list, function(x) which(x==q))
  i <- which(unlist(lapply(pos, function(x) length(x)))>0)

  if(length(i)==1){
    summary.fitted.values <- models.summary.fitted.values[[i]][pos[[i]],1:7]
    rownames(summary.fitted.values) <- q

    marginals.fitted.values <- models.marginals.fitted.values[[i]][pos[[i]]]
    names(marginals.fitted.values) <- q
  }else{
    if(merge.strategy=="original"){
      aux <- ID.group[which(rownames(ID.group)==q),"ID.group"]
      
      summary.fitted.values <- models.summary.fitted.values[[aux]][pos[[aux]],1:7]
      rownames(summary.fitted.values) <- q
      
      marginals.fitted.values <- models.marginals.fitted.values[[aux]][pos[[aux]]]
      names(marginals.fitted.values) <- q
    }else{
      post.mean <- mapply(function(x,y){x$mean[y]}, x=models.summary.fitted.values[i], y=pos[i])
      post.sd <- mapply(function(x,y){x$sd[y]}, x=models.summary.fitted.values[i], y=pos[i])
      post.cdf <- mapply(function(x,y){x$'1 cdf'[y]}, x=models.summary.fitted.values[i], y=pos[i])
      marginals <- mapply(function(x,y){x[[y]]}, x=models.marginals.fitted.values[i], y=pos[i], SIMPLIFY=FALSE)

      cpo <- mapply(function(x,y){x$cpo[y]}, x=models.cpo[i], y=pos[i])
      w <- cpo/sum(cpo)
      
      xx <- sort(unlist(lapply(marginals, function(x) x[,"x"])))
      at <- round(seq(1,length(xx),length.out=75))
      marginals.fitted.values <- matrix(0, nrow=0, ncol=2, dimnames=list(NULL, c("x","y")))
      for(j in xx[at]){
        aux <- unlist(lapply(marginals, function(x) INLA::inla.dmarginal(j,x)))
        marginals.fitted.values <- rbind(marginals.fitted.values,c(j,sum(aux*w)))
      }
      
      summary.fitted.values <- data.frame(sum(post.mean*w),
                                          sqrt(sum((post.sd^2+post.mean^2-sum(post.mean*w)^2)*w)),
                                          INLA::inla.qmarginal(0.025,marginals.fitted.values),
                                          INLA::inla.qmarginal(0.5,marginals.fitted.values),
                                          INLA::inla.qmarginal(0.975,marginals.fitted.values),
                                          INLA::inla.mmarginal(marginals.fitted.values),
                                          sum(post.cdf*w))
      colnames(summary.fitted.values) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode","1 cdf")
      rownames(summary.fitted.values) <- q
      
      marginals.fitted.values <- list(marginals.fitted.values)
      names(marginals.fitted.values) <- q
    }
  }
  return(list(summary.fitted.values,marginals.fitted.values))
}

riskSampleDeviance <- function(x){
  INLA::inla.rmarginal(n.sample,x)
}

compute.summary <- function(marginal){
  aux <- data.frame(INLA::inla.emarginal(function(x) x, marginal),
                    sqrt(INLA::inla.emarginal(function(x) x^2, marginal)-INLA::inla.emarginal(function(x) x, marginal)^2),
                    INLA::inla.qmarginal(0.025,marginal),
                    INLA::inla.qmarginal(0.5,marginal),
                    INLA::inla.qmarginal(0.975,marginal),
                    INLA::inla.mmarginal(marginal))
  colnames(aux) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode")
  aux
}

compute.CMC <- function(marginals,names){
  
  marginals.CMC <- vector("list",length(names))
  names(marginals.CMC) <- names
  
  for(i in names){
    pos <- grep(i,names(marginals))
    pos <- pos[!unlist(lapply(marginals[pos], function(x) all(is.na(x))))]
    
    w <- 1/unlist(lapply(marginals[pos], function(x) INLA::inla.emarginal(function(y) y^2,x)-INLA::inla.emarginal(function(y) y,x)^2))
    w <- w/sum(w)
    
    marginals.sample <- lapply(marginals[pos], function(x) INLA::inla.rmarginal(10000,x))
    marginals.density <- density(c(do.call(cbind,marginals.sample)%*%w), n=75, bw="SJ")
    marginals.matrix <- cbind(x=marginals.density$x, y=marginals.density$y)
    dimnames(marginals.matrix) <- list(NULL, c("x","y"))
    
    marginals.CMC[[i]] <- marginals.matrix
  }
  
  summary.CMC <- do.call(rbind,lapply(marginals.CMC, function(x) compute.summary(x)))
  summary.CMC <- summary.CMC[,c("mean","sd","0.025quant","0.5quant","0.975quant","mode")]
  
  return(list(marginals.CMC=marginals.CMC, summary.CMC=summary.CMC))
}

utils::globalVariables(c("ID.list",
                         "models.summary.fitted.values",
                         "models.marginals.fitted.values",
                         "models.cpo",
                         "n.sample",
                         "inla.seed",
                         "inla.posterior.sample",
                         "inla.posterior.sample.eval",
                         "inla.rmarginal",
                         "inla.dmarginal",
                         "inla.qmarginal",
                         "inla.mmarginal",
                         "inla.make.lincombs",
                         "inla",
                         "merge.strategy",
                         "ID.group"))