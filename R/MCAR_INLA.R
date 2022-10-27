#' Fit a (scalable) spatial multivariate Poisson mixed model to areal count data where dependence between spatial patterns of the diseases is addressed through the use of M-models \insertCite{botella2015unifying}{bigDM}.
#' 
#' @description Fit a spatial multivariate Poisson mixed model to areal count data. The linear predictor is modelled as \deqn{\log{r_{ij}}=\alpha_j + \theta_{ij}, \quad \mbox{for} \quad i=1,\ldots,n; \quad j=1,\ldots,J}
#' where \eqn{\alpha_j} is a disease-specific intercept and \eqn{\theta_{ij}} is the spatial main effect of area \eqn{i} for the \eqn{j}-th disease.
#' Following \insertCite{botella2015unifying;textual}{bigDM}, we rearrange the spatial effects into the matrix \eqn{\mathbf{\Theta} = \{ \theta_{ij}: i=1, \ldots, I; j=1, \ldots, J \}} whose columns are spatial random effects and its joint distribution specifies how dependence within-diseases and between-diseases is defined.
#' Several conditional autoregressive (CAR) prior distributions can be specified to deal with spatial dependence within-diseases, such as the intrinsic CAR prior \insertCite{besag1991}{bigDM}, the CAR prior proposed by \insertCite{leroux1999estimation;textual}{bigDM}, and the proper CAR prior distribution.
#' \cr\cr
#' As in the \code{\link{CAR_INLA}} function, three main modelling approaches can be considered:
#' \itemize{
#' \item the usual model with a global spatial random effect whose dependence structure is based on the whole neighbourhood graph of the areal units (\code{model="global"} argument)
#' \item a Disjoint model based on a partition of the whole spatial domain where independent spatial CAR models are simultaneously fitted in each partition (\code{model="partition"} and \code{k=0} arguments)
#' \item a modelling approach where \emph{k}-order neighbours are added to each partition to avoid border effects in the Disjoint model (\code{model="partition"} and \code{k>0} arguments).
#' }
#' For both the Disjoint and k-order neighbour models, parallel or distributed computation strategies can be performed to speed up computations by using the 'future' package \insertCite{bengtsson2020unifying}{bigDM}.
#'
#' Inference is conducted in a fully Bayesian setting using the integrated nested Laplace approximation (INLA; \insertCite{rue2009approximate;textual}{bigDM}) technique through the R-INLA package (\url{https://www.r-inla.org/}).
#' For the scalable model proposals \insertCite{orozco2020}{bigDM}, approximate values of the Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) can also be computed.
#'
#' @details For a full model specification and further details see the vignettes accompanying this package.
#' 
#' @references
#' \insertRef{bengtsson2020unifying}{bigDM}
#'
#' \insertRef{besag1991}{bigDM}
#' 
#' \insertRef{botella2015unifying}{bigDM}
#'
#' \insertRef{leroux1999estimation}{bigDM}
#'
#' \insertRef{rue2009approximate}{bigDM}
#'
#' \insertRef{vicente2022}{bigDM}
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}. This object must contain at least the variable with the identifiers of the spatial areal units specified in the argument \code{ID.area}.
#' @param data object of class \code{data.frame} that must contain the target variables of interest specified in the arguments \code{ID.area}, \code{ID.disease}, \code{O} and \code{E}.
#' @param ID.area character; name of the variable that contains the IDs of spatial areal units. The values of this variable must match those given in the \code{carto} and \code{data} variable.
#' @param ID.disease character; name of the variable that contains the IDs of the diseases.
#' @param ID.group character; name of the variable that contains the IDs of the spatial partition (grouping variable). Only required if \code{model="partition"}.
#' @param O character; name of the variable that contains the observed number of cases for each areal unit and disease.
#' @param E character; name of the variable that contains either the expected number of cases or the population at risk for each areal unit and disease.
#' @param W optional argument with the binary adjacency matrix of the spatial areal units. If \code{NULL} (default), this object is computed from the \code{carto} argument (two areas are considered as neighbours if they share a common border).
#' @param prior one of either \code{"intrinsic"} (default), \code{"Leroux"}, \code{"proper"}, or \code{"iid"} which specifies the prior distribution considered for the spatial random effect.
#' @param model one of either \code{"global"} or \code{"partition"} (default), which specifies the \emph{Global model} or one of the scalable model proposal's (\emph{Disjoint model} and \emph{k-order neighbourhood model}, respectively).
#' @param k numeric value with the neighbourhood order used for the partition model. Usually k=2 or 3 is enough to get good results. If k=0 (default) the \emph{Disjoint model} is considered. Only required if \code{model="partition"}.
#' @param strategy one of either \code{"gaussian"}, \code{"simplified.laplace"} (default), \code{"laplace"} or \code{"adaptive"}, which specifies the approximation strategy considered in the \code{inla} function.
#' @param seed numeric (default \code{NULL}); control the RNG of the \code{inla.qsample} function. See \code{help(inla.qsample)} for further information.
#' @param n.sample numeric; number of samples to generate from the posterior marginal distribution of the risks. Default to 1000.
#' @param compute.intercept logical value (default \code{FALSE}); if \code{TRUE} then disease-specific overall log-risk \eqn{\alpha_{j}} are computed.
#' Only works if \code{k=0} argument (\emph{Disjoint model}) is specified. CAUTION: This method might be very time consuming.
#' @param compute.DIC logical value; if \code{TRUE} (default) then approximate values of the Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) are computed.
#' @param save.models logical value (default \code{FALSE}); if \code{TRUE} then a list with all the \code{inla} submodels is saved in '/temp/' folder, which can be used as input argument for the \code{\link{mergeINLA}} function.
#' @param plan one of either \code{"sequential"} or \code{"cluster"}, which specifies the computation strategy used for model fitting using the 'future' package.
#' If \code{plan="sequential"} (default) the models are fitted sequentially and in the current R session (local machine). If \code{plan="cluster"} the models are fitted in parallel on external R sessions (local machine) or distributed in remote compute nodes.
#' @param workers character or vector (default \code{NULL}) containing the identifications of the local or remote workers where the models are going to be processed. Only required if \code{plan="cluster"}.
#' @param merge.strategy one of either \code{"mixture"} or \code{"original"} (default), which specifies the merging strategy to compute posterior marginal estimates of relative risks. See \code{\link{mergeINLA}} for further details.
#' 
#' @return This function returns an object of class \code{inla}. See the \code{\link{mergeINLA}} function for details.
#' 
#' @import crayon future Matrix parallel Rdpack
#' @importFrom fastDummies dummy_cols
#' @importFrom future.apply future_mapply
#' @importFrom sf st_as_sf st_set_geometry
#' @importFrom stats as.formula cov
#' @importFrom utils capture.output combn
#' 
#' @examples
#' \dontrun{
#' if(require("INLA", quietly=TRUE)){
#' 
#'   ## Load the sf object that contains the spatial polygons of the municipalities of Spain ##
#'   data(Carto_SpainMUN)
#'   str(Carto_SpainMUN)
#'
#'   ## Load the simulated cancer mortality data (three diseases) ##
#'   data(Data_MultiCancer)
#'   str(Data_MultiCancer)
#'
#'   ## Fit the Global model with an iCAR prior for the within-disease random effects ##
#'   Global <- MCAR_INLA(carto=Carto_SpainMUN, data=Data_MultiCancer,
#'                       ID.area="ID", ID.disease="disease", O="obs", E="exp",
#'                       prior="intrinsic", model="global", strategy="gaussian")
#'   summary(Global)
#'   
#'   ## Fit the Disjoint model with an iCAR prior for the within-disease random effects ##
#'   ## using 4 local clusters to fit the models in parallel                            ##
#'   Disjoint <- MCAR_INLA(carto=Carto_SpainMUN, data=Data_MultiCancer,
#'                         ID.area="ID", ID.disease="disease", O="obs", E="exp", ID.group="region", 
#'                         prior="intrinsic", model="partition", k=0, strategy="gaussian",
#'                         plan="cluster", workers=rep("localhost",4))
#'   summary(Disjoint)
#'
#'   ## 1st-order neighbourhood model with an iCAR prior for the within-disease random effects ##
#'   ## using 4 local clusters to fit the models in parallel                                   ##
#'   order1 <- MCAR_INLA(carto=Carto_SpainMUN, data=Data_MultiCancer,
#'                       ID.area="ID", ID.disease="disease", O="obs", E="exp", ID.group="region",
#'                       prior="intrinsic", model="partition", k=1, strategy="gaussian",
#'                       plan="cluster", workers=rep("localhost",4))
#'   summary(order1)
#' }
#' }
#'
#' @export
MCAR_INLA <- function(carto=NULL, data=NULL, ID.area=NULL, ID.disease=NULL, ID.group=NULL, O=NULL, E=NULL,
                      W=NULL, prior="intrinsic", model="partition", k=0, strategy="simplified.laplace",
                      seed=NULL, n.sample=1000, compute.intercept=FALSE, compute.DIC=TRUE,
                      save.models=FALSE, plan="sequential", workers=NULL, merge.strategy="original"){
  
  if(suppressPackageStartupMessages(requireNamespace("INLA", quietly=TRUE))){
    
    ## Check for errors ##
    if(is.null(carto))
      stop("the 'carto' argument is missing")
    if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
      stop("the 'carto' argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
    if(is.null(ID.area))
      stop("the 'ID.area' argument is missing")
    if(is.null(ID.disease))
      stop("the 'ID.disease' argument is missing")
    if(is.null(O))
      stop("the 'O' argument is missing")
    if(is.null(E))
      stop("the 'E' argument is missing")
    if(!(prior %in% c("intrinsic","Leroux","proper","iid")))
      stop("invalid 'prior' argument")
    if(!(model %in% c("global","partition")))
      stop("invalid 'model' argument")
    if(!(strategy %in% c("gaussian","simplified.laplace","laplace","adaptative")))
      stop("invalid 'strategy' argument")
    if(!(plan %in% c("sequential","cluster")))
      stop("invalid 'plan' argument")
    if(plan=="cluster" & is.null(workers))
      stop("argument 'workers' must be specified when using plan='cluster' computation strategy")
    if(!(merge.strategy %in% c("mixture","original")))
      stop("invalid 'merge.strategy' argument")
    
    cat("STEP 1: Pre-processing data\n")
    
    ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class
    carto <- sf::st_as_sf(carto)
    
    ## Order the data ##
    if(!ID.area %in% colnames(carto))
      stop(sprintf("'%s' variable not found in carto object",ID.area))
    if(!ID.area %in% colnames(data))
      stop(sprintf("'%s' variable not found in data object",ID.area))
    if(!ID.disease %in% colnames(data))
      stop(sprintf("'%s' variable not found in data object",ID.disease))
    if(!O %in% colnames(data))
      stop(sprintf("'%s' variable not found in carto object",O))
    if(!E %in% colnames(data))
      stop(sprintf("'%s' variable not found in carto object",E))
    
    carto <- carto[order(sf::st_set_geometry(carto, NULL)[,ID.area]),]
    data <- merge(data,carto[,c(ID.area,ID.group)])
    data$geometry <- NULL
    data[,ID.disease] <- paste(sprintf("%02d", as.numeric(as.character(data[,ID.disease]))))
    data <- data[order(data[,ID.disease],data[,ID.area]),]
    
    ## Merge disjoint connected subgraphs ##
    if(is.null(W)){
      invisible(utils::capture.output(aux <- connect_subgraphs(carto, ID.area)))
      carto.nb <- aux$nb
    }else{
      carto.nb <- spdep::mat2listw(W, style="B")$neighbours
      invisible(utils::capture.output(aux <- connect_subgraphs(carto, ID.area, nb=carto.nb)))
      carto.nb <- aux$nb
    }
    
    ## Define hyperprior distributions ##
    sdunif="expression:
          logdens=-log_precision/2;
          return(logdens)"
    
    lunif = "expression:
          a = 1;
          b = 1;
          beta = exp(theta)/(1+exp(theta));
          logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
          log_jacobian = log(beta*(1-beta));
          return(logdens+log_jacobian)"
    
    ## Formula for INLA model ##
    J <- length(unique(data[,ID.disease]))
    
    form <- "O ~ -1+"
    form <- paste(form, paste(paste0("I",1:J),collapse="+"),sep="")
    form <- paste(form, "+ f(idx, model=Mmodel, constr=FALSE, extraconstr=list(A=A.constr, e=rep(0,J)))")
    formula <- stats::as.formula(form)

    ## Auxiliary functions ##
    FitModels <- function(W, A.constr, data.INLA, d, D, initial.values){
      
      cat(sprintf("+ Model %d of %d",d,D),"\n")
      
      form <- "O ~ -1+"
      form <- paste(form, paste(paste0("I",1:J),collapse="+"),sep="")
      form <- paste(form, "+ f(idx, model=Mmodel, constr=FALSE, extraconstr=list(A=A.constr, e=rep(0,J)))")
      formula <- stats::as.formula(form)
      
      if(prior=="intrinsic"){
        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar, debug=FALSE, J=J, W=W, initial.values=initial.values)
      }
      if(prior=="Leroux"){
        Mmodel <- INLA::inla.rgeneric.define(Mmodel_lcar, debug=FALSE, J=J, W=W, initial.values=initial.values, alpha.min=0, alpha.max=1)
      }
      if(prior=="proper"){
        Mmodel <- INLA::inla.rgeneric.define(Mmodel_pcar, debug=FALSE, J=J, W=W, initial.values=initial.values, alpha.min=0, alpha.max=1)
      }
      if(prior=="iid"){
        Mmodel <- INLA::inla.rgeneric.define(Mmodel_iid, debug=FALSE, J=J, W=W, initial.values=initial.values)
      }

      models <- inla(formula, family="poisson", data=data.INLA, E=E,
                     control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                     control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                     control.inla=list(strategy=strategy))
      
      models$Mmodel <- list(model=model, prior=prior)
      
      Mmodel.compute <- Mmodel_compute_cor(models, n.sample=1000)
      models$summary.cor <- Mmodel.compute$summary.cor
      models$marginals.cor <- Mmodel.compute$marginals.cor
      models$summary.var <- Mmodel.compute$summary.var
      models$marginals.var <- Mmodel.compute$marginals.var
      
      return(models)
    }
    
    ## Global model ##
    if(model=="global"){
      cat("STEP 2: Fitting global model with INLA (this may take a while...)\n")
      
      W <- aux$W
      S <- nrow(W)

      A.constr <- kronecker(diag(J), matrix(1,1,S))
      
      data.INLA <- data.frame(O=data[,O], E=data[,E], Area=data[,ID.area], Disease=data[,ID.disease],
                              ID.area=rep(1:S,J), ID.disease=rep(1:J,each=S), idx=seq(1,J*S))
      
      intercepts <- fastDummies::dummy_cols(data.INLA$ID.disease)[,-1]
      intercepts[intercepts==0] <- NA
      colnames(intercepts) <- paste0("I",1:J)
      data.INLA <- cbind(data.INLA, intercepts)
      
      # aux <- log(data.INLA[,O]/data.INLA[,E])
      # aux[aux==-Inf] <- NA
      # aux[is.na(aux)] <- min(aux,na.rm=T)
      # Sigma <- cov(matrix(aux,S,J,byrow=F))
      Sigma <- cov(matrix(data.INLA[,"O"]/data.INLA[,"E"],S,J,byrow=F))
      N <- t(chol(Sigma))
      initial.values <- as.vector(c(log(diag(N)), N[lower.tri(N,diag=FALSE)]))
      
      if(prior=="intrinsic"){
        Mmodel <- INLA::inla.rgeneric.define(Mmodel_icar, debug=FALSE, J=J, W=W, initial.values=initial.values)
      }
      if(prior=="Leroux"){
        Mmodel <- INLA::inla.rgeneric.define(Mmodel_lcar, debug=FALSE, J=J, W=W, initial.values=initial.values, alpha.min=0, alpha.max=1)
      }
      if(prior=="proper"){
        Mmodel <- INLA::inla.rgeneric.define(Mmodel_pcar, debug=FALSE, J=J, W=W, initial.values=initial.values, alpha.min=0, alpha.max=1)
      }
      if(prior=="iid"){
        Mmodel <- INLA::inla.rgeneric.define(Mmodel_iid, debug=FALSE, J=J, W=W, initial.values=initial.values)
      }
      
      Model <- inla(formula, family="poisson", data=data.INLA, E=E,
                    control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                    control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                    control.inla=list(strategy=strategy))
      
      Model$Mmodel <- list(model=model, prior=prior)
      
      Mmodel.compute <- Mmodel_compute_cor(Model, n.sample=1000)
      Model$summary.cor <- Mmodel.compute$summary.cor
      Model$marginals.cor <- Mmodel.compute$marginals.cor
      Model$summary.var <- Mmodel.compute$summary.var
      Model$marginals.var <- Mmodel.compute$marginals.var
    }
    
    ## Partition model ##
    if(model=="partition"){
      if(is.null(ID.group))
        stop("the ID.group argument is missing")
      
      cat("STEP 2:",sprintf("Fitting partition (k=%d) model with INLA",k),"\n")
      
      carto.d <- divide_carto(carto, ID.group, k)
      data.d <- lapply(carto.d, function(x) data[data[,ID.area] %in% unlist(sf::st_set_geometry(x[,ID.area],NULL)),])
      
      fun <- function(){
        text <- sprintf("\n%d subdomains(s) have more than 50%% of areas with no observed cases for all the diseases.\nAre you sure that you want to continue fitting the model?\nPress any key to continue or [s] to stop: ",n.zero)
        answer <- readline(cat(red(text," ")))
        if(answer=="s"){
          stop("Stopped by the user.", call.=FALSE)
        }else{
          cat(red("Running...\n"))
        }
      }
      aux <- lapply(data.d, function(xx) aggregate(xx[,O], by=list(xx[,ID.area]), sum)$x)
      prop.zero <- unlist(lapply(aux, function(x) mean(x==0)))
      n.zero <- sum(prop.zero>0.5)
      if(n.zero>0) fun()
      
      invisible(utils::capture.output(aux <- lapply(carto.d, function(x) connect_subgraphs(x, ID.area))))
      Wd <- lapply(aux, function(x) x$W)
      nd <- lapply(Wd, function(x) nrow(x))
      D <- length(nd)
      
      A.constr <- lapply(nd, function(x) kronecker(diag(J), matrix(1,1,x)))
      data.INLA <- mapply(function(x,y){data.frame(O=x[,O], E=x[,E], Area=x[,ID.area], Disease=x[,ID.disease], ID.area=rep(1:y,J), ID.disease=rep(1:J,each=y), idx=seq(1,J*y), ID.group=x[,ID.group])}, x=data.d, y=nd, SIMPLIFY=FALSE)

      intercepts <- lapply(data.INLA, function(x){
        aux <- fastDummies::dummy_cols(x$ID.disease)[,-1]
        aux[aux==0] <- NA
        colnames(aux) <- paste0("I",1:J)
        return(aux)
      })
      
      data.INLA <- mapply(function(x,y){cbind(x,y)}, x=data.INLA, y=intercepts, SIMPLIFY=FALSE)

      Sigma <- lapply(data.INLA, function(x) cov(matrix(x[,"O"]/x[,"E"],ncol=J,byrow=F)))
      N <- lapply(Sigma, function(x) t(chol(x)))
      initial.values <- lapply(N, function(x) as.vector(c(log(diag(x)), x[lower.tri(x,diag=FALSE)])))
      
      if(plan=="sequential"){
        inla.models <- mapply(FitModels, W=Wd, A.constr=A.constr, data.INLA=data.INLA, d=seq(1,D), D=D, initial.values=initial.values, SIMPLIFY=FALSE)
      }
      
      if(plan=="cluster"){
        cl <- future::makeClusterPSOCK(workers, revtunnel=TRUE, outfile="")
        oplan <- future::plan(list(future::tweak(cluster, workers=workers), multisession))
        on.exit(future::plan(oplan))
        
        cpu.time <- system.time({
          inla.models <- future.apply::future_mapply(FitModels, W=Wd, A.constr=A.constr, data.INLA=data.INLA, d=seq(1,D), D=D, initial.values=initial.values, SIMPLIFY=FALSE, future.seed=TRUE)
        })
        
        stopCluster(cl)
      }
      
      if(save.models){
        cat("+ Saving all the inla submodels in '/temp/' folder\n")
        if(!file.exists("temp")) {
          dir.create(file.path(getwd(), "temp"))
        }
        models.dir <- paste("temp/INLAsubmodels_",format(Sys.time(),"%Y%m%d%H%M"),".Rdata",sep="")
        suppressWarnings(save("inla.models", file=models.dir))
      }else{
        models.dir <- NULL
      }

      cat("STEP 3: Merging the results\n")
      Model <- mergeINLA(inla.models=inla.models, ID.disease="Disease", k=k, seed=seed, n.sample=n.sample, compute.intercept=compute.intercept, compute.DIC=compute.DIC, merge.strategy=merge.strategy)

      if(plan=="cluster"){
        Model$cpu.used <- c(Running=as.numeric(cpu.time[3]), Merging=as.numeric(Model$cpu.used["Merging"]), Total=as.numeric(cpu.time[3]+Model$cpu.used["Merging"]))
      }

    }

    return(Model)
    
  }else{
    stop("\nINLA library is not installed!\nPlease use following command to install the stable version of the R-INLA package:\n\ninstall.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }  
}


#' Compute correlation coefficients between diseases
#'
#' @description This function takes a \code{inla} object fitted using the \code{\link{MCAR_INLA}} function and computes the correlation coefficients between diseases. See Details for more information.
#'
#' @param model object of class \code{inla} fitted using the \code{\link{MCAR_INLA}} function.
#' @param n.sample numeric; number of samples to generate from the approximated joint posterior for the hyperparameters (see \code{help(inla.hyperpar.sample)}). Default to 1000.
#' 
#' @return The input \code{inla} object with two additional elements:
#' \item{\code{summary.cor}}{A data.frame containing the mean, standard deviation, quantiles and mode of the correlation coefficients between diseases.}
#' \item{\code{marginals.cor}}{A list containing the posterior marginal densities of the correlation coefficients between diseases.}
#' \item{\code{summary.var}}{A data.frame containing the mean, standard deviation, quantiles and mode of the variances for each disease.}
#' \item{\code{marginals.var}}{A list containing the posterior marginal densities of the variances for each disease.}
#' 
#' @export
Mmodel_compute_cor <- function(model, n.sample=10000){
  
  if(suppressPackageStartupMessages(requireNamespace("INLA", quietly=TRUE))){
    
    if(is.null(model$Mmodel))
      stop("The inla model should be fitted using the MCAR_INLA() function.")
    
    o <- tryCatch.W.E({
      J <- length(unique(model$.args$data$ID.disease))
      hyperpar.sample <- INLA::inla.hyperpar.sample(n.sample, model, improve.marginals=TRUE)
    
      if(model$Mmodel$prior %in% c("intrinsic","iid")){
        hyperpar.sample[,1:J] <- exp(hyperpar.sample[,1:J])
        hyperpar.sample <- split(hyperpar.sample[,seq(J*(J+1)/2)], seq(nrow(hyperpar.sample)))
      }else{
        hyperpar.sample[,seq(J+1,2*J)]<- exp(hyperpar.sample[,seq(J+1,2*J)])
        hyperpar.sample <- split(hyperpar.sample[,seq(1+J,J+J*(J+1)/2)], seq(nrow(hyperpar.sample)))
      }
    
      param.sample <- lapply(hyperpar.sample, function(x){
        N <- diag(x[seq(J)]) 
        N[lower.tri(N, diag=FALSE)] <- x[-seq(J)]
        Sigma <- N %*% t(N)
        Rho <- cov2cor(Sigma)
        Rho.values <- Rho[lower.tri(Rho)]
        
        return(list(sigma=diag(Sigma),rho=Rho.values))
      })
    
      ## Between-disease correlations ##
      cor.sample <- do.call(rbind,lapply(param.sample, function(x) x$rho))
      cor.density <- apply(cor.sample, 2, function(x) density(x, n=75, bw="SJ", from=-1, to=1))
    
      marginals.cor <- lapply(cor.density, function(xx) cbind(x=xx$x, y=xx$y))
      names(marginals.cor) <- paste("rho",apply(combn(J,2), 2, function(x) paste0(x, collapse="")),sep="")
      
      summary.cor <- do.call(rbind,lapply(marginals.cor, compute.summary))

      
      ## Within-disease variances ##
      var.sample <- do.call(rbind,lapply(param.sample, function(x) x$sigma))
      var.density <- apply(var.sample, 2, function(x) density(x, n=75, bw="SJ", from=0))
      
      marginals.var <- lapply(var.density, function(xx) cbind(x=xx$x, y=xx$y))
      names(marginals.var) <- paste("var",1:J,sep="")
      
      summary.var <- do.call(rbind,lapply(marginals.var, compute.summary))
    })
    
    if(any(class(o[[1]])=="error")){
      summary.cor <- data.frame(rep(NA,ncol(combn(J,2))),rep(NA,ncol(combn(J,2))),rep(NA,ncol(combn(J,2))),rep(NA,ncol(combn(J,2))),rep(NA,ncol(combn(J,2))),rep(NA,ncol(combn(J,2))))
      colnames(summary.cor) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode")
      rownames(summary.cor) <- paste("rho",apply(combn(J,2), 2, function(x) paste0(x, collapse="")),sep="")
      
      marginals.cor <- as.list(rep(NA,ncol(combn(J,2))))
      names(marginals.cor) <- rownames(summary.cor)
      
      summary.var <- data.frame(rep(NA,J),rep(NA,J),rep(NA,J),rep(NA,J),rep(NA,J),rep(NA,J))
      colnames(summary.var) <- c("mean","sd","0.025quant","0.5quant","0.975quant","mode")
      rownames(summary.var) <- paste("var",1:J,sep="")
      
      marginals.var <- as.list(rep(NA,J))
      names(marginals.var) <- rownames(summary.var)
    }
    
    Mmodel.compute <- list(summary.cor=summary.cor, marginals.cor=marginals.cor,
                           summary.var=summary.var, marginals.var=marginals.var)
    return(Mmodel.compute)
    
  }else{
    stop("\nINLA library is not installed!\nPlease use following command to install the stable version of the R-INLA package:\n\ninstall.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }
}


############################################
## Auxiliary function to deal with errors ##
############################################
tryCatch.W.E <- function(expr){
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
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

# utils::globalVariables(c("combn"))