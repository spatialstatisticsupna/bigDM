#' Fit a (scalable) spatial Poisson mixed model to areal count data, where several CAR prior distributions can be specified for the spatial random effect.
#'
#' @description Fit a spatial Poisson mixed model to areal count data. The linear predictor is modelled as \deqn{\log{r_{i}}=\alpha+\mathbf{x_i}^{'}\mathbf{\beta} + \xi_i, \quad \mbox{for} \quad i=1,\ldots,n;}
#' where \eqn{\alpha} is a global intercept, \eqn{\mathbf{x_i}^{'}=(x_{i1},\ldots,x_{ip})} is a p-vector of standardized covariates in the i-th area,
#' \eqn{\mathbf{\beta}=(\beta_1,\ldots,\beta_p)} is the p-vector of fixed effects coefficients, and \eqn{\xi_i} is a spatially structured random effect.
#' Several conditional autoregressive (CAR) prior distributions can be specified for the spatial random effect, such as the intrinsic CAR prior \insertCite{besag1991}{bigDM}, the convolution or BYM prior \insertCite{besag1991}{bigDM},
#' the CAR prior proposed by \insertCite{leroux1999estimation;textual}{bigDM}, and the reparameterization of the BYM model given by \insertCite{dean2001detecting;textual}{bigDM} named BYM2 \insertCite{riebler2016intuitive}{bigDM}.
#'
#' If covariates are included in the model, two different approaches can be used to address the potential confounding issues between the fixed effects and the spatial random effects of the model: restricted regression and the use of orthogonality constraints.
#' At the moment, only continuous covariates can be included in the model as potential risk factors, which are automatically standardized before fitting the model. See \insertCite{adin2021alleviating;textual}{bigDM} for further details.
#' \cr\cr
#' Three main modelling approaches can be considered:
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
#' The function allows also to use the new hybrid approximate method that combines the Laplace method with a low-rank Variational Bayes correction to the posterior mean \insertCite{vanNiekerk2023}{bigDM} by including the \code{inla.mode="compact"} argument.
#'
#' @details For a full model specification and further details see the vignettes accompanying this package.
#'
#' @references
#' \insertRef{adin2021alleviating}{bigDM}
#'
#' \insertRef{bengtsson2020unifying}{bigDM}
#'
#' \insertRef{besag1991}{bigDM}
#'
#' \insertRef{dean2001detecting}{bigDM}
#'
#' \insertRef{leroux1999estimation}{bigDM}
#'
#' \insertRef{riebler2016intuitive}{bigDM}
#'
#' \insertRef{rue2009approximate}{bigDM}
#'
#' \insertRef{orozco2020}{bigDM}
#'
#' \insertRef{vanNiekerk2023}{bigDM}
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
#' This object must contain at least the target variables of interest specified in the arguments \code{ID.area}, \code{O} and \code{E}.
#' @param ID.area character; name of the variable that contains the IDs of spatial areal units.
#' @param ID.group character; name of the variable that contains the IDs of the spatial partition (grouping variable).
#' Only required if \code{model="partition"}.
#' @param O character; name of the variable that contains the observed number of disease cases for each areal units.
#' @param E character; name of the variable that contains either the expected number of disease cases or the population at risk for each areal unit.
#' @param X a character vector containing the names of the covariates within the \code{carto} object to be included in the model as fixed effects,
#' or a matrix object playing the role of the fixed effects design matrix. For the latter case, the row names must match with the IDs of the spatial units defined by the \code{ID.area} variable.
#' If \code{X=NULL} (default), only a global intercept is included in the model as fixed effect.
#' @param confounding one of either \code{NULL}, \code{"restricted"} (restricted regression) or \code{"constraints"} (orthogonal constraints), which specifies the estimation method used to alleviate spatial confounding between fixed and random effects.
#' If only an intercept is considered in the model (\code{X=NULL}), the default value \code{confounding=NULL} will be set.
#' At the moment, it only works for the \emph{Global model} (specified through the \code{model="global"} argument).
#' @param W optional argument with the binary adjacency matrix of the spatial areal units.  If \code{NULL} (default), this object is computed from the \code{carto} argument (two areas are considered as neighbours if they share a common border).
#' @param prior one of either \code{"Leroux"} (default), \code{"intrinsic"}, \code{"BYM"} or \code{"BYM2"},
#' which specifies the prior distribution considered for the spatial random effect.
#' @param model one of either \code{"global"} or \code{"partition"} (default), which specifies the \emph{Global model}
#' or one of the scalable model proposal's (\emph{Disjoint model} and \emph{k-order neighbourhood model}, respectively).
#' @param k numeric value with the neighbourhood order used for the partition model. Usually k=2 or 3 is enough to get good results.
#' If k=0 (default) the \emph{Disjoint model} is considered. Only required if \code{model="partition"}.
#' @param strategy one of either \code{"gaussian"}, \code{"simplified.laplace"} (default), \code{"laplace"} or \code{"adaptive"},
#' which specifies the approximation strategy considered in the \code{inla} function.
#' @param PCpriors logical value (default \code{FALSE}); if \code{TRUE} then penalised complexity (PC) priors are used for the precision parameter of the spatial random effect.
#' It only works if arguments \code{prior="intrinsic"} or \code{prior="BYM2"} are specified.
#' @param merge.strategy one of either \code{"mixture"} or \code{"original"} (default), which specifies the merging strategy to compute posterior marginal estimates of the linear predictor. See \code{\link{mergeINLA}} for further details.
#' @param compute.intercept CAUTION! This argument is deprecated from version 0.5.2.
#' @param compute.DIC logical value; if \code{TRUE} (default) then approximate values of the Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) are computed.
#' @param n.sample numeric; number of samples to generate from the posterior marginal distribution of the linear predictor when computing approximate DIC/WAIC values. Default to 1000.
#' @param compute.fitted.values logical value (default \code{FALSE}); if \code{TRUE} transforms the posterior marginal distribution of the linear predictor to the exponential scale (risks or rates).
#' @param save.models logical value (default \code{FALSE}); if \code{TRUE} then a list with all the \code{inla} submodels is saved in '/temp/' folder, which can be used as input argument for the \code{\link{mergeINLA}} function.
#' @param plan one of either \code{"sequential"} or \code{"cluster"}, which specifies the computation strategy used for model fitting using the 'future' package.
#' If \code{plan="sequential"} (default) the models are fitted sequentially and in the current R session (local machine). If \code{plan="cluster"} the models are fitted in parallel on external R sessions (local machine) or distributed in remote computing nodes.
#' @param workers character or vector (default \code{NULL}) containing the identifications of the local or remote workers where the models are going to be processed. Only required if \code{plan="cluster"}.
#' @param inla.mode one of either \code{"classic"} (default) or \code{"compact"}, which specifies the approximation method used by INLA. See \code{help(inla)} for further details.
#' @param num.threads maximum number of threads the inla-program will use. See \code{help(inla)} for further details.
#'
#' @return This function returns an object of class \code{inla}. See the \code{\link{mergeINLA}} function for details.
#'
#' @import crayon Matrix parallel Rdpack
#' @importFrom parallelly makeClusterPSOCK
#' @importFrom future cluster plan
#' @importFrom future.apply future_mapply
#' @importFrom sf st_as_sf st_set_geometry
#' @importFrom stats as.formula
#' @importFrom utils capture.output
#'
#' @examples
#' \dontrun{
#'
#' if(require("INLA", quietly=TRUE)){
#'
#'   ## Load the Spain colorectal cancer mortality data ##
#'   data(Carto_SpainMUN)
#'
#'   ## Global model with a Leroux CAR prior distribution ##
#'   Global <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", O="obs", E="exp",
#'                      prior="Leroux", model="global", strategy="gaussian")
#'
#'   summary(Global)
#'
#'   ## Disjoint model with a Leroux CAR prior distribution  ##
#'   ## using 4 local clusters to fit the models in parallel ##
#'   Disjoint <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
#'                        prior="Leroux", model="partition", k=0, strategy="gaussian",
#'                        plan="cluster", workers=rep("localhost",4))
#'   summary(Disjoint)
#'
#'   ## 1st-order neighbourhood model with a Leroux CAR prior distribution ##
#'   ## using 4 local clusters to fit the models in parallel               ##
#'   order1 <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
#'                      prior="Leroux", model="partition", k=1, strategy="gaussian",
#'                      plan="cluster", workers=rep("localhost",4))
#'   summary(order1)
#'
#'   ## 2nd-order neighbourhood model with a Leroux CAR prior distribution ##
#'   ## using 4 local clusters to fit the models in parallel               ##
#'   order2 <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
#'                      prior="Leroux", model="partition", k=2, strategy="gaussian",
#'                      plan="cluster", workers=rep("localhost",4))
#'   summary(order2)
#' }
#' }
#'
#' @export
CAR_INLA <- function(carto=NULL, ID.area=NULL, ID.group=NULL, O=NULL, E=NULL, X=NULL, confounding=NULL,
                     W=NULL, prior="Leroux", model="partition", k=0, strategy="simplified.laplace",
                     PCpriors=FALSE, merge.strategy="original", compute.intercept=NULL,
                     compute.DIC=TRUE, n.sample=1000, compute.fitted.values=FALSE,
                     save.models=FALSE, plan="sequential", workers=NULL,
                     inla.mode="classic", num.threads=NULL){

  if(suppressPackageStartupMessages(requireNamespace("INLA", quietly=TRUE))){

        ## Set the 'num.threads' argument ##
        if(is.null(num.threads)) num.threads <- INLA::inla.getOption("num.threads")

        ## Check for errors ##
        if(is.null(carto))
                stop("the 'carto' argument is missing")
        if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
                stop("the 'carto' argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
        if(is.null(ID.area))
                stop("the 'ID.area' argument is missing")
        if(is.null(O))
                stop("the 'O' argument is missing")
        if(is.null(E))
                stop("the 'E' argument is missing")
        if(!is.null(confounding)){
                if(!(confounding %in% c("restricted","constraints")))
                        stop("invalid 'confounding' argument")
        }
        if(!(prior %in% c("Leroux","intrinsic","BYM","BYM2")))
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

        if(!missing(compute.intercept)){
                warning("CAUTION! The 'compute.intercept' argument is deprecated from version 0.5.2\n", immediate.=TRUE)
        }

        cat("STEP 1: Pre-processing data\n")

        ## Transform 'SpatialPolygonsDataFrame' object to 'sf' class
        carto <- sf::st_as_sf(carto)

        ## Add the covariates defined in the X argument ##
        if(!is.null(X)){
                if(is.matrix(X)){
                        if(!isTRUE(all.equal(rownames(X),as.character(sf::st_set_geometry(carto, NULL)[,ID.area])))){
                                stop(sprintf("row names of 'X' must match with the IDs of the spatial units defined by the '%s' variable",ID.area))
                        }else{
                                if(is.null(colnames(X))) colnames(X) <- paste("X",seq(ncol(X)),sep="")
                                carto <- cbind(carto,X)
                                X <- colnames(X)
                        }
                }
                if(!all(X %in% names(carto))){
                        stop(sprintf("'%s' variable not found in carto object",X[!X %in% names(carto)]))
                }else{
                        carto[,X] <- scale(sf::st_set_geometry(carto, NULL)[,X])
                }
        }
        data <- sf::st_set_geometry(carto, NULL)

        ## Order the data ##
        if(!ID.area %in% colnames(data))
                stop(sprintf("'%s' variable not found in carto object",ID.area))
        if(!O %in% colnames(data))
                stop(sprintf("'%s' variable not found in carto object",O))
        if(!E %in% colnames(data))
                stop(sprintf("'%s' variable not found in carto object",E))

        if(!all(order(data[,ID.area])==seq(1,nrow(data)))){
                carto <- carto[order(data[,ID.area]),]
                data <- sf::st_set_geometry(carto, NULL)

                order.data <- TRUE
        }else{
                order.data <- FALSE
        }
        rownames(data) <- NULL

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
        form <- "O ~ "
        if(length(X)>0){
                form <- paste(form,paste0(X,collapse="+"),"+")
        }
        if(prior=="Leroux") {
                form <- paste(form,"f(ID.area, model='generic1', Cmatrix=Rs.Leroux, constr=TRUE, hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))", sep="")
        }
        if(prior=="intrinsic" & !PCpriors) {
                form <- paste(form,"f(ID.area, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
        }
        if(prior=="intrinsic" & PCpriors) {
                form <- paste(form,"f(ID.area, model='besag', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
        }
        if(prior=="BYM") {
                form <- paste(form,"f(ID.area, model='bym', graph=Rs, constr=TRUE, hyper=list(theta1=list(prior=sdunif), theta2=list(prior=sdunif)))", sep="")
        }
        if(prior=="BYM2" & !PCpriors) {
                form <- paste(form,"f(ID.area, model='bym2', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif),phi=list(prior=lunif, initial=0)))", sep="")
        }
        if(prior=="BYM2" & PCpriors) {
                form <- paste(form,"f(ID.area, model='bym2', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)),phi=list(prior='pc', param=c(0.5,0.5))))", sep="")
        }

        formula <- stats::as.formula(form)

        ## Auxiliary functions ##
        FitModels <- function(Rs, Rs.Leroux, data.INLA, d, D, ...){

                cat(sprintf("+ Model %d of %d",d,D),"\n")

                Rs <- as(Rs,"Matrix")
                Rs.Leroux <- as(Rs.Leroux,"Matrix")
                # Rs <- inla.as.sparse(Rs)
                # Rs.Leroux <- inla.as.sparse(Rs.Leroux)

                form <- "O ~ "
                if(length(X)>0){
                        form <- paste(form,paste0(X,collapse="+"),"+")
                }
                if(prior=="Leroux") {
                        form <- paste(form,"f(ID.area, model='generic1', Cmatrix=Rs.Leroux, constr=TRUE, hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))", sep="")
                }
                if(prior=="intrinsic" & !PCpriors) {
                        form <- paste(form,"f(ID.area, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
                }
                if(prior=="intrinsic" & PCpriors) {
                        form <- paste(form,"f(ID.area, model='besag', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
                }
                if(prior=="BYM") {
                        form <- paste(form,"f(ID.area, model='bym', graph=Rs, constr=TRUE, hyper=list(theta1=list(prior=sdunif), theta2=list(prior=sdunif)))", sep="")
                }
                if(prior=="BYM2" & !PCpriors) {
                        form <- paste(form,"f(ID.area, model='bym2', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif),phi=list(prior=lunif, initial=0)))", sep="")
                }
                if(prior=="BYM2" & PCpriors) {
                        form <- paste(form,"f(ID.area, model='bym2', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)),phi=list(prior='pc', param=c(0.5,0.5))))", sep="")
                }

                formula <- stats::as.formula(form)

                models <- inla(formula, family="poisson", data=data.INLA, E=E,
                               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                               control.inla=list(strategy=strategy), ...)

                return(models)
        }

        ## Global model ##
        if(model=="global"){
                cat("STEP 2: Fitting global model with INLA (this may take a while...)\n")

                W <- aux$W
                n <- nrow(W)
                Rs <- as(Diagonal(n,colSums(W))-W,"Matrix")
                Rs.Leroux <- as(Diagonal(n)-Rs,"Matrix")
                # Rs <- inla.as.sparse(Diagonal(n,colSums(W))-W)
                # Rs.Leroux <- inla.as.sparse(Diagonal(n)-Rs)

                data.INLA <- data.frame(O=data[,O], E=data[,E], Area=data[,ID.area], ID.area=seq(1,n), data[,X])
                names(data.INLA)[grep("^data...",names(data.INLA))] <- X

                Model <- inla(formula, family="poisson", data=data.INLA, E=E,
                              control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                              control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                              control.inla=list(strategy=strategy),
                              inla.mode=inla.mode, num.threads=num.threads)

                ## Alleviate spatial confounding ##
                if(!is.null(confounding)){
                        time <- Model$cpu.used[4]

                        if(confounding=="restricted"){
                                cat("      + Alleviating spatial confounding: restricted regression model\n")
                                p <- length(X)
                                W <- Diagonal(n, Model$summary.fitted.values$`0.5quant`*data.INLA$E)
                                W.sqrt <- Diagonal(n, sqrt(diag(W)))

                                cat("         * Computing eigen decomposition... ")
                                t.eigen <- system.time({
                                        XX <- cbind(rep(1,n),as.matrix(data.INLA[,X]))
                                        Pc <- Diagonal(n)-W.sqrt%*%XX%*%solve(t(XX)%*%W%*%XX)%*%t(XX)%*%W.sqrt
                                        eigen.Pc <- eigen(Pc)
                                        L <- eigen.Pc$vectors[,eigen.Pc$values>1e-12]
                                        Z.area <- solve(W.sqrt)%*%L%*%t(L)%*%W.sqrt
                                })
                                cat("Elapsed time",as.numeric(t.eigen[3]),"\n")

                                M0 <- solve(t(XX)%*%XX)%*%t(XX)
                                beta.lc = inla.make.lincombs(Predictor=M0, ID.area=-M0%*%Z.area)
                                names(beta.lc) <- paste("X",as.character(0:p),sep="")

                                cat("         * Fitting INLA model... ")
                                Model <- inla(formula, family="poisson", data=data.INLA, E=E,
                                              control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                                              control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                                              lincomb=beta.lc,
                                              control.inla=list(strategy=strategy),
                                              inla.mode=inla.mode, num.threads=num.threads)
                                cat("Done \n")

                                names <- rownames(Model$summary.fixed)
                                Model$summary.fixed <- Model$summary.lincomb.derived[,-1]
                                rownames(Model$summary.fixed) <- names
                                Model$summary.lincomb.derived <- NULL

                                names <- names(Model$marginals.lincomb.derived)
                                Model$marginals.fixed <- Model$marginals.lincomb.derived
                                names(Model$marginals.fixed) <- names
                                Model$marginals.lincomb.derived <- NULL
                        }
                        if(confounding=="constraints"){
                                cat("      + Alleviating spatial confounding: orthogonality constraints model\n")
                                p <- length(X)
                                W <- Diagonal(n, Model$summary.fitted.values$`0.5quant`*data.INLA$E)
                                Bs <- rbind(matrix(1,1,n)%*%W, t(as.matrix(data.INLA[,X]))%*%W)
                                Bs <- as(Bs,"matrix")

                                formula.char <- Reduce(paste,deparse(formula))
                                formula.char <- gsub("constr = TRUE",sprintf("constr = FALSE, rankdef=%d, extraconstr=list(A=Bs, e=rep(0,nrow(Bs)))",p+1),formula.char)
                                formula <- as.formula(formula.char)

                                Model <- inla(formula, family="poisson", data=data.INLA, E=E,
                                              control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                                              control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                                              control.inla=list(strategy=strategy),
                                              inla.mode=inla.mode, num.threads=num.threads)

                                t.eigen <- 0
                        }

                        Model$cpu.used <- c(time+Model$cpu.used[4],t.eigen,t.eigen+time+Model$cpu.used[4])
                        names(Model$cpu.used) <- c("INLA.time","eigen.time","Total")
                }
        }

        ## Partition model ##
        if(model=="partition"){
                if(is.null(ID.group))
                        stop("the ID.group argument is missing")

                cat("STEP 2:",sprintf("Fitting partition (k=%d) model with INLA",k),"\n")

                carto.d <- divide_carto(carto, ID.group, k)
                data.d <- lapply(carto.d, function(x) sf::st_set_geometry(x, NULL))

                fun <- function(){
                        text <- sprintf("\n%d subdomains(s) have more than 70%% of areas with no observed cases.\nAre you sure that you want to continue fitting the model?\nPress any key to continue or [s] to stop: ",n.zero)
                        answer <- readline(cat(red(text," ")))
                        if(answer=="s"){
                                stop("Stopped by the user.", call.=FALSE)
                        }else{
                                cat(red("Running...\n"))
                        }
                }
                prop.zero <- unlist(lapply(data.d, function(x) mean(x[,O]==0)))
                n.zero <- sum(prop.zero>0.7)
                if(n.zero>0) fun()

                invisible(utils::capture.output(aux <- lapply(carto.d, function(x) connect_subgraphs(x, ID.area))))
                Wd <- lapply(aux, function(x) x$W)
                nd <- lapply(Wd, function(x) nrow(x))
                Rs <- mapply(function(x,y){Diagonal(x,colSums(y))-y}, x=nd, y=Wd)
                Rs.Leroux <- mapply(function(x,y){Diagonal(x)-y}, x=nd, y=Rs)

                data.INLA <- mapply(function(x,y){
                        aux <- data.frame(O=x[,O], E=x[,E], Area=x[,ID.area], ID.area=seq(1,y), ID.group=x[,ID.group], x[,X])
                        names(aux)[grep("^x...",names(aux))] <- X
                        aux
                }, x=data.d, y=nd, SIMPLIFY=FALSE)
                D <- length(data.INLA)

                if(plan=="sequential"){
                        inla.models <- mapply(FitModels, Rs=Rs, Rs.Leroux=Rs.Leroux, data.INLA=data.INLA, d=seq(1,D), D=D, num.threads=num.threads, inla.mode=inla.mode, SIMPLIFY=FALSE)
                }

                if(plan=="cluster"){
                        cl <- makeClusterPSOCK(workers, revtunnel=TRUE, outfile="")
                        oplan <- plan(cluster, workers=cl)
                        # oplan <- plan(list(tweak(cluster, workers=workers), multisession))
                        on.exit(plan(oplan))

                        cpu.time <- system.time({
                                inla.models <- future_mapply(FitModels, Rs=Rs, Rs.Leroux=Rs.Leroux, data.INLA=data.INLA, d=seq(1,D), D=D, num.threads=num.threads, inla.mode=inla.mode, future.seed=TRUE, SIMPLIFY=FALSE)
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
                Model <- mergeINLA(inla.models=inla.models, k=k, n.sample=n.sample, compute.DIC=compute.DIC, merge.strategy=merge.strategy, compute.fitted.values=compute.fitted.values)

                if(plan=="cluster"){
                        Model$cpu.used <- c(Running=as.numeric(cpu.time[3]), Merging=as.numeric(Model$cpu.used["Merging"]), Total=as.numeric(cpu.time[3]+Model$cpu.used["Merging"]))
                }
        }

        if(order.data) warning("CAUTION: The input dataset has been sorted by 'ID.area' variable to ensure correct model fitting. Please, check the model$.args$data object for details.\n", call. = FALSE)

        return(Model)

  }else{
        stop("\nINLA library is not installed!\nPlease use the following command to install the stable version of the R-INLA package:\n\ninstall.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }
}

# utils::globalVariables(c("inla.as.sparse"))
utils::globalVariables(c("inla.setOption"))
