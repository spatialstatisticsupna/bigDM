#' Fit a (scalable) spatio-temporal Poisson mixed model to areal count data.
#'
#' @description Fit a spatio-temporal Poisson mixed model to areal count data, where several CAR prior distributions for the spatial random effects, first and second order random walk priors for the temporal random effects, and different types of spatio-temporal interactions described in \insertCite{knorrheld2000;textual}{bigDM} can be specified.
#' The linear predictor is modelled as \deqn{\log{r_{it}}=\alpha+\xi_i+\gamma_t+\delta_{it}, \quad \mbox{for} \quad i=1,\ldots,n; \quad t=1,\ldots,T}
#' where \eqn{\alpha} is a global intercept, \eqn{\xi_i} is a spatially structured random effect, \eqn{\gamma_t} is a temporally structured random effect, and \eqn{\delta_{it}} is the space-time interaction effect. If the interaction term is dropped, an additive model is obtained.
#' To ensure model identifiability, sum-to-zero constraints are imposed over the random effects of the model. Details on the derivation of these constraints can be found in \insertCite{goicoa2018spatio;textual}{bigDM}.
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
#' For the scalable model proposals \insertCite{orozco2022}{bigDM}, approximate values of the Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) can also be computed.
#'
#' @details For a full model specification and further details see the vignettes accompanying this package.
#'
#' @references
#' \insertRef{goicoa2018spatio}{bigDM}
#'
#' \insertRef{knorrheld2000}{bigDM}
#'
#' \insertRef{orozco2020}{bigDM}
#'
#' \insertRef{orozco2022}{bigDM}
#'
#' @param carto object of class \code{SpatialPolygonsDataFrame} or \code{sf}. This object must contain at least the variable with the identifiers of the spatial areal units specified in the argument \code{ID.area}.
#' @param data object of class \code{data.frame} that must contain the target variables of interest specified in the arguments \code{ID.area}, \code{ID.year}, \code{O} and \code{E}.
#' @param ID.area character; name of the variable that contains the IDs of spatial areal units. The values of this variable must match those given in the \code{carto} and \code{data} variable.
#' @param ID.year character; name of the variable that contains the IDs of time points.
#' @param ID.group character; name of the variable that contains the IDs of the spatial partition (grouping variable). Only required if \code{model="partition"}.
#' @param O character; name of the variable that contains the observed number of disease cases for each areal and time point.
#' @param E character; name of the variable that contains either the expected number of disease cases or the population at risk for each areal unit and time point.
#' @param W optional argument with the binary adjacency matrix of the spatial areal units. If \code{NULL} (default), this object is computed from the \code{carto} argument (two areas are considered as neighbours if they share a common border).
#' @param spatial one of either \code{"Leroux"} (default), \code{"intrinsic"}, \code{"BYM"} or \code{"BYM2"}, which specifies the prior distribution considered for the spatial random effect.
#' @param temporal one of either \code{"rw1"} (default) or \code{"rw2"}, which specifies the prior distribution considered for the temporal random effect.
#' @param interaction one of either \code{"none"}, \code{"TypeI"}, \code{"TypeII"}, \code{"TypeIII"} or \code{"TypeIV"} (default), which specifies the prior distribution for the space-time interaction random effect.
#' @param model one of either \code{"global"} or \code{"partition"} (default), which specifies the \emph{Global model} or one of the scalable model proposal's (\emph{Disjoint model} and \emph{k-order neighbourhood model}, respectively).
#' @param k numeric value with the neighbourhood order used for the partition model. Usually k=2 or 3 is enough to get good results. If k=0 (default) the \emph{Disjoint model} is considered. Only required if \code{model="partition"}.
#' @param strategy one of either \code{"gaussian"}, \code{"simplified.laplace"} (default), \code{"laplace"} or \code{"adaptive"}, which specifies the approximation strategy considered in the \code{inla} function.
#' @param PCpriors logical value (default \code{FALSE}); if \code{TRUE} then penalised complexity (PC) priors are used for the precision parameter of the spatial random effect.
#' It only works if arguments \code{spatial="intrinsic"} or \code{spatial="BYM2"} are specified.
#' @param seed numeric (default \code{NULL}); control the RNG of the \code{inla.qsample} function. See \code{help(inla.qsample)} for further information.
#' @param n.sample numeric; number of samples to generate from the posterior marginal distribution of the risks. Default to 1000.
#' @param compute.intercept logical value (default \code{FALSE}); if \code{TRUE} then the overall log-risk \eqn{\alpha} is computed.
#' It only works if \code{k=0} argument (\emph{Disjoint model}) is specified. CAUTION: This method might be very time consuming.
#' @param compute.DIC logical value; if \code{TRUE} (default) then approximate values of the Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) are computed.
#' @param merge.strategy one of either \code{"mixture"} or \code{"original"} (default), which specifies the merging strategy to compute posterior marginal estimates of relative risks. See \code{\link{mergeINLA}} for further details.
#' @param save.models logical value (default \code{FALSE}); if \code{TRUE} then a list with all the \code{inla} submodels is saved in '/temp/' folder, which can be used as input argument for the \code{\link{mergeINLA}} function.
#' @param plan one of either \code{"sequential"} or \code{"cluster"}, which specifies the computation strategy used for model fitting using the 'future' package.
#' If \code{plan="sequential"} (default) the models are fitted sequentially and in the current R session (local machine). If \code{plan="cluster"} the models are fitted in parallel on external R sessions (local machine) or distributed in remote computing nodes.
#' @param workers character or vector (default \code{NULL}) containing the identifications of the local or remote workers where the models are going to be processed. Only required if \code{plan="cluster"}.
#' @param inla.mode one of either \code{"classic"} (default) or \code{"compact"}, which specifies the approximation method used by INLA. See \code{help(inla)} for further details.
#'
#' @return This function returns an object of class \code{inla}. See the \code{\link{mergeINLA}} function for details.
#'
#' @import crayon future Matrix parallel
#' @importFrom future.apply future_mapply
#' @importFrom MASS ginv
#' @importFrom sf st_as_sf st_set_geometry
#' @importFrom stats as.formula
#' @importFrom utils capture.output
#' @importFrom methods as
#'
#' @examples
#' \dontrun{
#' if(require("INLA", quietly=TRUE)){
#'
#'   ## Load the sf object that contains the spatial polygons of the municipalities of Spain ##
#'   data(Carto_SpainMUN)
#'   str(Carto_SpainMUN)
#'
#'   ## Create province IDs ##
#'   Carto_SpainMUN$ID.prov <- substr(Carto_SpainMUN$ID,1,2)
#'
#'   ## Load simulated data of lung cancer mortality data during the period 1991-2015 ##
#'   data("Data_LungCancer")
#'   str(Data_LungCancer)
#'
#'   ## Disjoint model with a BYM2 spatial random effect, RW1 temporal random effect and      ##
#'   ## Type I interaction random effect using 4 local clusters to fit the models in parallel ##
#'   Disjoint <- STCAR_INLA(carto=Carto_SpainMUN, data=Data_LungCancer,
#'                          ID.area="ID", ID.year="year", O="obs", E="exp", ID.group="ID.prov",
#'                          spatial="BYM2", temporal="rw1", interaction="TypeI",
#'                          model="partition", k=0, strategy="gaussian",
#'                          plan="cluster", workers=rep("localhost",4))
#'   summary(Disjoint)
#'
#'   ## 1st-order nb. model with a BYM2 spatial random effect, RW1 temporal random effect and ##
#'   ## Type I interaction random effect using 4 local clusters to fit the models in parallel ##
#'   order1 <- STCAR_INLA(carto=Carto_SpainMUN, data=Data_LungCancer,
#'                        ID.area="ID", ID.year="year", O="obs", E="exp", ID.group="ID.prov",
#'                        spatial="BYM2", temporal="rw1", interaction="TypeI",
#'                        model="partition", k=1, strategy="gaussian",
#'                        plan="cluster", workers=rep("localhost",4))
#'   summary(order1)
#' }
#' }
#'
#' @export
STCAR_INLA <- function(carto=NULL, data=NULL, ID.area=NULL, ID.year=NULL, ID.group=NULL, O=NULL, E=NULL,
                       W=NULL, spatial="Leroux", temporal="rw1", interaction="TypeIV",
                       model="partition", k=0, strategy="simplified.laplace",
                       PCpriors=FALSE, seed=NULL, n.sample=1000, compute.intercept=FALSE, compute.DIC=TRUE,
                       save.models=FALSE, plan="sequential", workers=NULL, merge.strategy="original",
                       inla.mode="classic"){

  if(suppressPackageStartupMessages(requireNamespace("INLA", quietly=TRUE))){

        ## Set the "inla.mode" argument ##
        inla.setOption(inla.mode=inla.mode)

        ## Check for errors ##
        if(is.null(carto))
                stop("the 'carto' argument is missing")
        if(!any(class(carto) %in% c("SpatialPolygonsDataFrame","sf")))
                stop("the 'carto' argument must be of class 'SpatialPolygonsDataFrame' or 'sf'")
        if(is.null(ID.area))
                stop("the 'ID.area' argument is missing")
        if(is.null(ID.year))
                stop("the 'ID.year' argument is missing")
        if(is.null(O))
                stop("the 'O' argument is missing")
        if(is.null(E))
                stop("the 'E' argument is missing")
        if(!(spatial %in% c("Leroux","intrinsic","BYM","BYM2")))
                stop("invalid 'spatial' argument")
        if(!(temporal %in% c("rw1","rw2")))
                stop("invalid 'temporal' argument")
        if(!(interaction %in% c("none","TypeI","TypeII","TypeIII","TypeIV")))
                stop("invalid 'interaction' argument")
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
        if(!ID.year %in% colnames(data))
                stop(sprintf("'%s' variable not found in data object",ID.year))
        if(!O %in% colnames(data))
                stop(sprintf("'%s' variable not found in carto object",O))
        if(!E %in% colnames(data))
                stop(sprintf("'%s' variable not found in carto object",E))

        carto <- carto[order(sf::st_set_geometry(carto, NULL)[,ID.area]),]
        data <- merge(data,carto[,c(ID.area,ID.group)])
        data$geometry <- NULL
        data[,ID.year] <- paste(sprintf("%02d", as.numeric(as.character(data[,ID.year]))))
        data <- data[order(data[,ID.year],data[,ID.area]),]

        ## Merge disjoint connected subgraphs ##
        if(is.null(W)){
                invisible(utils::capture.output(aux <- connect_subgraphs(carto, ID.area)))
                carto.nb <- aux$nb
        }else{
                carto.nb <- spdep::mat2listw(W, style="B")$neighbours
                invisible(utils::capture.output(aux <- connect_subgraphs(carto, ID.area, nb=carto.nb)))
                carto.nb <- aux$nb
        }
        S <- length(carto.nb)

        ## Temporal structure matrix ##
        T <- length(unique(data[,ID.year]))
        if(temporal=="rw1") dif <- 1
        if(temporal=="rw2") dif <- 2
        D <- diff(diag(T), differences=dif)
        Rt <- as(t(D)%*%D, "TsparseMatrix")
        # Rt <- inla.as.sparse(t(D)%*%D)

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

        ## Identifiability constraints ##
        constraints <- function(Rs, Rt) {
                S <- nrow(Rs)
                T <- nrow(Rt)

                if(interaction %in% c("none","TypeI")){
                        R <- NULL
                        r.def <- NULL
                        A.constr <- NULL
                }
                if(interaction=="TypeII"){
                        R <- kronecker(Rt,Diagonal(S))
                        if(temporal=="rw1") r.def <- S
                        if(temporal=="rw2") r.def <- 2*S
                        A.constr <- kronecker(matrix(1,1,T),Diagonal(S))
                        A.constr <- as(A.constr[-1,],"matrix")

                        if(PCpriors){
                                sigma.ref <- exp(mean(log(diag(MASS::ginv(as.matrix(Rt))))))
                                R <- R*sigma.ref
                        }
                }
                if(interaction=="TypeIII"){
                        R <- kronecker(Diagonal(T),Rs)
                        r.def <- T
                        A.constr <- kronecker(Diagonal(T),matrix(1,1,S))
                        A.constr <- as(A.constr[-1,],"matrix")

                        if(PCpriors){
                                sigma.ref <- exp(mean(log(diag(MASS::ginv(as.matrix(Rs))))))
                                R <- R*sigma.ref
                        }
                }
                if(interaction=="TypeIV"){
                        R <- kronecker(Rt,Rs)
                        if(temporal=="rw1") r.def <- S+T-1
                        if(temporal=="rw2") r.def <- 2*S+T-2
                        A1 <- kronecker(matrix(1,1,T),Diagonal(S))
                        A2 <- kronecker(Diagonal(T),matrix(1,1,S))
                        A.constr <- as(rbind(A1[-1,],A2[-1,]),"matrix")
                        if(PCpriors){
                                sigma.ref1 <- exp(mean(log(diag(MASS::ginv(as.matrix(Rs))))))
                                sigma.ref2 <- exp(mean(log(diag(MASS::ginv(as.matrix(Rt))))))
                                R <- R*sigma.ref1*sigma.ref2
                        }
                }

                return(list(R=R,r.def=r.def,A.constr=A.constr))
        }

        ## Formula for INLA model ##
        form <- "O ~ "

        if(spatial=="Leroux") {
                form <- paste(form,"f(ID.area, model='generic1', Cmatrix=Rs.Leroux, constr=TRUE, hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))", sep="")
        }
        if(spatial=="intrinsic" & !PCpriors) {
                form <- paste(form,"f(ID.area, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
        }
        if(spatial=="intrinsic" & PCpriors) {
                form <- paste(form,"f(ID.area, model='besag', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
        }
        if(spatial=="BYM") {
                form <- paste(form,"f(ID.area, model='bym', graph=Rs, constr=TRUE, hyper=list(theta1=list(prior=sdunif), theta2=list(prior=sdunif)))", sep="")
        }
        if(spatial=="BYM2" & !PCpriors) {
                form <- paste(form,"f(ID.area, model='bym2', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif),phi=list(prior=lunif, initial=0)))", sep="")
        }
        if(spatial=="BYM2" & PCpriors) {
                form <- paste(form,"f(ID.area, model='bym2', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)),phi=list(prior='pc', param=c(0.5,0.5))))", sep="")
        }

        if(temporal=="rw1" & !PCpriors) {
                form <- paste(form,"+ f(ID.year, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
        }
        if(temporal=="rw1" & PCpriors) {
                form <- paste(form,"+ f(ID.year, model='rw1', constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
        }
        if(temporal=="rw2" & interaction!="TypeII" & !PCpriors) {
                form <- paste(form,"+ f(ID.year, model='rw2', constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
        }
        if(temporal=="rw2" & interaction!="TypeII" & PCpriors) {
                form <- paste(form,"+ f(ID.year, model='rw2', constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
        }
        if(temporal=="rw2" & interaction=="TypeII" & !PCpriors) {
                form <- paste(form,"+ f(ID.year, model='rw2', constr=TRUE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=matrix(1:T,1,T),e=0))", sep="")
        }
        if(temporal=="rw2" & interaction=="TypeII" & PCpriors) {
                form <- paste(form,"+ f(ID.year, model='rw2', constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))), extraconstr=list(A=matrix(1:T,1,T),e=0))", sep="")
        }

        if(interaction=="TypeI" & temporal=="rw1" & !PCpriors) {
                form <- paste(form,"+ f(ID.area.year, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
        }
        if(interaction=="TypeI" & temporal=="rw1" & PCpriors) {
                form <- paste(form,"+ f(ID.area.year, model='iid', constr=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
        }
        if(interaction=="TypeI" & temporal=="rw2" & !PCpriors) {
                form <- paste(form,"+ f(ID.area.year, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=matrix(rep(1:S,each=T),1,S*T),e=0))", sep="")
        }
        if(interaction=="TypeI" & temporal=="rw2" & PCpriors) {
                form <- paste(form,"+ f(ID.area.year, model='iid', constr=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))), extraconstr=list(A=matrix(rep(1:S,each=T),1,S*T),e=0))", sep="")
        }
        if((interaction %in% c("TypeII","TypeIII","TypeIV")) & !PCpriors) {
                form <- paste(form,"+ f(ID.area.year, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))", sep="")
        }
        if((interaction %in% c("TypeII","TypeIII","TypeIV")) & PCpriors) {
                form <- paste(form,"+ f(ID.area.year, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))), extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))", sep="")
        }

        formula <- stats::as.formula(form)

        ## Auxiliary functions ##
        FitModels <- function(Rs, Rs.Leroux, R, r.def, A.constr, data.INLA, d, D){

                cat(sprintf("+ Model %d of %d",d,D),"\n")

                Rs <- as(Rs,"TsparseMatrix")
                Rs.Leroux <- as(Rs.Leroux,"TsparseMatrix")
                # Rs <- inla.as.sparse(Rs)
                # Rs.Leroux <- inla.as.sparse(Rs.Leroux)
                S <- nrow(Rs)

                form <- "O ~ "

                if(spatial=="Leroux") {
                        form <- paste(form,"f(ID.area, model='generic1', Cmatrix=Rs.Leroux, constr=TRUE, hyper=list(prec=list(prior=sdunif),beta=list(prior=lunif, initial=0)))", sep="")
                }
                if(spatial=="intrinsic" & !PCpriors) {
                        form <- paste(form,"f(ID.area, model='besag', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
                }
                if(spatial=="intrinsic" & PCpriors) {
                        form <- paste(form,"f(ID.area, model='besag', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
                }
                if(spatial=="BYM") {
                        form <- paste(form,"f(ID.area, model='bym', graph=Rs, constr=TRUE, hyper=list(theta1=list(prior=sdunif), theta2=list(prior=sdunif)))", sep="")
                }
                if(spatial=="BYM2" & !PCpriors) {
                        form <- paste(form,"f(ID.area, model='bym2', graph=Rs, constr=TRUE, hyper=list(prec=list(prior=sdunif),phi=list(prior=lunif, initial=0)))", sep="")
                }
                if(spatial=="BYM2" & PCpriors) {
                        form <- paste(form,"f(ID.area, model='bym2', graph=Rs, constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01)),phi=list(prior='pc', param=c(0.5,0.5))))", sep="")
                }

                if(temporal=="rw1" & !PCpriors) {
                        form <- paste(form,"+ f(ID.year, model='rw1', constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
                }
                if(temporal=="rw1" & PCpriors) {
                        form <- paste(form,"+ f(ID.year, model='rw1', constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
                }
                if(temporal=="rw2" & interaction!="TypeII" & !PCpriors) {
                        form <- paste(form,"+ f(ID.year, model='rw2', constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
                }
                if(temporal=="rw2" & interaction!="TypeII" & PCpriors) {
                        form <- paste(form,"+ f(ID.year, model='rw2', constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
                }
                if(temporal=="rw2" & interaction=="TypeII" & !PCpriors) {
                        form <- paste(form,"+ f(ID.year, model='rw2', constr=TRUE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=matrix(1:T,1,T),e=0))", sep="")
                }
                if(temporal=="rw2" & interaction=="TypeII" & PCpriors) {
                        form <- paste(form,"+ f(ID.year, model='rw2', constr=TRUE, scale.model=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))), extraconstr=list(A=matrix(1:T,1,T),e=0))", sep="")
                }

                if(interaction=="TypeI" & temporal=="rw1" & !PCpriors) {
                        form <- paste(form,"+ f(ID.area.year, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)))", sep="")
                }
                if(interaction=="TypeI" & temporal=="rw1" & PCpriors) {
                        form <- paste(form,"+ f(ID.area.year, model='iid', constr=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))))", sep="")
                }
                if(interaction=="TypeI" & temporal=="rw2" & !PCpriors) {
                        form <- paste(form,"+ f(ID.area.year, model='iid', constr=TRUE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=matrix(rep(1:S,each=T),1,S*T),e=0))", sep="")
                }
                if(interaction=="TypeI" & temporal=="rw2" & PCpriors) {
                        form <- paste(form,"+ f(ID.area.year, model='iid', constr=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))), extraconstr=list(A=matrix(rep(1:S,each=T),1,S*T),e=0))", sep="")
                }
                if((interaction %in% c("TypeII","TypeIII","TypeIV")) & !PCpriors) {
                        form <- paste(form,"+ f(ID.area.year, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior=sdunif)), extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))", sep="")
                }
                if((interaction %in% c("TypeII","TypeIII","TypeIV")) & PCpriors) {
                        form <- paste(form,"+ f(ID.area.year, model='generic0', Cmatrix=R, rankdef=r.def, constr=TRUE, hyper=list(prec=list(prior='pc.prec', param=c(1,0.01))), extraconstr=list(A=A.constr, e=rep(0,nrow(A.constr))))", sep="")
                }

                formula <- stats::as.formula(form)

                models <- inla(formula, family="poisson", data=data.INLA, E=E,
                               control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                               control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                               control.inla=list(strategy=strategy))
                return(models)
        }

        ## Global model ##
        if(model=="global"){
                cat("STEP 2: Fitting global model with INLA (this may take a while...)\n")

                W <- aux$W
                Rs <- as(Diagonal(S,colSums(W))-W, "TsparseMatrix")
                Rs.Leroux <- as(Diagonal(S)-Rs, "TsparseMatrix")
                # Rs <- inla.as.sparse(Diagonal(S,colSums(W))-W)
                # Rs.Leroux <- inla.as.sparse(Diagonal(S)-Rs)

                constr <- constraints(Rs,Rt)
                R <- constr$R
                r.def <- constr$r.def
                A.constr <- constr$A.constr

                data.INLA <- data.frame(O=data[,O], E=data[,E], Area=data[,ID.area], Year=data[,ID.year],
                                        ID.area=rep(1:S,T), ID.year=rep(1:T,each=S), ID.area.year=seq(1,T*S))

                Model <- inla(formula, family="poisson", data=data.INLA, E=E,
                              control.predictor=list(compute=TRUE, link=1, cdf=c(log(1))),
                              control.compute=list(dic=TRUE, cpo=TRUE, waic=TRUE, config=TRUE, return.marginals.predictor=TRUE),
                              control.inla=list(strategy=strategy))
        }

        ## Partition model ##
        if(model=="partition"){
                if(is.null(ID.group))
                        stop("the ID.group argument is missing")

                cat("STEP 2:",sprintf("Fitting partition (k=%d) model with INLA",k),"\n")

                carto.d <- divide_carto(carto, ID.group, k)
                data.d <- lapply(carto.d, function(x) data[data[,ID.area] %in% unlist(sf::st_set_geometry(x[,ID.area],NULL)),])

                fun <- function(){
                        text <- sprintf("\n%d subdomains(s) have more than 50%% of areas with no observed cases during the whole study period.\nAre you sure that you want to continue fitting the model?\nPress any key to continue or [s] to stop: ",n.zero)
                        answer <- readline(cat(red(text," ")))
                        if(answer=="s"){
                                stop("Stopped by the user.", call.=FALSE)
                        }else{
                                cat(red("Running...\n"))
                        }
                }
                aux <- lapply(data.d, function(xx) aggregate(xx[,O], by=list(xx[,ID.area]), sum, na.rm=T)$x)
                prop.zero <- unlist(lapply(aux, function(x) mean(x==0)))
                n.zero <- sum(prop.zero>0.5)
                if(n.zero>0) fun()

                invisible(utils::capture.output(aux <- lapply(carto.d, function(x) connect_subgraphs(x, ID.area))))
                Wd <- lapply(aux, function(x) x$W)
                nd <- lapply(Wd, function(x) nrow(x))
                Rs <- mapply(function(x,y){Diagonal(x,colSums(y))-y}, x=nd, y=Wd)
                Rs.Leroux <- mapply(function(x,y){Diagonal(x)-y}, x=nd, y=Rs)

                constr <- lapply(Rs, function(x) constraints(x,Rt))
                R <- lapply(constr, function(x) x$R)
                r.def <- lapply(constr, function(x) x$r.def)
                A.constr <- lapply(constr, function(x) x$A.constr)

                data.INLA <- mapply(function(x,y){data.frame(O=x[,O], E=x[,E], Area=x[,ID.area], Year=x[,ID.year], ID.group=x[,ID.group], ID.area=rep(1:y,T), ID.year=rep(1:T,each=y), ID.area.year=seq(1,T*y))}, x=data.d, y=nd, SIMPLIFY=FALSE)
                D <- length(data.INLA)

                if(plan=="sequential"){
                        inla.models <- mapply(FitModels, Rs=Rs, Rs.Leroux=Rs.Leroux, R=R, r.def=r.def, A.constr=A.constr, data.INLA=data.INLA, d=seq(1,D), D=D, SIMPLIFY=FALSE)
                }

                if(plan=="cluster"){
                        cl <- future::makeClusterPSOCK(workers, revtunnel=TRUE, outfile="")
                        oplan <- future::plan(list(future::tweak(cluster, workers=workers), multisession))
                        on.exit(future::plan(oplan))

                        cpu.time <- system.time({
                                inla.models <- future.apply::future_mapply(FitModels, Rs=Rs, Rs.Leroux=Rs.Leroux, R=R, r.def=r.def, A.constr=A.constr, data.INLA=data.INLA, d=seq(1,D), D=D, SIMPLIFY=FALSE, future.seed=TRUE)
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
                Model <- mergeINLA(inla.models=inla.models, ID.year="Year", k=k, seed=seed, n.sample=n.sample, compute.intercept=compute.intercept, compute.DIC=compute.DIC, merge.strategy=merge.strategy)

                if(plan=="cluster"){
                        Model$cpu.used <- c(Running=as.numeric(cpu.time[3]), Merging=as.numeric(Model$cpu.used["Merging"]), Total=as.numeric(cpu.time[3]+Model$cpu.used["Merging"]))
                }
        }

        return(Model)

  }else{
        stop("INLA library is not installed! Please use following command to install the stable version of the R-INLA package:\n install.packages('INLA', repos=c(getOption('repos'), INLA='https://inla.r-inla-download.org/R/stable'), dep=TRUE)")
  }
}

# utils::globalVariables(c("inla.as.sparse"))
utils::globalVariables(c("inla.setOption"))
