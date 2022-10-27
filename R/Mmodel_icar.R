#' Intrinsic multivariate CAR latent effect
#'
#' @description M-model implementation of the intrinsic multivariate CAR latent effect using the \code{rgeneric} model of INLA.
#'
#' @details This function considers an intrinsic CAR prior for the spatial latent effects of the different diseases and introduces correlation between them using the M-model proposal of \insertCite{botella2015unifying;textual}{bigDM}.
#' Putting the spatial latent effects for each disease in a matrix, the between disease dependence is introduced through the M matrix as \eqn{\Theta=\Phi M}, where the columns of \eqn{\Phi} follow an intrinsic CAR prior distribution (within-disease correlation).
#' A Wishart prior for the between covariance matrix \eqn{M'M} is considered using the Bartlett decomposition.
#' \cr\cr
#' The following arguments are required to be defined before calling the functions:
#' \itemize{
#' \item \code{W}: binary adjacency matrix of the spatial areal units
#' \item \code{J}: number of diseases
#' \item \code{initial.values}: initial values defined for the cells of the M-matrix
#' }
#'
#' @references
#' \insertRef{botella2015unifying}{bigDM}
#'
#' @param cmd Internal functions used by the \code{rgeneric} model to define the latent effect.
#' @param theta Vector of hyperparameters.
#'
#' @return This is used internally by the \code{INLA::inla.rgeneric.define()} function.
#' 
#' @import Matrix
#' 
#'
#' @export
########################################################################
## Mmodels - iCAR - FE (BARTLETT DECOMPOSITION)
########################################################################
Mmodel_icar <- function(cmd=c("graph","Q","mu","initial","log.norm.const","log.prior","quit"), theta=NULL){
  
  envir <- parent.env(environment())
  if(!exists("cache.done", envir=envir)){
    BlockIW <- kronecker(Diagonal(J),Diagonal(x=colSums(W))-W)
    
    assign("BlockIW", BlockIW, envir=envir)
    assign("cache.done", TRUE, envir=envir)
  }
  
  ########################################################################
  ## theta
  ########################################################################
  interpret.theta <- function(){
    diag.N <- sapply(theta[as.integer(1:J)], function(x) { exp(x) })
    no.diag.N <- theta[as.integer(J+1:(J*(J-1)/2))]
    
    N <- diag(diag.N,J) 
    N[lower.tri(N, diag=FALSE)] <- no.diag.N
    
    Covar <- N %*% t(N)
    
    e <- eigen(Covar)
    M <- t(e$vectors %*% diag(sqrt(e$values)))
    # S <- svd(Covar)
    # M <- t(S$u %*% diag(sqrt(S$d)))

    return(list(Covar=Covar, M=M))
  }
  
  ########################################################################
  ## Graph of precision function; i.e., a 0/1 representation of precision matrix
  ########################################################################
  graph <- function(){ return(Q()) }
  
  ########################################################################
  ## Precision matrix
  ########################################################################
  Q <- function(){
    param <- interpret.theta()
    M.inv <- solve(param$M)
    MI <- kronecker(M.inv, Diagonal(nrow(W)))
    Q <- (MI %*% BlockIW) %*% t(MI)
    return(Q)
  }
  
  ########################################################################
  ## Mean of model
  ########################################################################
  mu <- function(){ return(numeric(0)) }
  
  ########################################################################
  ## log.norm.const
  ########################################################################
  log.norm.const <- function(){
    val <- numeric(0)
    return(val)
  }
  
  ########################################################################
  ## log.prior: return the log-prior for the hyperparameters
  ########################################################################
  log.prior <- function(){
    param <- interpret.theta()
    
    ## n^2_jj ~ chisq(J-j+1) ##
    val <- J*log(2) + 2 *sum(theta[1:J]) + sum(dchisq(exp(2*theta[1:J]), df=(J+2)-1:J+1, log=TRUE))
    
    ## n_ji ~ N(0,1) ##
    val <- val + sum(dnorm(theta[as.integer(J+1:(J*(J-1)/2))], mean=0, sd=1, log=TRUE))
    
    return(val)
  }

  ########################################################################
  ## initial: return initial values
  ########################################################################
  initial <- function(){ return(as.vector(initial.values)) }
  
  ########################################################################
  ########################################################################
  quit <- function(){ return(invisible()) }
  
  if(!length(theta)) theta <- initial()
  val <- do.call(match.arg(cmd), args=list())
  
  return(val)
}

utils::globalVariables(c("J","W","interpret.theta","Q","log.prior","initial",
                         "dchisq","dnorm","initial.values"))