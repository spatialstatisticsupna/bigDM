#' Proper multivariate CAR latent effect
#'
#' @description M-model implementation of the proper multivariate CAR latent effect with different spatial autocorrelation parameters using the \code{rgeneric} model of INLA.
#'
#' @details This function considers a proper CAR prior (denoted as pCAR) for the spatial latent effects of the different diseases and introduces correlation between them using the M-model proposal of \insertCite{botella2015unifying;textual}{bigDM}.
#' Putting the spatial latent effects for each disease in a matrix, the between disease dependence is introduced through the M matrix as \eqn{\Theta=\Phi M}, where the columns of \eqn{\Phi} follow a pCAR prior distribution (within-disease correlation).
#' A Wishart prior for the between covariance matrix \eqn{M'M} is considered using the Bartlett decomposition.
#' Uniform prior distributions on the interval [\code{alpha.min}, \code{alpha.max}] are considered for all the spatial autocorrelation parameters.
#' \cr\cr
#' The following arguments are required to be defined before calling the functions:
#' \itemize{
#' \item \code{W}: binary adjacency matrix of the spatial areal units
#' \item \code{J}: number of diseases
#' \item \code{initial.values}: initial values defined for the cells of the M-matrix
#' \item \code{alpha.min}: lower limit defined for the uniform prior distribution of the spatial smoothing parameters
#' \item \code{alpha.max}: upper limit defined for the uniform prior distribution of the spatial smoothing parameters
#' }
#'
#' @note The M-model implementation of this model using R-INLA requires the use of \eqn{J \times (J+3)/2} hyperparameters. So, the results must be carefully checked.
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
## Mmodels - pCAR - FE (BARTLETT DECOMPOSITION)
########################################################################
Mmodel_pcar <- function(cmd=c("graph","Q","mu","initial","log.norm.const","log.prior","quit"), theta=NULL){

  envir <- parent.env(environment())
  if(!exists("cache.done", envir=envir)){
    DiagD <- Matrix::Diagonal(x=colSums(W))

    assign("DiagD", DiagD, envir=envir)
    assign("cache.done", TRUE, envir=envir)
  }

  ########################################################################
  ## theta
  ########################################################################
  interpret.theta <- function(){
    alpha <- alpha.min + (alpha.max-alpha.min)/(1+exp(-theta[as.integer(1:J)]))

    diag.N <- sapply(theta[as.integer(J+1:J)], function(x) { exp(x) })
    no.diag.N <- theta[as.integer(2*J+1:(J*(J-1)/2))]

    N <- diag(diag.N,J)
    N[lower.tri(N, diag=FALSE)] <- no.diag.N

    Covar <- N %*% t(N)

    e <- eigen(Covar)
    M <- t(e$vectors %*% diag(sqrt(e$values)))
    # S <- svd(Covar)
    # M <- t(S$u %*% diag(sqrt(S$d)))

    return(list(alpha=alpha, Covar=Covar, M=M))
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
    MI <- kronecker(M.inv, Matrix::Diagonal(nrow(W)))
    BlockIW <- kronecker(diag(J),DiagD)-kronecker(diag(param$alpha),W)
    Q <- (MI %*% BlockIW) %*% Matrix::t(MI)
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

    ## Uniform prior in (alpha.min, alpha.max) on model scale ##
    val <-  sum(-theta[as.integer(1:J)] - 2*log(1+exp(-theta[as.integer(1:J)])))

    ## n^2_jj ~ chisq(J-j+1) ##
    val <- val + J*log(2) + 2*sum(theta[J+1:J]) + sum(dchisq(exp(2*theta[J+1:J]), df=(J+2)-1:J+1, log=TRUE))

    ## n_ji ~ N(0,1) ##
    val <- val + sum(dnorm(theta[as.integer((2*J)+1:(J*(J-1)/2))], mean=0, sd=1, log=TRUE))

    return(val)
  }

  ########################################################################
  ## initial: return initial values
  ########################################################################
  initial <- function(){
    p <- (0.9-alpha.min)/(alpha.max-alpha.min)

    return(c(rep(log(p/(1-p)),J), as.vector(initial.values)))
  }

  ########################################################################
  ########################################################################
  quit <- function(){ return(invisible()) }

  if(!length(theta)) theta <- initial()
  val <- do.call(match.arg(cmd), args=list())

  return(val)
}

utils::globalVariables(c("alpha.min","alpha.max"))
