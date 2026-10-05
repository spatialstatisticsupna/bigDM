#' BYM2 multivariate CAR latent effect
#'
#' @description M-model implementation of the BYM2 multivariate CAR latent effect using the \code{rgeneric} model of INLA.
#'
#' @details This function considers an BYM2 CAR prior \insertCite{riebler2016intuitive}{bigDM} for the spatial latent effects of the different diseases and introduces correlation between them using the M-model proposal of \insertCite{botella2015unifying;textual}{bigDM}.
#' Putting the spatial latent effects for each disease in a matrix, the between disease dependence is introduced through the M matrix as \eqn{\Theta=\Phi M}, where the columns of \eqn{\Phi} follow an intrinsic BYM2 prior distribution (within-disease correlation).
#' A Wishart prior for the between covariance matrix \eqn{M'M} is considered using the Bartlett decomposition.
#' Uniform prior distributions on the interval [\code{alpha.min}, \code{alpha.max}] are considered for all the spatial mixing parameters.
#' The joint precision matrix of the spatial latent effects is given by \deqn{Q=(M^{-1}\otimes U)\operatorname{Blockdiag}(Q_1,\ldots,Q_J)(M^{-1}\otimes U)^\top,}
#' where \eqn{Q_j} denotes the diagonal BYM2 precision matrix in the spectral basis of the intrinsic CAR precision matrix, with \eqn{Q_{\mathrm{iCAR}}=U\Lambda U^\top}.
#' Since the resulting precision matrix is not sparse, this formulation is computationally more demanding than the corresponding CAR models.
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
#' @references
#' \insertRef{botella2015unifying}{bigDM}
#'
#' \insertRef{riebler2016intuitive}{bigDM}
#'
#' @param cmd Internal functions used by the \code{rgeneric} model to define the latent effect.
#' @param theta Vector of hyperparameters.
#'
#' @return This is used internally by the \code{INLA::inla.rgeneric.define()} function.
#'
#' @import Matrix
#' @importFrom utils getFromNamespace
#'
#' @seealso
#' \code{\link{Mmodel_icar}}, \code{\link{Mmodel_lcar}} and \code{\link{Mmodel_pcar}} for alternative multivariate CAR prior distributions.
#'
#' @export
########################################################################
## Mmodels - BYM2 (BARTLETT DECOMPOSITION)
########################################################################
Mmodel_bym2 <- function(cmd=c("graph","Q","mu","initial","log.norm.const","log.prior","quit"), theta=NULL){

  envir <- parent.env(environment())
  if(!exists("cache.done", envir=envir)){
          DW <- Matrix::Diagonal(x=colSums(W))-W
          U <- eigen(DW)$vectors
          # var.scale <- inla.ginv.diag(DW)
          var.scale <- exp(mean(log(getFromNamespace("inla.ginv.diag", "bigDM")(DW))))
          lambda <- eigen(DW)$values[-nrow(DW)]

          assign("U", U, envir = envir)
          assign("var.scale", var.scale, envir = envir)
          assign("lambda", lambda, envir = envir)
          assign("cache.done", TRUE, envir = envir)
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
          MI <- kronecker(M.inv, U)

          BlockIW <-
                  Matrix::bdiag(lapply(1:J, function(i) {
                          Matrix::Diagonal(x=c(var.scale*lambda/(param$alpha[i] + lambda*var.scale*(1-param$alpha[i])), 1/(1-param$alpha[i])))
                          }))
          Q <- (MI %*% BlockIW) %*% Matrix::t(MI)
          Q <- INLA::inla.as.sparse(Q)
          return (Q)
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


##################################################################################
## INLA-based function for computing marginal variances from a precision matrix ##
## i.e, diagonal elements of the Moore-Penrose inverse                          ##
##################################################################################
inla.ginv.diag <- function(Q, constr = NULL, eps = sqrt(.Machine$double.eps)) {

        marg.var <- rep(0, nrow(Q))
        Q <- INLA::inla.as.sparse(Q)
        g <- INLA::inla.read.graph(Q)

        if(is.null(constr)) constr <- list(A=matrix(1,1,nrow(Q)), e=0)

        for (k in seq_len(g$cc$n)) {
                i <- g$cc$nodes[[k]]
                n <- length(i)
                QQ <- Q[i, i, drop = FALSE]
                if (n == 1) {
                        QQ[1, 1] <- 1
                        marg.var[i] <- 1
                }
                else {
                        cconstr <- constr
                        if (!is.null(constr)) {
                                cconstr$A <- constr$A[, i, drop = FALSE]
                                eeps <- eps
                        }
                        else {
                                eeps <- 0
                        }
                        idx.zero <- which(rowSums(abs(cconstr$A)) == 0)
                        if (length(idx.zero) > 0) {
                                cconstr$A <- cconstr$A[-idx.zero, , drop = FALSE]
                                cconstr$e <- cconstr$e[-idx.zero]
                        }
                        res <- INLA::inla.qinv(QQ + Matrix::Diagonal(n) * max(diag(QQ)) * eeps, constr = cconstr)

                        # fac <- exp(mean(log(diag(res))))
                        # QQ <- fac * QQ
                        # marg.var[i] <- diag(res)/fac

                        marg.var[i] <- diag(res)
                }
                Q[i,i] <- QQ
        }
        return(marg.var)
}

utils::globalVariables(c("alpha.min","alpha.max"))
utils::globalVariables(c("J","W","interpret.theta","Q","log.prior","initial",
                         "dchisq","dnorm","initial.values"))
