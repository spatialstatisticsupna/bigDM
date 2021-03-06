% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CAR_INLA.R
\name{CAR_INLA}
\alias{CAR_INLA}
\title{Fit a (scalable) spatial Poisson mixed model to areal count data, where several CAR prior distributions can be specified for the spatial random effect.}
\usage{
CAR_INLA(
  carto = NULL,
  ID.area = NULL,
  ID.group = NULL,
  O = NULL,
  E = NULL,
  W = NULL,
  prior = "Leroux",
  model = "partition",
  k = 0,
  strategy = "simplified.laplace",
  PCpriors = FALSE,
  seed = NULL,
  n.sample = 1000,
  compute.fixed = FALSE,
  compute.DIC = TRUE,
  save.models = FALSE,
  plan = "sequential",
  workers = NULL
)
}
\arguments{
\item{carto}{object of class \code{SpatialPolygonsDataFrame} or \code{sf}.
This object must contain at least the target variables of interest specified in the arguments \code{ID.area}, \code{O} and \code{E}.}

\item{ID.area}{character; name of the variable which contains the IDs of spatial areal units.}

\item{ID.group}{character; name of the variable which contains the IDs of the spatial partition (grouping variable).
Only required if \code{model="partition"}.}

\item{O}{character; name of the variable which contains the observed number of disease cases for each areal units.}

\item{E}{character; name of the variable which contains either the expected number of disease cases or the population at risk for each areal unit.}

\item{W}{optional argument with the binary adjacency matrix of the spatial areal units.  If \code{NULL} (default), this object is computed from the \code{carto} argument (two areas are considered as neighbours if they share a common border).}

\item{prior}{one of either \code{"Leroux"} (default), \code{"intrinsic"}, \code{"BYM"} or \code{"BYM2"},
which specifies the prior distribution considered for the spatial random effect.}

\item{model}{one of either \code{"global"} or \code{"partition"} (default), which specifies the \emph{Global model}
or one of the scalable model proposal's (\emph{Disjoint model} and \emph{k-order neighbourhood model}, respectively).}

\item{k}{numeric value with the neighbourhood order used for the partition model. Usually k=2 or 3 is enough to get good results.
If k=0 (default) the \emph{Disjoint model} is considered. Only required if \code{model="partition"}.}

\item{strategy}{one of either \code{"gaussian"}, \code{"simplified.laplace"} (default), \code{"laplace"} or \code{"adaptive"},
which specifies the approximation strategy considered in the \code{inla} function.}

\item{PCpriors}{logical value (default \code{FALSE}); if \code{TRUE} then penalised complexity (PC) priors are used for the precision parameter of the spatial random effect.
Only works if arguments \code{prior="intrinsic"} or \code{prior="BYM2"} are specified.}

\item{seed}{numeric (default \code{NULL}); control the RNG of the \code{inla.qsample} function. See \code{help(inla.qsample)} for further information.}

\item{n.sample}{numeric; number of samples to generate from the posterior marginal distribution of the risks. Default to 1000.}

\item{compute.fixed}{logical value (default \code{FALSE}); if \code{TRUE} then the overall log-risk \eqn{\alpha} is computed.
Only works if \code{k=0} argument (\emph{disjoint model}) is specified.}

\item{compute.DIC}{logical value; if \code{TRUE} (default) then approximate values of the Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) are computed.}

\item{save.models}{logical value (default \code{FALSE}); if \code{TRUE} then a list with all the \code{inla} submodels is saved in '/temp/' folder, which can be used as input argument for the \code{\link{mergeINLA}} function.}

\item{plan}{one of either \code{"sequential"} or \code{"cluster"}, which specifies the computation strategy used for model fitting using the 'future' package.
If \code{plan="sequential"} (default) the models are fitted sequentially and in the current R session (local machine). If \code{plan="cluster"} the models are fitted in parallel on external R sessions (local machine) or distributed in remote compute nodes.}

\item{workers}{character or vector (default \code{NULL}) containing the identifications of the local or remote workers where the models are going to be processed. Only required if \code{plan="cluster"}.}
}
\value{
This function returns an object of class \code{inla}. See the \code{\link{mergeINLA}} function for details.
}
\description{
Fit a spatial Poisson mixed model to areal count data. The linear predictor is modelled as the sum of a global intercept and a spatially structured random effect.
For the latter, several conditional autoregressive (CAR) prior distributions can be specified, such as the intrinsic CAR prior \insertCite{besag1991}{bigDM}, the convolution or BYM prior \insertCite{besag1991}{bigDM},
the CAR prior proposed by \insertCite{leroux1999estimation;textual}{bigDM}, and the reparameterization of the BYM model given by \insertCite{dean2001detecting;textual}{bigDM} named BYM2.
\cr\cr
Three main modeling approaches can be considered:
\itemize{
\item the usual model with a global spatial random effect whose dependence structure is based on the whole neighbourhood graph of the areal units (\code{model="global"} argument)
\item a disjoint model based on a partition of the whole spatial domain where independent spatial CAR models are simultaneously fitted in each partition (\code{model="partition"} and \code{k=0} arguments)
\item a modelling approach where \emph{k}-order neighbours are added to each partition to avoid border effects in the disjoint model (\code{model="partition"} and \code{k>0} arguments).
}
For both the disjoint and k-order neighbour models, parallel or distributed computation strategies can be performed to speed up computations by using the 'future' package \insertCite{bengtsson2020unifying}{bigDM}.

Inference is conducted in a fully Bayesian setting using the integrated nested Laplace approximation (INLA; \insertCite{rue2009approximate;textual}{bigDM}) technique through the R-INLA package (\url{http://www.r-inla.org/}).
For the scalable model proposals \insertCite{orozco2020}{bigDM}, approximate values of the Deviance Information Criterion (DIC) and Watanabe-Akaike Information Criterion (WAIC) can also be computed.
}
\details{
For a full model specification and further details see the vignettes accompanying this package.
}
\examples{
## load the Spain colorectal cancer mortality data ##
data(Carto_SpainMUN)

\dontrun{
## fit the global model with a Leroux CAR prior distribution ##
Global <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", O="obs", E="exp",
                   prior="Leroux", model="global", strategy="gaussian")

summary(Global)

## fit the disjoint model with a Leroux CAR prior distribution ##
Disjoint <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
                     prior="Leroux", model="partition", k=0, strategy="gaussian")
summary(Disjoint)

## fit the 1st order neighbourhood model with a Leroux CAR prior distribution ##
order1 <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
                   prior="Leroux", model="partition", k=1, strategy="gaussian")
summary(order1)

## fit the 2nd order neighbourhood model with a Leroux CAR prior distribution ##
order2 <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
                   prior="Leroux", model="partition", k=2, strategy="gaussian")
summary(order2)
}

}
\references{
\insertRef{bengtsson2020unifying}{bigDM}

\insertRef{besag1991}{bigDM}

\insertRef{dean2001detecting}{bigDM}

\insertRef{leroux1999estimation}{bigDM}

\insertRef{rue2009approximate}{bigDM}

\insertRef{orozco2020}{bigDM}
}
