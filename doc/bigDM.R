## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE, warning=FALSE---------------------------------------------
#### import to spatialpolygonsdataframe ####
library(rgdal)
zones_spdf <-
  readOGR(system.file("shape/Carto_SpainMUN.shp", package = "bigDM"),
          layer = "Carto_SpainMUN")

## ----message=FALSE, warning=FALSE---------------------------------------------
#### import to sf####
library(sf)
zones_sf <-
  st_read(system.file("shape/Carto_SpainMUN.shp", package = "bigDM"))

## ----Merge data and polygons--------------------------------------------------
#### data with expected and observed cases in csv file ####
library(readr)
data <- read_csv(system.file("csv/data.csv", package = "bigDM"), 
                 col_types = cols(ID = col_character(),
                                  exp = col_number(), 
                                  obs = col_number()))
zones_sf<- merge (zones_sf, data, by.x="ID", by.y="ID")
head(zones_sf)

## ----Data and partition-------------------------------------------------------
#### libraries and data####
library(bigDM)
data("Carto_SpainMUN")
head(Carto_SpainMUN)

#### disjoint partition ####
carto.d <- divide_carto(carto=Carto_SpainMUN, ID.group="region", k=0)
data.d <- lapply(carto.d, function(x) sf::st_set_geometry(x, NULL))

## ----fig.align='center', fig.width=7,fig.height=5,fig.retina=2----------------
plot(carto.d[[14]]$geometry, main = "Navarre")

## ----Q matrix, message=FALSE, results=FALSE-----------------------------------
#### list of object’s neighbours and spatial neighbourhood matrix ####
library(Matrix)
aux <- lapply(carto.d, function(x) connect_subgraphs(x, ID.area = "ID"))
Wd <- lapply(aux, function(x) x$W)
nd <- lapply(Wd, function(x) nrow(x))
Rd <- mapply(function(x,y){Diagonal(x,colSums(y))-y}, x=nd, y=Wd)
Rd.Leroux <- mapply(function(x,y){Diagonal(x)-y}, x=nd, y=Rd)

## ----Data Partition-----------------------------------------------------------
#### List of datasets for INLA ####
data.INLA <-
  mapply(function(x, y) {
    data.frame(
      O = x[, "obs"],
      E = x[, "exp"],
      Area = x[, "ID"],
      ID.area = seq(1, y)
    )
  },
  x = data.d,
  y = nd,
  SIMPLIFY = FALSE)
D <- length(data.INLA)

## ----Priors-------------------------------------------------------------------
#### Priors distributions of hiperparameters ####
sdunif = "expression: logdens=-log_precision/2; return(logdens)"
lunif = "expression: a = 1; b = 1; beta = exp(theta)/(1+exp(theta));
        logdens = lgamma(a+b)-lgamma(a)-lgamma(b)+(a-1)*log(beta)+(b-1)*log(1-beta);
        log_jacobian = log(beta*(1-beta));
        return(logdens+log_jacobian)"

## ----Fit models disjoint partitions, message=FALSE, results=FALSE-------------
#### Fit of the models ####
library(INLA)
Models <- vector("list",D)
for(i in 1:D) {
  cat(sprintf("+ Model %d of %d", i, D), "\n")
  Rs <- Rd[[i]]
  Rs.Leroux <- Rd.Leroux[[i]]
  formula <- O ~ f(
    ID.area,
    model = "generic1",
    Cmatrix = Rs.Leroux,
    constr = TRUE,
    hyper = list(
      prec = list(prior = sdunif),
      beta = list(prior = lunif, initial = 0)
    )
  )
  Models[[i]] <-
    inla(
      formula,
      family = "poisson",
      data = data.INLA[[i]],
      E = E,
      control.predictor = list(compute = TRUE,
                               cdf = c(log(1))),
      control.compute = list(
        dic = TRUE,
        cpo = TRUE,
        waic = TRUE,
        config = TRUE
      ),
      control.inla = list(strategy = "gaussian")
    )
}

## ----Fit disjoint MI no alpha-------------------------------------------------
#### Merge models ####
disjoint <-
  mergeINLA(
    inla.models = Models,
    k=0,
    ID.area = "Area",
    O = "O",
    E = "E",
    seed=1234,
    compute.fixed = FALSE,
    compute.DIC = TRUE
  )

#### DIC ####
disjoint$dic$dic

#### WAIC ####
disjoint$waic$waic

## ----Fit disjoint MI alpha, eval=FALSE, echo=TRUE-----------------------------
#  #### Merge models with \alpha ####
#  disjoint <- mergeINLA(
#    inla.models = Models,
#    k = 0,
#    ID.area = "muni",
#    O = "O",
#    E = "E",
#    seed = 1234,
#    n.sample = 1000,
#    compute.fixed = TRUE,
#    compute.DIC = TRUE
#  )
#  summary(disjoint)

## ----Fit disjoint CI----------------------------------------------------------
disjoint_1 <- CAR_INLA(
  carto = Carto_SpainMUN,
  ID.area = "ID",
  ID.group = "region",
  O = "obs",
  E = "exp",
  prior = "Leroux",
  model = "partition",
  k = 0,
  seed = 1234,
  strategy = "gaussian",
  PCpriors = FALSE,
  compute.fixed = FALSE,
  compute.DIC = TRUE,
  save.models = FALSE
)

#### DIC ####
disjoint_1$dic$dic

#### WAIC ####
disjoint_1$waic$waic

## ----Carto korder, fig.align='center', fig.height=5, fig.retina=2, fig.width=7----
#### 1st order neighbours ####
a <-
  divide_carto(
    carto = Carto_SpainMUN,
    ID.group = "region",
    k = 1,
    plot = FALSE
  )

#### 2nd order neighbours ####
b <-
  divide_carto(
    carto = Carto_SpainMUN,
    ID.group = "region",
    k = 2,
    plot = FALSE
  )

#### 3rd order neighbours ####
c <-
  divide_carto(
    carto = Carto_SpainMUN,
    ID.group = "region",
    k = 3,
    plot = FALSE
  )

#### Navarre with 1st, 2nd and 3rd order neighbours ####
plot(c[[14]]$geometry, border = "blue",main="Navarre with 1st, 2nd and 3rd order neighbours")
plot(b[[14]]$geometry, border = "red", add = TRUE)
plot(a[[14]]$geometry, border = "green", add = TRUE)
plot(carto.d[[14]]$geometry, add = TRUE)

## ----Fit korder models, eval=FALSE, echo=TRUE---------------------------------
#  #### 1st order neighbours model ####
#  order1<-CAR_INLA(
#  carto = Carto_SpainMUN,
#  ID.area = "ID",
#  ID.group = "region",
#  O = "obs",
#  E = "exp",
#  prior = "Leroux",
#  model = "partition",
#  k = 1,
#  strategy = "gaussian")
#  
#  #### 2nd order neighbours model ####
#  order2<-CAR_INLA(
#  carto = Carto_SpainMUN,
#  ID.area = "ID",
#  ID.group = "region",
#  O = "obs",
#  E = "exp",
#  prior = "Leroux",
#  model = "partition",
#  k = 2,
#  strategy = "gaussian")
#  
#  #### 3rd order neighbours model ####
#  order3<-CAR_INLA(
#  carto = Carto_SpainMUN,
#  ID.area = "ID",
#  ID.group = "region",
#  O = "obs",
#  E = "exp",
#  prior = "Leroux",
#  model = "partition",
#  k = 3,
#  strategy = "gaussian")

## ----Random Partition, fig.align='center', fig.height=5, fig.retina=2, fig.width=7----
library(tmap)

## define a random partition based on a 2x2 regular grid ##
carto_random <- random_partition(carto=Carto_SpainMUN, rows=2, columns=2)

## plot of the grouping variable 'ID.group' ## 
tm_shape(carto_random) +
tm_polygons(col="ID.group") + 
  tm_layout(legend.outside=TRUE)

