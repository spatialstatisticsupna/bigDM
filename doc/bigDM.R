## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval=FALSE, include=TRUE-------------------------------------------------
#  # Install devtools package from CRAN repository
#  install.packages("devtools")
#  
#  # Load devtools library
#  library(devtools)
#  
#  # Install the R-INLA package
#  install.packages("INLA", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#  
#  # Install bigDM from GitHub repositoy
#  install_github("spatialstatisticsupna/bigDM")

## -----------------------------------------------------------------------------
library(bigDM)

head(Carto_SpainMUN)

## ---- message=FALSE, warning=FALSE, class.source="indent"---------------------
library(rgdal)
  
carto_spdf <- readOGR(system.file("shape/Carto_SpainMUN.shp", package="bigDM"))

## ----message=FALSE, warning=FALSE, class.source="indent"----------------------
library(sf)

carto_sf <- st_read(system.file("shape/Carto_SpainMUN.shp", package="bigDM"))

## ---- class.source="indent"---------------------------------------------------
library(readr)

data <- read_csv(system.file("csv/data.csv", package="bigDM"))
carto <- merge(carto_sf, data, by="ID")

head(carto)

## ---- eval=FALSE--------------------------------------------------------------
#  ## Not run:
#  Global <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", O="obs", E="exp", model="global")

## -----------------------------------------------------------------------------
carto.new <- random_partition(carto=Carto_SpainMUN, rows=3, columns=3)

## ---- label="randomPartition", fig.cap=paste("Random partition for the Spanish municipalities based on a 3x3 regular grid."), fig.height=4.5----
library(tmap)

tm_shape(carto.new) +
  tm_polygons(col="ID.group") +
  tm_layout(legend.outside=TRUE)

## ----warning=FALSE------------------------------------------------------------
Disjoint <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
                     model="partition", k=0, seed=1234, strategy="gaussian", save.models=TRUE)

summary(Disjoint)

## ----warning=FALSE------------------------------------------------------------
## 1st order neighbourhood model ##
order1 <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
                   model="partition", k=1, seed=1234, strategy="gaussian", save.models=TRUE)
summary(order1)

## 2nd order neighbourhood model ##
order2 <- CAR_INLA(carto=Carto_SpainMUN, ID.area="ID", ID.group="region", O="obs", E="exp",
                   model="partition", k=2, seed=1234, strategy="gaussian", save.models=TRUE)
summary(order2)

## ----fig.cap=paste("Disjoint partition for the Spanish municipalities based on the $D=15$ autonomous regions."), fig.height=4.5, label="disjointPartition", echo=FALSE----
tm_shape(Carto_SpainMUN) +
  tm_polygons(col="region") +
  tm_layout(legend.outside=TRUE)

## -----------------------------------------------------------------------------
## disjoint partition ##
carto.k0 <- divide_carto(carto=Carto_SpainMUN, ID.group="region", k=0)

## 1st order partition ##
carto.k1 <- divide_carto(carto=Carto_SpainMUN, ID.group="region", k=1)

## 2nd order partition ##
carto.k2 <- divide_carto(carto=Carto_SpainMUN, ID.group="region", k=2)

## ----fig.cap=paste("Spatial polygons for the autonomous region of Castilla y León with 1st and 2nd order neighbours in the border."), fig.height=4.5, label="korderPartition", echo=FALSE----
par(mai=c(0,0.1,0.3,0.1))
plot(carto.k2$`Castilla y Leon`$geometry, col="dodgerblue4", main="Castilla-León")
plot(carto.k1$`Castilla y Leon`$geometry, col="dodgerblue", add=TRUE)
plot(carto.k0$`Castilla y Leon`$geometry, col="lightgrey", add=TRUE)

## ----include=FALSE------------------------------------------------------------
files <- list.files("temp/", full.names=TRUE)
file.rename(files[1],"temp/INLAsubmodels_Disjoint.Rdata")
file.rename(files[2],"temp/INLAsubmodels_order1.Rdata")
file.rename(files[3],"temp/INLAsubmodels_order2.Rdata")

## ----warning=FALSE------------------------------------------------------------
## Disjoint model ##
load("temp/INLAsubmodels_Disjoint.Rdata")

Disjoint <- mergeINLA(inla.models=inla.models, k=0, seed=1234, n.sample=1000,
                      compute.fixed=TRUE, compute.DIC=TRUE)

summary(Disjoint)

## 1st order neighbourhood model ##
load("temp/INLAsubmodels_order1.Rdata")

order1 <- mergeINLA(inla.models=inla.models, k=1, seed=1234, n.sample=1000,
                      compute.fixed=FALSE, compute.DIC=TRUE)

summary(order1)

## 2nd order neighbourhood model ##
load("temp/INLAsubmodels_order2.Rdata")

order2 <- mergeINLA(inla.models=inla.models, k=2, seed=1234, n.sample=1000,
                      compute.fixed=FALSE, compute.DIC=TRUE)

summary(order2)

