## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message = FALSE----------------------------------------------------------
library(BoundaryStats)
library(magrittr)
library(raster)

## ----fig.width = 8, fig.height = 6--------------------------------------------
data(ecoregions)
plot(ecoregions)

## ----fig.width = 8, fig.height = 6--------------------------------------------
data(A.delicatus)
plot(A.delicatus)

## ----fig.width = 8, fig.height = 6--------------------------------------------
A.delicatus_ecoregion <- resample(ecoregions, A.delicatus) %>%
  crop(., A.delicatus) %>%
  mask(., A.delicatus) %>%
  categorical_boundary(., projection = 4210)
A.delicatus <- crop(A.delicatus, A.delicatus_ecoregion) %>%
  mask(., A.delicatus_ecoregion)
A.delicatus <- categorical_boundary(A.delicatus, projection = 4210)
plot(A.delicatus)

