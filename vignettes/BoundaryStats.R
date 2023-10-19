## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message = FALSE----------------------------------------------------------
library(BoundaryStats)
library(terra)
library(magrittr)

## ----fig.width = 8, fig.height = 6--------------------------------------------
data(ecoregions)
ecoregions <- terra::rast(ecoregions_matrix, crs = ecoregions_crs)
ext(ecoregions) <- ecoregions_ext
plot(ecoregions)

## ----fig.width = 8, fig.height = 6--------------------------------------------
data(L.flavomaculatus)
L.flavomaculatus <- terra::rast(L.flavomaculatus_matrix, crs = L.flavomaculatus_crs)
ext(L.flavomaculatus) <- L.flavomaculatus_ext
plot(L.flavomaculatus)

## -----------------------------------------------------------------------------
crs(ecoregions) <- crs(L.flavomaculatus)
ecoregions <- resample(ecoregions, L.flavomaculatus) %>%
  crop(., L.flavomaculatus) %>%
  mask(., L.flavomaculatus)
L.flavomaculatus <- crop(L.flavomaculatus, ecoregions) %>%
  mask(., ecoregions)

## ----fig.width = 8, fig.height = 6--------------------------------------------
ecoregions_boundaries <- categorical_boundary(ecoregions)
L.flavomaculatus_boundaries <- categorical_boundary(L.flavomaculatus)

## ----warning = FALSE, fig.width = 8, fig.height = 6---------------------------
plot_boundary(L.flavomaculatus_boundaries, ecoregions_boundaries, trait_names = c('A. delicatus genetic group', 'Ecoregion'))

