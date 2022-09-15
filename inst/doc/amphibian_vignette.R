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

## ----warning = FALSE, fig.width = 8, fig.height = 6---------------------------
plot_boundary(A.delicatus, A.delicatus_ecoregion, trait_names = c('A. delicatus genetic group', 'Ecoregion'))

## -----------------------------------------------------------------------------
A.deli_bound.null <- boundary_null_distrib(A.delicatus, cat = T, n_iterations = 100, projection = 4210, model = 'random_cluster', progress = F)

## -----------------------------------------------------------------------------
n_subgraph(A.delicatus, A.deli_bound.null)
max_subgraph(A.delicatus, A.deli_bound.null, projection = 4210)

## -----------------------------------------------------------------------------
A.deli_overlap.null <- overlap_null_distrib(A.delicatus, A.delicatus_ecoregion, x_cat = T, y_cat = T, n_iterations = 100, projection = 4210, x_model = 'random_cluster', y_model = 'random_cluster', progress = F)

## -----------------------------------------------------------------------------
Odirect(A.delicatus, A.delicatus_ecoregion, A.deli_overlap.null)
Ox(A.delicatus, A.delicatus_ecoregion, A.deli_overlap.null)
Oxy(A.delicatus, A.delicatus_ecoregion, A.deli_overlap.null)

