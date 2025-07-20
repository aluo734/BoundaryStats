
# BoundaryStats

BoundaryStats was designed to test for the presence of geographic
boundaries in ecological variables and overlap between such boundaries.
Users can calculate boundary and boundary overlap statistics with raster
data. BoundaryStats can create null distributions for the statistics
based on various neutral landscape models that are parameterized on the
empirical data. The primary functions are statistical tests for the
presence of spatial boundaries of a variable and significant overlap
between the spatial boundaries of two variables.

[![DOI](https://zenodo.org/badge/534683960.svg)](https://zenodo.org/badge/latestdoi/534683960)

## Installation

You can install BoundaryStats with either:

``` r
install.packages('BoundaryStats')
remotes::install_github("aluo734/BoundaryStats")
```

## Statistical Tests

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">

<thead>

<tr>

<th style="text-align:left;">

Function
</th>

<th style="text-align:left;">

Category
</th>

<th style="text-align:left;">

Description
</th>

</tr>

</thead>

<tbody>

<tr>

<td style="text-align:left;">

n_boundaries
</td>

<td style="text-align:left;">

Boundary
</td>

<td style="text-align:left;">

The number of subgraphs, or sets of contiguous boundary elements, in the
data.
</td>

</tr>

<tr>

<td style="text-align:left;">

longest_boundary
</td>

<td style="text-align:left;">

Boundary
</td>

<td style="text-align:left;">

The length of the longest subgraph.
</td>

</tr>

<tr>

<td style="text-align:left;">

n_overlap_boundaries
</td>

<td style="text-align:left;">

Boundary Overlap
</td>

<td style="text-align:left;">

The number of directly overlapping boundary elements, or raster cells
labeled as part of a boundary, of two traits.
</td>

</tr>

<tr>

<td style="text-align:left;">

average_min_x_to_y
</td>

<td style="text-align:left;">

Boundary Overlap
</td>

<td style="text-align:left;">

The average minimum distance between each boundary element in raster x
and the nearest boundary element in raster y. Uses Euclidean distance.
The boundaries of trait x depend on the boundaries of trait y.
</td>

</tr>

<tr>

<td style="text-align:left;">

average_min_distance
</td>

<td style="text-align:left;">

Boundary Overlap
</td>

<td style="text-align:left;">

The average minimum distance between boundary elements in two raster
layers. Uses Euclidean distance. Boundaries for each trait affect one
another reciprocally (x affects y and y affects x).
</td>

</tr>

</tbody>

</table>

## Example

``` r
library(BoundaryStats)
library(tidyverse)

data(ecoregions)
ecoregions <- terra::rast(ecoregions_matrix, crs = ecoregions_crs)
terra::ext(ecoregions) <- ecoregions_ext

data(L.flavomaculatus)
L.flavomaculatus <- terra::rast(L.flavomaculatus_matrix, crs = L.flavomaculatus_crs)
terra::ext(L.flavomaculatus) <- L.flavomaculatus_ext

terra::crs(ecoregions) <- terra::crs(L.flavomaculatus)
ecoregions <- terra::resample(ecoregions, L.flavomaculatus) |>
  terra::crop(L.flavomaculatus) |>
  terra::mask(L.flavomaculatus)
L.flavomaculatus <- terra::crop(L.flavomaculatus, ecoregions) |>
  terra::mask(ecoregions)

L.flavomaculatus_boundary <- define_boundary(L.flavomaculatus, cat = TRUE)
ecoregions_boundary <- define_boundary(ecoregions, cat = T)
plot_boundary(L.flavomaculatus, ecoregions)

Tcrist_ovlp_null <- overlap_null_distrib(L.flavomaculatus, ecoregions, rand_both = FALSE, x_cat = T, n_iterations = 100, x_model = 'random_cluster')

n_overlap_boundaries(L.flavomaculatus, ecoregions, Tcrist_ovlp_null)
average_min_x_to_y(L.flavomaculatus, ecoregions, Tcrist_ovlp_null)
average_min_distance(L.flavomaculatus, ecoregions, Tcrist_ovlp_null)
```

Data source: Cox, Karen; Schepers, Robbert; Van Breusegem, An;
Speybroeck, Jeroen (2023), The common ground in landscape effects on
gene flow in two newt species in an agroecosystem, Dryad, Dataset,
<https://doi.org/10.5061/dryad.bk3j9kdhz>.
