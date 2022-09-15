
# BoundaryStats

BoundaryStats was designed to test for the presence of geographic
boundaries in ecological variables and overlap between such boundaries.
Users can calculate boundary and boundary overlap statistics with raster
data. BoundaryStats can create null distributions for the statistics
based on various neutral landscape models that are parameterized on the
empirical data. The primary functions are statistical tests for the
presence of spatial boundaries of a variable and significant overlap
between the spatial boundaries of two variables.

## Installation

You can install BoundaryStats with:

``` r
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
n_subgraph
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
max_subgraph
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
Odirect
</td>
<td style="text-align:left;">
Boundary Overlap
</td>
<td style="text-align:left;">
The number of directly overlapping boundary elements, or raster cells
labelled as part of a boundary, of two traits.
</td>
</tr>
<tr>
<td style="text-align:left;">
Ox
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
Oxy
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
library(raster)

data(ecoregions)
data(A.delicatus)

A.delicatus_ecoregion <- resample(ecoregions, A.delicatus) %>%
  crop(., A.delicatus) %>%
  mask(., A.delicatus) %>%
  categorical_boundary(., projection = 4210)
A.delicatus <- crop(A.delicatus, A.delicatus_ecoregion) %>%
  mask(., A.delicatus_ecoregion)
A.delicatus <- categorical_boundary(A.delicatus, projection = 4210)

overlap.null <- overlap_null_distrib(A.delicatus, A.delicatus_ecoregion, x_cat = T, y_cat = T, n_iterations = 100, projection = 4210, x_model = 'random_cluster', y_model = 'random_cluster')

Odirect(A.delicatus, A.delicatus_ecoregion, overlap.null)
Ox(A.delicatus, A.delicatus_ecoregion, overlap.null)
Oxy(A.delicatus, A.delicatus_ecoregion, overlap.null)
```

Data from: Barratt, C.D., Bwong, B.A., Jehle, R., Liedtke, C.H., Nagel,
P., Onstein, R.E., Portik, D.M., Streicher, J.W. & Loader, S.P. (2018)
Vanishing refuge? Testing the forest refuge hypothesis in coastal East
Africa using genome‐wide sequence data for seven amphibians. Molecular
Ecology, 27, 4289-4308.
