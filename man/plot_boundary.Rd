% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_boundary.R
\name{plot_boundary}
\alias{plot_boundary}
\title{Map the boundary elements of two raster layers}
\usage{
plot_boundary(x, y, color = NA, trait_names = NA)
}
\arguments{
\item{x}{A RasterLayer object with boundary elements.}

\item{y}{A RasterLayer object with boundary elements.}

\item{color}{Optional. A character vector of up to three colors (x boundary, y boundary, and overlapping elements).}

\item{trait_names}{Optional. A character vector with up to two elements (legend name for x and legend name for y).}
}
\value{
A ggplot object.
}
\description{
This is a wrapper function for ggplot2 that will produce a map of boundary
elements for two traits and show where boundary elements intersect.
}
\examples{
plot_boundary(soil_boundary, genetic_boundary, col = c('pink', 'orange', 'red'),
              trait_names = c('Soil Type', 'Genetic Group'))

}
\author{
Amy Luo
}
