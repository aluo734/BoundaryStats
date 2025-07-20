## version 2.3.0

---

### Updates

- length of boundaries is calculated as distance from farthest points through edges of a subgraph of adjacent boundary elements, instead of straight line between points
- Gaussian neutral model uses autocorrelation range based on variogram, rather than LISA clusters
- categorical_boundary function deprecated, and functionality merged with define_boundary function
- renamed statistical test functions to improve clarity
- plot_boundary now optionally returns a SpatRaster object with boundary elements and boundary element overlaps
- removed dependencies on sf and pdqr
- added dependency on gstat

## version 2.2.0

---

### Updates

- removed dependency on rgeoda

## version 2.1.0

---

### Updates
- reduced the number of iterations in the vignette to save checktime

## version 2.1.0

---


### Updates
- fixed algorithm for modified random clusters
- fixed bug for overlap null distribution function
- updated vignette to match above changes

## version 2.0.1

---


### Updates

- Edits to upload to CRAN repository


## version 2.0.0

---


### Updates

- removed dependencies on raster and sp, in anticipation of the retirement of rgeos and rgdal
- improved computational efficiency for null models
- removed simultaneous autoregressive model option for null models
- updated vignette and examples