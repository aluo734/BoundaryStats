#' @name simulate_rasters
#' @title Simulates rasters with boundary elements
#' @description Internal function that simulates a raster of the same extent and resolution as
#' the input raster, then defines boundaries
#' of the simulated raster.
#'
#' @param x A RasterLayer object.
#' @param model Neutral model to use. Options: 'random', 'gaussian_random' (Gaussian random field),
#' 'autoregressive' (simultaneous autoregressive model), and 'random_cluster' (modified random clusters method)
#' @param projection Numeric. EPSG code of input raster layer.
#'
#' @return A RasterLayer object with simulated values based on input data.
#'
#' @references
#' Saura, S. & Martínez-Millán, J. (2000) Landscape patterns simulation with a modified random clusters method. Landscape Ecology, 15, 661 – 678.
#' James, P. M. A., Fleming, R.A., & Fortin, M.-J. (2010) Identifying significant scale-specific spatial boundaries using wavelets and null models: Spruce budworm defoliation in Ontario, Canada as a case study. Landscape Ecology, 6, 873-887.
#'
#' @details This is an internal function.
simulate_rasters <- function (x, model, projection) {
  if (model == 'random') {
    #### RANDOM ####
    values <- x@data@values %>%
      na.omit(.) %>%
      sample(., length(.))
    x_sim <- x
    count = 1
    for (i in which(!is.na(x@data@values), arr.ind = T)) {
      x_sim@data@values[i] = values[count]
      count = count + 1
    }

  } else if (model == 'gaussian_random') {
    #### GAUSSIAN RANDOM FIELD ####
    # train variogram and model with input data
    input <- which(!is.na(x@data@values), arr.ind = T) %>%
      raster::xyFromCell(x, .) %>%
      cbind(na.omit(x@data@values)) %>%
      as.data.frame(.)
    colnames(input) <- c('lon', 'lat', 'value')
    sp::coordinates(input) <- ~ lon + lat

    vgm <- gstat::variogram(value ~ lon + lat, data = input) %>%
      gstat::fit.variogram(., model = gstat::vgm('Gau'), fit.kappa = T) # fit variogram to empirical data
    g.dummy <- gstat::gstat(formula = value ~ lon + lat, data = input, model = vgm, nmax = 100, dummy = T, beta = 0)

    # simulate new raster
    x_sim <- expand.grid(1:ncol(x), 1:nrow(x))
    names(x_sim) <- c('lon', 'lat')
    sp::gridded(x_sim) <- ~ lon + lat
    x_sim <- stats::predict(g.dummy, x_sim, nsim = 1) %>%
      as.data.frame(.) %>%
      raster::rasterFromXYZ(.)

    # change value range and crop spatial extent to match input data
    x_sim@data@values <- x_sim@data@values * diff(range(na.omit(x@data@values)))/diff(range(x_sim@data@values))
    a <- min(na.omit(x@data@values)) - min(x_sim@data@values)
    x_sim@data@values <- a + x_sim@data@values
    if (all(na.omit(x@data@values) %% 1 == 0)) {x_sim@data@values <- round(x_sim@data@values)}
    for (i in 1:length(x@data@values)) {
      if (is.na(x@data@values[i])) {x_sim@data@values[i] <- NA}
    }

  } else if (model == 'autoregressive') {
    #### SIMULTANEOUS AUTOREGRESSIVE MODEL ####
    # fit model parameters (intercept and rho) on input data
    input <- which(!is.na(x@data@values), arr.ind = T) %>%
      raster::xyFromCell(x, .) %>%
      cbind(na.omit(x@data@values)) %>%
      as.data.frame(.)
    colnames(input) <- c('lon', 'lat', 'value')
    sp::coordinates(input) <- ~ lon + lat
    adj <- spdep::cell2nb(nrow(x), ncol(x), 'queen') %>%
      spdep::subset.nb(., !is.na(x@data@values))
    listw <- spdep::nb2listw(adj, style = 'W')

    mod <- spatialreg::lagsarlm(value ~ 1, data = input, listw = listw)

    # use parameters to simulate new raster
    x_sim <- which(!is.na(x@data@values), arr.ind = T) %>%
      raster::xyFromCell(x, .) %>%
      cbind(., rep(NA,  nrow(.))) %>%
      as.data.frame(.)
    colnames(x_sim) <- c('lon', 'lat', 'value')
    sp::coordinates(x_sim) <- ~ lon + lat

    for (i in 1:nrow(x_sim@data)) {
      neigh <- adj[[i]][adj[[i]] < i]
      if (length(neigh) < 1 & length(na.omit(x_sim@data[neigh,])) < 1) {
        x_sim@data[i,] = sample(na.omit(x@data@values), 1)
      } else {
        lag <- mean(na.omit(x_sim@data[neigh,]))
        x_sim@data[i,] =mod$rho*lag + rnorm(1)
      }
    }

    x_sim <- raster::rasterFromXYZ(as.data.frame(x_sim))

    # if input values are all non-negative, make all simulated values non-negative
    if (min(na.omit(x@data@values)) == 0 & min(na.omit(x_sim@data@values)) < 0) {x_sim@data@values <- x_sim@data@values + abs(min(na.omit(x_sim@data@values)))}

    #i input values are all integers, make all simulated values integers
    if (all(na.omit(x@data@values) %% 1 == 0)) {x_sim@data@values <- round(x_sim@data@values)}

  } else if (model == 'random_cluster') {
    #### MODIFIED RANDOM CLUSTER ####
    # A: make percolated raster, where proportion of filled cells = p_cluster
    x_sim <- matrix(nrow = raster::nrow(x), ncol = raster::ncol(x))
    for (i in 1:length(x_sim)) {
      if (stats::runif(1) <= 0.5) {x_sim[i] = 1}
    }

    opp <- matrix(nrow = raster::nrow(x), ncol = raster::ncol(x)) # matrix with opposite cells marked for clumping
    opp[is.na(x_sim[])] <- 1

    # B: find clusters in raster
    x_sim <- raster::raster(x_sim, crs = paste0('+init=EPSG:', projection),
                            xmn = raster::extent(x)[1], xmx = raster::extent(x)[2],
                            ymn = raster::extent(x)[3], ymx = raster::extent(x)[4]) %>%
      raster::clump(., directions = 4)

    opp <- raster::raster(opp, crs = paste0('+init=EPSG:', projection),
                            xmn = raster::extent(x)[1], xmx = raster::extent(x)[2],
                            ymn = raster::extent(x)[3], ymx = raster::extent(x)[4]) %>%
      raster::clump(., directions = 4)
    opp@data@values <- opp@data@values + 0.5

    x_sim <- raster::merge(x_sim, opp)
    x_sim@data@values  <- raster::values(x) %>%
      na.omit(.) %>%
      max(.) + x_sim@data@values

    # C: assign clusters to categories
    prop <- raster::freq(x) %>% # proportions of cells in category
      as.data.frame(.) %>%
      subset(., value != 'NA' & value != 'NaN') %>%
      tibble::add_column(., p = .$count/sum(.$count)) %>%
      .[order(-.$p),]

    clump_selection_order <- sample(unique(x_sim@data@values))
    areas <- raster::area(x_sim)@data@values
    total_area <- sum(raster::area(x_sim)@data@values)
    row = 1
    for (i in clump_selection_order) {
      max_prop = prop[row, 3]; assignment = prop[row, 1]
      x_sim@data@values[x_sim@data@values == i] <- assignment

      current_prop <- which(x_sim@data@values == assignment) %>%
        areas[.] %>%
        sum(.)/total_area
      if (current_prop >= max_prop) {row = row + 1}

      if (row > nrow(prop)) {break}
    }

    # # D: fill in any remaining unclassified cells
    # # part of process from Saura and Martinex-Millan, but the previous steps as coded don't
    # # leave unclassified cells
    # for (i in 1:terra::ncell(x_sim)) {
    #
    #   if (is.na(x_sim[i])) {
    #     neighbors <- raster::adjacent(x_sim, i, directions = 8, pairs = F) %>%
    #       x_sim[.] %>%
    #       na.omit()
    #
    #     if (length(neighbors > 0)) {
    #       u <- unique(neighbors)
    #       x_sim[i] <- match(neighbors, u) %>%
    #         tabulate(.) %>%
    #         which.max(.) %>%
    #         u[.]
    #     } else {
    #       a <- prop[,c(1,3)]
    #       colnames(a) <- c('x', 'prob')
    #       a <- pdqr::new_r(a, 'discrete')
    #       x_sim[i] <- a(1)
    #     }
    #
    #   }
    # }

    # crop extent of filled values to input data range
    for (i in 1:length(x@data@values)) {
      if (is.na(x@data@values[i])) {x_sim@data@values[i] <- NA}
    }

  } else {
    stop('null model not available')
  }

  raster::extent(x_sim) <- raster::extent(x)
  raster::crs(x_sim) <- paste0('+init=EPSG:', projection)

  return(x_sim)

}
