# number of subgraphs ----
#' @name n_boundaries
#' @title Number of boundaries
#' @description Statistical test the for number of subgraphs, or sets of contiguous boundary elements, in the data.
#'
#' @param x A SpatRaster object with boundary elements.
#' @param null_distrib A list of probability functions output from boundary_null_distrib().
#'
#' @return The number of subgraphs in the raster and a p-value.
#' @examples \donttest{
#' data(T.cristatus)
#' T.cristatus <- terra::rast(T.cristatus_matrix, crs = T.cristatus_crs)
#' terra::ext(T.cristatus) <- T.cristatus_ext
#'
#' T.crist_boundaries <- define_boundary(T.cristatus, cat = TRUE)
#' T.crist_bound_null <- boundary_null_distrib(T.cristatus, cat = TRUE,
#' n_iterations = 100, model = 'random_cluster')
#'
#' n_boundaries(T.crist_boundaries, T.crist_bound_null)
#' }
#' 
#' @author Amy Luo
#' @references Jacquez, G.M., Maruca,I S. & Fortin M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' @export
n_boundaries <- function(x, null_distrib) {
  n_boundary <- terra::patches(x, directions = 8, zeroAsNA = TRUE) %>%
    terra::values() %>%
    max(., na.rm = TRUE)

  p <- null_distrib$n_boundary(n_boundary) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)

  names(n_boundary) <- 'n_boundary'
  names(p) <-'p-value'
  return(c(n_boundary, p))
}

#length of the longest subgraph ----
#' @name longest_boundary
#' @title Length of the longest boundary
#' @description Statistical test for the length of the longest subgraph, or set of contiguous boundary elements.
#'
#' @param x A SpatRaster object with boundary elements.
#' @param null_distrib A list of probability functions output from boundary_null_distrib().
#'
#' @return The length of the longest subgraph and a p-value.
#' @examples \donttest{
#' data(T.cristatus)
#' T.cristatus <- terra::rast(T.cristatus_matrix, crs = T.cristatus_crs)
#' terra::ext(T.cristatus) <- T.cristatus_ext
#'
#' Tcrist_boundaries <- define_boundary(T.cristatus, cat = FALSE)
#' T.crist_bound_null <- boundary_null_distrib(T.cristatus, cat = TRUE, n_iterations = 100,
#'   model = 'random_cluster')
#'
#' longest_boundary(Tcrist_boundaries, T.crist_bound_null)
#' }
#' 
#' @author Amy Luo
#' @references Jacquez, G.M., Maruca,I S. & Fortin M.-J. (2000) From fields to objects: A review of geographic boundary analysis. Journal of Geographical Systems, 3, 221, 241.
#' @export
longest_boundary <- function(x, null_distrib) {
  
  # organize boundary elements into separate boundaries based on which boundary elements are adjacent
  boundary_patches <- terra::patches(x, directions = 8, zeroAsNA = TRUE)
  patch_ids <- terra::values(boundary_patches) %>%
    unique %>%
    na.omit %>%
    as.vector

  # for each boundary/subgraph
  lengths <- c()
  for (i in patch_ids) {
    # find the farthest two cells in the boundary
    cell_coord <- terra::xyFromCell(boundary_patches, terra::cells(boundary_patches, i)$patches)
    rownames(cell_coord) <- terra::cells(boundary_patches, i)$patches
    dist.matrix <- as.matrix(terra::distance(cell_coord, lonlat = T))
    rownames(dist.matrix) <- rownames(cell_coord)
    colnames(dist.matrix) <- rownames(cell_coord)
    far_points <- which(dist.matrix == max(dist.matrix), arr.ind = T)[1,] %>%
      as.vector %>%
      cell_coord[., ] %>%
      rownames
    
    # convert the raster cells into a subgraph with nodes, edges, and edge weights (weights are geographic distance)
    neighbor_cells <- terra::adjacent(boundary_patches, directions = 'queen', cells = terra::cells(boundary_patches, i)$patches)
    if(length(intersect(rownames(cell_coord), rownames(neighbor_cells))) < nrow(cell_coord)) {
      rownames(neighbor_cells) <- as.numeric(rownames(neighbor_cells)) + 1
      neighor_cells = neighbor_cells + 1
    }
    
    edges <- data.frame(focal = character(), neighbor = character())
    for (j in 1:nrow(neighbor_cells)) {
      edges <- rownames(neighbor_cells)[j] %>%
        rep(8) %>%
        data.frame(focal = ., neighbor = neighbor_cells[j,]) %>%
        rbind(edges, .)
    }
    edges <- subset(edges, focal != neighbor & focal %in% rownames(neighbor_cells) & neighbor %in% rownames(neighbor_cells)) %>%
      as.matrix %>%
      trimws
    
    weights <- c()
    for (j in 1:nrow(edges)) {weights <- append(weights, dist.matrix[edges[j,1], edges[j,2]])}
    
    graph <- igraph::graph_from_edgelist(edges)
    igraph::edge.attributes(graph)$weight <- weights
    
    # and find the shortest distance along the subgraph between the furthest nodes
    lengths <- append(lengths, igraph::distances(graph, v = far_points[1], to = far_points[2]))
  }
  
  # identify the longest boundary
  longest_boundary <- max(lengths)

  # statistical test
  p <- null_distrib$longest_boundary(longest_boundary) %>%
    ifelse(. > 0.5, 1 - ., .) %>%
    as.numeric(.)

  names(p) <-'p-value'
  names(longest_boundary) <- 'longest_boundary'
  return(c(longest_boundary, p))
}
