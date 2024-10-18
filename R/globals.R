#' @importFrom stats na.omit rnorm median
#' @importFrom utils setTxtProgressBar txtProgressBar head capture.output
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @import igraph
NULL

#' @importFrom utils globalVariables
utils::globalVariables(c('.', 'long', 'lat', 'region', 'patches', 'nb', 'wt', 'p_ii', 'geometry', 'values', 'local_moran', 'cluster'))
