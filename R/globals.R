#' @importFrom stats na.omit median
#' @importFrom utils setTxtProgressBar txtProgressBar head capture.output
#' @importFrom methods is
#' @importFrom magrittr %>%
#' @import igraph
NULL

#' @importFrom utils globalVariables
utils::globalVariables(c('.', 'long', 'lat', 'region', 'patches', 'p', 'focal', 'neighbor', 'values'))
