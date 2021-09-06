#' proportion mapping function
#'
#' this is changing the score to percentage ranking
#'
#'
#' @param x this is the input score vector
#' @param from a numerical value, this is the lower bound
#' @param to a numerical value, this is the upper bound
#' @return another vector with percentage ranking score
#'
#' @examples
#' data("rawdata", package = "simKAP")
#'
#'
#' new_score <- linMap(rawdata$recip_epts,1,100);
#' @importFrom stats na.omit
#' @export


linMap <- function(x, from, to) {
    (x - min(stats::na.omit(x))) / max(stats::na.omit(x) - min(stats::na.omit(x))) * (to - from) + from
}
