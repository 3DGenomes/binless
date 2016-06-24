#' test stan
#'
#' @useDynLib test, .registration = TRUE
#' @param obs a vector of integers 
#'
#' @return the stan fit
#' @export
#'
#' @examples
a=function(obs) {
  optimizing(stanmodels$nb, data=list(N=length(obs),observed=obs))
}
