#' test stan
#'
#' @useDynLib test, .registration = TRUE
#' @param obs a vector of integers 
#'
#' @return the stan fit
#' @export
#'
#' @examples
b=function(obs) {
  optimizing(stanmodels$gauss, data=list(N=length(obs),observed=obs))
}
