#' test stan
#'
#' @param obs a vector of integers 
#'
#' @return the stan fit
#' @export
#'
#' @examples
a=function(obs) {
  optimize(stanmodels$nb, data=list(N=length(obs),observed=obs))
}
