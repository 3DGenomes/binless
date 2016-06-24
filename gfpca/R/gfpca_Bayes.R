#' gfpca_Bayes
#' 
#' Implements a Bayesian approach to generalized functional principal components analysis for 
#' sparsely observed binary curves
#'
#' @useDynLib gfpca, .registration = TRUE  
#' @param data
#' @param npc
#' @param pve
#' @param grid
#' @param type
#' @param basis
#' @param nbasis
#' 
#' @references
#' Gertheiss, J., Goldsmith, J., and Staicu, A.-M. (2016).
#' A note on modeling sparse exponential-family functional response curves. 
#' \emph{Under Review}.
#' 
#' @author Jan Gertheiss \email{jan.gertheiss@@agr.uni-goettingen.de}
#' 
#' @seealso \code{\link{gfpca_TwoStep}}, \code{\link{gfpca_Bayes}}.
#'   
#' @import mgcv
#' @import gamm4
#' @import refund 
#' @export
gfpcaBayes <- function(data, npc=3, grid = NULL, nbasis=10, iter=1000, warmup=400){
  
  ## implement some data checks
  
  if(is.null(grid)){ grid = sort(unique(data['.index'][[1]])) }
  
  I = length(unique(data['.id'][[1]]))
  D = length(grid)
  
  X.des = matrix(1, nrow = I, ncol = 1)
  p = 1
  
  BS.pen = bs(grid, df=nbasis, intercept=TRUE, degree=3)
  BS = bs(data['.index'][[1]], df=nbasis, intercept=TRUE, degree=3)
  
  alpha = .1
  diff0 = diag(1, D, D)
  diff2 = matrix(rep(c(1,-2,1, rep(0, D-2)), D-2)[1:((D-2)*D)], D-2, D, byrow = TRUE)
  P0 = t(BS.pen) %*% t(diff0) %*% diff0 %*% BS.pen
  P2 = t(BS.pen) %*% t(diff2) %*% diff2 %*% BS.pen
  P.mat = alpha * P0 + (1-alpha) * P2
  
  
  ## format sparse data for stan fitting
  Y.vec.obs = data['.value'][[1]]
  subject.obs = data['.id'][[1]]
  n.total = length(Y.vec.obs)
  
  ## fit model using STAN
  dat = list(Y = Y.vec.obs, X = X.des, BS = BS,
             subjId = subject.obs,
             N = n.total, I = I, D = D, p = p, Kt = nbasis, Kp = npc, 
             PenMat = P.mat)
  stanfit <- stanmodels$gfpca
  GenFPCA.fit = sampling(stanfit,
                         data=dat, iter = iter, warmup = warmup,
                         control = list(adapt_delta = .65), 
                         chains = 1, verbose = FALSE)
  
  ## post process to obtain point estimates
  beta.post = extract(GenFPCA.fit, "beta")$beta
  beta_psi.post = extract(GenFPCA.fit, "beta_psi")$beta_psi
  c.post = extract(GenFPCA.fit, "c")$c
  
  betaHat.post = array(NA, dim = c(p, D, dim(c.post)[1]))
  for(i in 1:dim(c.post)[1]) {
    betaHat.post[,,i] = (beta.post[i,,] %*% t(BS.pen))
  }
  
  y.post = z.post = array(NA, dim = c(I, D, dim(c.post)[1]))
  for(i in 1:dim(c.post)[1]) {
    y.post[,,i] = X.des %*% (beta.post[i,,] %*% t(BS.pen)) + c.post[i,,] %*% (beta_psi.post[i,,] %*% t(BS.pen))
    z.post[,,i] = c.post[i,,] %*% (beta_psi.post[i,,] %*% t(BS.pen))
  }
  Zstan = apply(z.post, c(1,2), mean)
  W.bayes = apply(y.post, c(1,2), mean)
  
  ret = list(apply(betaHat.post, 2, mean), Zstan, W.bayes)
  names(ret) = c("mu", "z", "yhat")
  ret

}
  
