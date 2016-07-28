  //design matrix and projector for dense SCAM cubic spline
  //needs declaration of the following
  //int<lower=1> N
  //real<lower=0> dmin
  //real<lower=dmin> dmax
  //vector<lower=dmin,upper=dmax>[N] dist;
  //matrix[N,Kdiag] Xdiag;
  //row_vector[Kdiag] pdiag;
  //row_vector[N] diag_weights;

  //diagonal SCAM spline, dense, exact and mean field model
  {
    diag_weights <- diag_weights/mean(diag_weights);
    Xdiag <- bspline(log(dist), Kdiag, splinedegree(), log(dmin), log(dmax));
    //projector for diagonal (SCAM)
    pdiag <- diag_weights * Xdiag;
    pdiag <- pdiag / (pdiag * rep_vector(1,Kdiag));
  }
