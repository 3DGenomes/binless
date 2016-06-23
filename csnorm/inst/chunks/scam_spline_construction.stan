  //design matrix and projector for dense SCAM cubic spline
  //needs declaration of the following
  //int<lower=1> N
  //real<lower=0> dmin
  //real<lower=dmin> dmax
  //vector<lower=dmin,upper=dmax>[N] dist;
  //matrix[N,Kdiag] Xdiag;
  //row_vector[Kdiag] pdiag;

  //diagonal SCAM spline, dense, exact and mean field model
  {
    Xdiag <- bspline(log(dist), Kdiag, splinedegree(), log(dmin), log(dmax));
    //projector for diagonal (SCAM)
    pdiag <- rep_row_vector(1,N) * Xdiag;
    pdiag <- pdiag / (pdiag * rep_vector(1,Kdiag));
  }
