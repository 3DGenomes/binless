////
// Spline base construction in various ways
////
functions {
  ///////////
  //BEGIN internal use
  //
  //generate [imin, imin+1, ... imax]
  vector range(int imin, int imax) {
    return cumulative_sum(rep_vector(1,imax-imin+1))-1+imin;
  }
  //Compute cumulative histogram of x at left cut points q
  //i.e. index of first x value that falls in bins defined by cutpoints defined in q
  //assumes x is ordered.
  //x[hist[i]] is the first value in x that is >= q[i]
  int[] cumulative_hist(vector x, row_vector q) {
    int indices[cols(q)];
    int ix;
    int N;
    ix <- 1;
    N <- rows(x);
    for (iq in 1:cols(q)) {
      while (ix < N && x[ix] < q[iq]) ix <- ix + 1;
      indices[iq] <- ix;
    }
    return indices;
  }
  int splinedegree() {return 3;} //set 3 for cubic spline
  int difforder() {return 2;} //set 2 for 2nd order difference penalty
  matrix bspl_gen(vector x, real dx, row_vector t, int q) {
    int N;
    int K;
    K <- cols(t);
    N <- rows(x);
    {
      int r[K];
      matrix[N, K] T;
      matrix[N, K] X;
      matrix[N, K] P;
      matrix[N, K] B;
      for (i in 2:K) r[i-1] <- i;
      r[K] <- 1;
      T <- rep_matrix(t, N);
      X <- rep_matrix(x, K);
      P <- (X - T) / dx;
      for (i in 1:N)
        for (j in 1:K)
          B[i,j] <- (T[i,j] <= X[i,j]) && (X[i,j] < T[i,j]+dx);
      for (k in 1:q)
        B <- ( P .* B + (k+1-P) .* B[,r]) / k;
      return B;
    }
  }
  //END internal use
  ///////////
  
  // Evaluate b-spline basis functions of degree q at the values
  // in x within support [xmin,xmax] and with K basis functions.
  // Use q=3 for cubic spline.
  // adapted from Eilers and Marx, Stat Sci, 1996
  matrix bspline(vector x, int K, int q) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    //
    dx <- 1.01*(max(x)-min(x))/(K-q); //make it slightly larger
    t <- min(x) - dx*0.01 + dx * range(-q,K-q-1)';
    return bspl_gen(x, dx, t, q);
  }

  matrix bspline_piecewise(vector x, int K, int q) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    int x_n[K-q+1]; //cumulative histogram of x values
    matrix[rows(x),K] Xs; //spline base
    print("bspline calculation with K=",K," and q=",q);
    //
    dx <- 1.01*(max(x)-min(x))/(K-q); //make it slightly larger
    t <- min(x) - dx*0.01 + dx * range(-q,K-q-1)';
    //get indices of x which cross to the next segment and build x_n.
    x_n[1] <- 1;
    x_n[2:(K-q)] <- cumulative_hist(x, t[(q+2):]);
    x_n[K-q+1] <- rows(x)+1;
    print("cut points: ",t);
    print("histogram: ",x_n);
    //
    //build spline, interval per interval
    Xs <- rep_matrix(0,rows(x),K);
    for (i in 1:(K-q)) {
      //at any x, there are q+1 nonzero b-splines. Compute them.
      int xbegin;
      int xend;
      xbegin <- x_n[i];
      xend <- x_n[i+1]-1;
      Xs[xbegin:xend,i:(i+q)] <- bspl_gen(x[xbegin:xend], dx, t[i:(i+q)], q);
      print("   loop iteration ",i, " x=", xbegin, ":",xend, " K=", i, ":", i+q);
    }
    return Xs;
  }
  
  //difference penalty for P-splines
  //K is the size of the parameter vector
  //d is the order of the difference
  //Returns P, which is K-d x K
  //D=t(P)P is the difference matrix, P*beta is the difference vector
  matrix difference_matrix_sqrt(int K, int d) {
    matrix[K,K] P;
    P <- diag_matrix(rep_vector(1, K));
    for (i in 1:d) {
      matrix[K-i,K] tmp; //stan 2.9.0 issues compile error for deep_copy
      tmp <- P[2:(K-i+1),]-P[1:(K-i),];
      P[1:(K-i),] <- tmp;
    }
    return P[1:(K-d),];
  }
}
////////////////////////////////////////////////////////////////
data {
  int<lower=1> N; //number of data points
  int<lower=1> K; //number of parameters
  vector[N] y; //dependent variable
  vector[N] x; //independent variable
}
parameters {
}
model {
}
generated quantities {
  matrix[N,K] Xw; //design matrix
  matrix[N,K] Xp; //design matrix
  Xw <- bspline(x,K,splinedegree());
  Xp <- bspline_piecewise(x,K,splinedegree());
}
