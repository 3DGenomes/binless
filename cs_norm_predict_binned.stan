////
// Cut-site normalization model: predict mean for binned matrices
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
  
  int splinedegree() {return 3;} //set 3 for cubic spline
  
  matrix bspline(vector x, int K, int q, real xmin, real xmax) {
    real dx; //interval length
    row_vector[K] t; //knot locations (except last)
    //
    dx <- 1.01*(xmax-xmin)/(K-q); //make it slightly larger
    t <- xmin - dx*0.01 + dx * range(-q,K-q-1)';
    return bspl_gen(x, dx, t, q);
  }
  
  //find the index in x such as x[i]<=z and x[i+1]>z
  //startpos is a suggested starting position for i
  int bisect(real z, int startpos, vector x) {
    int pos;
    int left;
    int right;
    int nmax;
    int ctr;
    left <- 1;
    right <- rows(x);
    ctr <- 1;
    nmax <- rows(x);
    pos <- startpos;
    //return edges if queried
    if (z<x[1]) return 1;
    if (z>x[nmax]) return nmax;
    //shortcut if z has same pos or the one after
    if (pos < right && x[pos] <= z && x[pos+1] > z) return pos;
    if (pos < right-1 && x[pos+1] <= z && x[pos+2] > z) return pos;
    //bisection algorithm
    while (ctr <= nmax && left != right) {
      ctr <- ctr + 1;
      print(ctr, " left=", left, " (",x[left], ") right=", right, " (", x[right], ") pos=", pos, " (", x[pos], ")");
      if ( (x[right]-z)*(x[pos]-z) > 0 ) {right <- pos;} else {left <- pos;}
      if (right-left==1 && x[right] > z && x[left] <= z) return left;
      pos <- (right+left)/2;
    }
    return pos;
  }
}
////////////////////////////////////////////////////////////////
data {
  //spline parameters
  int Kdiag; //number of functions in spline base for diagonal decay
  int npoints; //number of points to compute
  int circularize; //if >0 assume genome is circular and has this size
  //biases
  int<lower=1> S1; //number of cut sites
  int<lower=1> S2;
  vector<lower=0>[S1] cutsites1; //cut site locations, better be sorted
  vector<lower=0>[S2] cutsites2;
  //distance bounds
  real<lower=0> dmin;
  real<lower=dmin> dmax;
  //counts
  int<lower=0> N; //number of nonzero counts to bin
  int<lower=1> counts[N]; //value of the count
  vector<lower=0>[N] cdist; //distance of that count  (used for normalized matrix only)
  vector<lower=0>[N] cmean; //posterior mean of the count  (used for normalized matrix only)
  //estimated parameters
  real eC;  //exposure for counts
  vector[S1] log_nu1; //take nu and delta directly to avoid base reconstruction
  vector[S2] log_nu2;
  vector[S1] log_delta1; 
  vector[S2] log_delta2; 
  vector[Kdiag] beta_diag_centered; //need to build spline base
  //binned matrix info
  int<lower=1> B1; //binned matrix dimensions
  int<lower=1> B2;
  int<lower=1,upper=B1> bbins1[S1]; //in which bin the biases fall
  int<lower=1,upper=B2> bbins2[S2];
  int<lower=1,upper=B1> cbins1[N]; //in which bin the counts fall
  int<lower=1,upper=B2> cbins2[N];
}
transformed data {
  if (circularize>0) {
    if (max(cutsites2)-min(cutsites1) > circularize) {
      reject("circular genome size smaller than maximum distance between cutsites!");
    }
  }
}
parameters {}
model {}
generated quantities {
  //decay
  vector[npoints] log_dist;
  vector[npoints] log_decay;
  //matrices
  matrix[B1,B2] ncounts; //number of possible counts in bin
  matrix[B1,B2] observed; //summed counts per bin
  matrix[B1,B2] expected; //posterior mean of negative binomial per bin
  matrix[B1,B2] normalized; // (sum_i observed_i * decay_i / expected_i) * ncounts
  
  //decay
  {
    matrix[npoints,Kdiag] Xdiag;
    log_dist <- range(0,npoints-1)*(log(dmax)-log(dmin))/(npoints-1)+log(dmin);
    Xdiag <- bspline(log_dist, Kdiag, splinedegree(), log(dmin), log(dmax));
    log_decay <- Xdiag * beta_diag_centered;
  }

  //observed and normalized
  observed <- rep_matrix(0,B1,B2);
  normalized <- rep_matrix(0,B1,B2);
  for (i in 1:N) { //do not vectorize to avoid aliasing issues
    int b1;
    int b2;
    int k;
    b1<-cbins1[i];
    b2<-cbins2[i];
    observed[b1,b2] <- observed[b1,b2] + counts[i];
    k <- bisect(log(cdist[i]), npoints/2, log_dist);
    normalized[b1,b2] <- normalized[b1,b2] + counts[i]/cmean[i]*exp(log_decay[k]);
  }
  
  //expected and ncounts
  ncounts <- rep_matrix(0,B1,B2);
  expected <- rep_matrix(0,B1,B2);
  for (i in 1:S1) {
    real pos1;
    int k;
    pos1 <- cutsites1[i];
    k <- 1;
    for (j in 1:S2) {
      real pos2;
      pos2 <- cutsites2[j];
      if (pos2 > pos1) {
        int b1;
        int b2;
        b1<-bbins1[i];
        b2<-bbins2[j];
        ncounts[b1,b2] <- ncounts[b1,b2]+4; //1 for each count type
        if (circularize>0) {
          if (pos2-pos1>circularize/2) {
            k <- bisect(log(circularize+1-(pos2-pos1)), k, log_dist);
          } else {
            k <- bisect(log(pos2-pos1), k, log_dist);
          }
        } else {
          k <- bisect(log(pos2-pos1), k, log_dist);
        }
        expected[b1,b2] <- expected[b1,b2] +
            exp(eC + log_nu1[i] + log_nu2[j] + log_decay[k])*2*cosh(log_delta1[i])*2*cosh(log_delta2[j]);
      }
    }
  }
  normalized <- normalized ./ ncounts;
}