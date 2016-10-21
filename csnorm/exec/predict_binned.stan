////
// Cut-site normalization model: predict mean for binned matrices
////
functions {
  #include "common_functions.stan"
}
////////////////////////////////////////////////////////////////
data {
  //biases
  int<lower=1> S1; //number of cut sites
  int<lower=1> S2;
  vector<lower=0>[S1] cutsites1; //cut site locations, better be sorted
  vector<lower=0>[S2] cutsites2;
  //counts
  int<lower=0> N; //number of nonzero counts to bin
  int cidx[2,N]; //indices of rsite pairs, pointing to within cutsites1 and cutsites2
  int<lower=0> counts_close[N]; //value of the count
  int<lower=0> counts_far[N];
  int<lower=0> counts_up[N];
  int<lower=0> counts_down[N];
  //binned matrix info
  int<lower=1> B1; //binned (sub)matrix dimensions
  int<lower=1> B2;
  int<lower=1,upper=B1> bbins1[S1]; //in which bin the biases fall
  int<lower=1,upper=B2> bbins2[S2];
  int<lower=1,upper=B1> cbins1[N]; //in which bin the counts fall
  int<lower=1,upper=B2> cbins2[N];
  //spline parameters
  int Kdiag; //number of functions in spline base for diagonal decay
  int npoints; //number of points to compute the decay on
  int circularize; //if >0 assume genome is circular and has this size
  real<lower=0> dmin;
  real<lower=dmin> dmax;
  //estimated parameters
  real eC;  //exposure for counts
  vector[S1] log_iota1; //take iota and rho directly to avoid base reconstruction
  vector[S2] log_iota2;
  vector[S1] log_rho1; 
  vector[S2] log_rho2; 
  vector[Kdiag] beta_diag_centered; //need to build spline base for diagonal decay
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
  matrix[B1,B2] decaymat; // average decay per bin
  
  //decay
  {
    matrix[npoints,Kdiag] Xdiag;
    log_dist = range(0,npoints-1)*(log(dmax)-log(dmin))/(npoints-1)+log(dmin);
    Xdiag = bspline(log_dist, Kdiag, splinedegree(), log(dmin), log(dmax));
    log_decay = Xdiag * beta_diag_centered;
  }

  //observed
  observed = rep_matrix(0,B1,B2);
  for (i in 1:N) { //do not vectorize to avoid aliasing issues
    int b1;
    int b2;
    b1=cbins1[i];
    b2=cbins2[i];
    observed[b1,b2] = observed[b1,b2] + counts_close[i] + counts_far[i] + counts_up[i] + counts_down[i];
  }
  
  //expected, ncounts and decaymat
  ncounts = rep_matrix(0,B1,B2);
  expected = rep_matrix(0,B1,B2);
  decaymat = rep_matrix(0,B1,B2);
  for (i in 1:S1) {
    real pos1;
    int k;
    pos1 = cutsites1[i];
    k = 1;
    for (j in 1:S2) {
      real pos2;
      real dist;
      pos2 = cutsites2[j];
      dist = pos2-pos1;
      if (dist >= dmin) {
        int b1;
        int b2;
        b1=bbins1[i];
        b2=bbins2[j];
        ncounts[b1,b2] = ncounts[b1,b2]+4; //1 for each count type
        if (circularize>0) {
          if (dist>circularize/2) {
            k = bisect(log(circularize+1-dist), k, log_dist);
          } else {
            k = bisect(log(dist), k, log_dist);
          }
        } else {
          k = bisect(log(dist), k, log_dist);
        }
        expected[b1,b2] = expected[b1,b2] +
            exp(eC + log_decay[k])*(exp(log_iota1[i]) + exp(log_rho[i]))*(exp(log_iota1[j]) + exp(log_rho[j]));
        decaymat[b1,b2] = decaymat[b1,b2] + exp(4*log_decay[k]);
      }
    }
  }
  decaymat = decaymat ./ ncounts;
}
