functions {
  //find the index in x such as x[i]<=z and x[i+1]>z
  //startpos is a suggested starting position for i
  real bisect(real z, int startpos, vector x) {
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
  int N;
  vector[N] values; //sorted
  real query;
  int<lower=1,upper=N> initial;
}
parameters {}
model {}
generated quantities {
  real x;
  x <- bisect(query, initial, values);
}
