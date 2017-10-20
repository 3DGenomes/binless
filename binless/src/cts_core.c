#include <math.h>
//#include <stdio.h>
#include "cts_core.h"

void cts_to_signal_mat_core(int N, int* cts_bin1, int* cts_bin2, double* count,
        double* lmu_nosig, double* weight, double* log_decay, int nbins, double dispersion,
        double* phi, double eCprime, double* phihat, double* phihat_var, double* phihat_var_nodecay,
        double* ncounts, int* bin1, int* bin2)
{
  //walk through cts
  int nbetas = nbins*(nbins+1)/2; //size of fused lasso problem
  int i;
  for (i = 0; i < N; ++i) {
    //pos(b1,b2) starts at 0 and is the index in the 2D triangle grid
    //of the cell at coordinates (b1,b2), where b1 and b2 start at 1 and b2>=b1
    int b1 = cts_bin1[i];
    int b2 = cts_bin2[i];
    int pos = (b1-1)*(nbins+1) - (b1-1)*b1/2 + b2 - b1;
    double mu = exp( lmu_nosig[i] + phi[pos] + eCprime);
    double mu_nodecay = exp( lmu_nosig[i] - log_decay[i] + phi[pos] + eCprime);
    double z = count[i]/mu-1;
    double var = 1/mu + 1/dispersion;
    double var_nodecay = 1/mu_nodecay + 1/dispersion;
    double w2v = weight[i]/(2*var);
    double w2v_nodecay = weight[i]/(2*var_nodecay);
    ncounts[pos] += weight[i];
    phihat_var[pos] += w2v;
    phihat_var_nodecay[pos] += w2v_nodecay;
    phihat[pos] += (z+phi[pos])*w2v;
  }
  
  //finish mat
  int b1 = 1;
  int b2 = 1;
  for (i = 0; i < nbetas; ++i) {
    phihat_var[i] = 1/phihat_var[i];
    phihat_var_nodecay[i] = 1/phihat_var_nodecay[i];
    if (ncounts[i]>0) {
      phihat[i] = phihat[i]*phihat_var[i];
    } else {
      phihat[i] = 0;
    }
    bin1[i] = b1;
    bin2[i] = b2++;
    if (b2 > nbins) b2 = ++b1;
  }
}

