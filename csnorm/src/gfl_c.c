#include <math.h>
//#include <stdio.h>
#include "gfl_c.h"

void cts_to_signal_mat_core(int N, int* cts_bin1, int* cts_bin2, double* count,
        double* lmu_nosig, double* weight, int nbins, double dispersion,
        double* phi, double eCprime, double* phihat, double* phihat_var, double* ncounts,
        int* bin1, int* bin2, int diag_rm)
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
    double z = count[i]/mu-1;
    double var = 1/mu + 1/dispersion;
    double w2v = weight[i]/(2*var);
    ncounts[pos] += weight[i];
    phihat_var[pos] += w2v;
    phihat[pos] += (z+phi[pos])*w2v;
    /*if (pos==108) printf("count=%f lmu_nosig=%f phi=%f mu=%f z=%f var=%f w2v=%f pos=%d\n",
        count[i], lmu_nosig[i], phi[pos], mu,z,var,w2v,pos);*/
  }
  
  //finish mat
  int b1 = 1;
  int b2 = 1;
  for (i = 0; i < nbetas; ++i) {
    phihat_var[i] = 1/phihat_var[i];
    if (ncounts[i]>0) {
      phihat[i] = phihat[i]*phihat_var[i];
    } else {
      phihat[i] = 0;
    }
    bin1[i] = b1;
    bin2[i] = b2++;
    if (b2 > nbins) b2 = ++b1;
    if (bin2[i]-bin1[i] <= diag_rm) phihat_var[i] = INFINITY;
    //if (bin1[i]==31 & bin2[i]==31) printf("phihat=%f phihat.var=%f i=%d\n",phihat[i],phihat_var[i],i);
  }
}

