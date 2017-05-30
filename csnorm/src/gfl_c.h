#ifndef GFL_C_H
#define GFL_C_H

#ifdef __cplusplus
extern "C" {
#endif
  
// This call is in C because it was meant to be called from within GFL, 
// but in the end it was not needed anymore. Could be ported to cpp
void cts_to_signal_mat_core(int N, int* cts_bin1, int* cts_bin2, double* count,
        double* lmu_nosig, double* weight, int nbins, double dispersion,
        double* phi, double eCprime, double* phihat, double* phihat_var, double* ncounts,
        int* bin1, int* bin2, int diag_rm);

#ifdef __cplusplus
}
#endif
    
#endif

