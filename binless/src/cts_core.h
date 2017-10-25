#ifndef CTS_CORE_H
#define CTS_CORE_H

#ifdef __cplusplus
extern "C" {
#endif
  
// This call is in C because it was meant to be called from within GFL, 
// but in the end it was not needed anymore. Could be ported to cpp
void cts_to_signal_mat_core(int N, int* cts_bin1, int* cts_bin2, double* count,
        double* lmu_nosig, double* weight, double* log_decay, int nbins, double dispersion,
        double* phi, double* phihat, double* phihat_var, double* phihat_var_nodecay,
        double* ncounts, int* bin1, int* bin2);

#ifdef __cplusplus
}
#endif
    
#endif

