#ifndef GFL_C_H
#define GFL_C_H

#ifdef __cplusplus
extern "C" {
#endif
  
void cts_to_signal_mat_core(int N, int* cts_bin1, int* cts_bin2, double* count,
        double* lmu_nosig, double* weight, int nbins, double dispersion,
        double* phi, double* phihat, double* phihat_var, double* ncounts,
        int* bin1, int* bin2, int diag_rm);

#ifdef __cplusplus
}
#endif
    
#endif

