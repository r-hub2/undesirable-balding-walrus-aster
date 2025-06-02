#ifndef ASTER_MATOPS2_H
#define ASTER_MATOPS2_H

void affine(double *a, double *b, double *c, int *nrow, int *ncol);

void matvecmult(double *a, double *b, int *nrow, int *ncol, double *result);

void matmatmult(double *a, double *b, int *nrowa, int *ncola, int *ncolb,
    double *c);

void vecadd(double *a, double *b, int len);

void veccopy(double *a, double *b, int len);

double log_det_mat_plus_identity(double *a, int n);

double log_det_mat_plus_identity_fussy(double *a, int n);

double log_det_mat_plus_identity_deriv(double *zwz, int n,
    double *sigma, int *nrand, int len_nu, int k);

void mat_plus_identity_inverse(double *a, int n, double *result);

double log_det_mat_times_diag_plus_identity_obsolete(double *zwz, double *dee, int *nin);

void grad_log_det_mat_times_diag_plus_identity(double *zwz, double *nu,
   int *dee_idx, int *n_in, int *len_nu_in, double *result);

void LU_mat_times_diag_plus_identity(double *zwz, double *dee, int *nin,
    double *lu, int *ipiv);

double log_det_mat_times_diag_plus_identity(double *lu, int *nin);

void inverse_mat_times_diag_plus_identity(double *lu, int *ipiv, int *nin);

double partial_log_det_mat_times_diag_plus_identity(double *lu,
    int *dee_idx, double *zwz, int *n_in, int *k_in);

double second_partial_log_det_mat_times_diag_plus_identity(double *lu,
    int *dee_idx, double *zwz, int *n_in, int *k1_in, int *k2_in);

#include <R.h>
#include <Rinternals.h>

SEXP pos_def_mat_inv(SEXP a);

#endif /* ASTER_MATOPS2_H */
