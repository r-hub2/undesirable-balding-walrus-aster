// See section 6.6.1 of Writing R Extensions (Fortran character strings)
#define USE_FC_LEN_T

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include "matops2.h"

#ifndef FCONE
#define FCONE
#endif

// A^T b + c and leave result in c
void affine(double *a, double *b, double *c, int *nrow, int *ncol)
{
    double one = 1.0;
    int ione = 1;
    F77_CALL(dgemv)("n", nrow, ncol, &one, a, nrow, b, &ione, &one, c,
        &ione FCONE);
}

// A^T b
void matvecmult(double *a, double *b, int *nrow, int *ncol, double *result)
{
    double zero = 0.0;
    double one = 1.0;
    int ione = 1;
    F77_CALL(dgemv)("n", nrow, ncol, &one, a, nrow, b, &ione, &zero, result,
        &ione FCONE);
}

// A B
void matmatmult(double *a, double *b, int *nrowa, int *ncola, int *ncolb,
    double *c)
{
    double one = 1.0;
    double zero = 0.0;
    F77_CALL(dgemm)("n", "n", nrowa, ncolb, ncola, &one, a, nrowa, b, ncola,
        &zero, c, nrowa FCONE FCONE);
}

// a + b
void vecadd(double *a, double *b, int len)
{
    double one = 1.0;
    int ione = 1;
    F77_CALL(daxpy)(&len, &one, a, &ione, b, &ione);
}

void veccopy(double *a, double *b, int len)
{
    int ione = 1;
    F77_CALL(dcopy)(&len, a, &ione, b, &ione);
}

// log(det(A + I))
// A is real symmetric
// use fact that eigenvalues of A + I are 1 + eigenvalues of A 
// use log1p to deal with log(1 + . )
double log_det_mat_plus_identity_fussy(double *a, int n)
{
    // Check A for Nan or Inf
    for (int i = 0, nsq = n * n; i < nsq; i++)
        if (! R_finite(a[i]))
           return R_PosInf; // Probably the Wrong Thing, but what better?
    // Ask for optimal size of work arrays
    double foo = 0.0; // Not actually referenced
    int ifoo = 0; // Not actually referenced
    double abstol = -1; // means figure out a reasonable tolerance by itself
    int nfound;
    double *values = (double *) R_alloc(n, sizeof(double));
    double z = 0.0; // Not actually referenced
    int lda_z = 1; // Not actually referenced
    int *isuppz = (int *) R_alloc(2 * n, sizeof(int));
    int lwork = -1, liwork = -1; // Will be set by first call to dsyerv
    double tmp; int itmp; // Not actually referenced
    int info;
    F77_CALL(dsyevr)("N", "A", "L", &n, a, &n, &foo, &foo, &ifoo, &ifoo,
        &abstol, &nfound, values, &z, &lda_z, isuppz, &tmp, &lwork, &itmp,
        &liwork, &info FCONE FCONE FCONE);
    if (info != 0)
        error("error code %d from Lapack routine dsyevr", info);
    lwork = (int) tmp;
    liwork = itmp;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    int *iwork = (int *) R_alloc(liwork, sizeof(int));
    F77_CALL(dsyevr)("N", "A", "L", &n, a, &n, &foo, &foo, &ifoo, &ifoo,
        &abstol, &nfound, values, &z, &lda_z, isuppz, work, &lwork, iwork,
        &liwork, &info FCONE FCONE FCONE);
    if (info != 0)
        error("error code %d from Lapack routine dsyevr", info);
    // REVISED DOWN TO HERE
    double result = 0.0;
    for (int i = 0; i < nfound; i++) {
        double foo = values[i];
        if (foo > 0.0)
            result += log1p(foo);
    }
    return result;
}

// log(det(A + I))
// A is real symmetric
// use Cholesky of A + I, avoid eigenvalues of A as too slow
double log_det_mat_plus_identity_obsolete(double *a_in, int n)
{
    // Be nondestructive
    double *a = (double *) R_alloc(n * n, sizeof(double));
    memcpy(a, a_in, n * n * sizeof(double));

    // Check A for Nan or Inf
    for (int i = 0, nsq = n * n; i < nsq; i++)
        if (! R_finite(a[i]))
           return R_PosInf; // Probably the Wrong Thing, but what better?
    // Add one to diagonal of A
    for (int i = 0; i < n; i++)
        a[i * (n + 1)] += 1.0;

    int *piv = R_Calloc(n, int);
    int rank;
    double tol = -1.0; // figure it out yourself
    double *work = R_Calloc(2 * n, double);
    int info;
    F77_CALL(dpstrf)("L", &n, a, &n, piv, &rank, &tol, work, &info FCONE);
    if (info < 0)
        error("argument %d of Lapack routine dsptrf had invalid value", -info);
    if (info == 0)
        rank = n;

    // result is sum of squares of diagonal of A (which is L after dpstrf);
    double result = 0.0;
    for (int i = 0; i < rank; i++)
        result += 2.0 * log(a[i * (n + 1)]);

    R_Free(piv);
    R_Free(work);
    return result;
}

#include <Rinternals.h>

// A^{-1}
SEXP pos_def_mat_inv(SEXP a)
{
    if (! isMatrix(a))
        error("argument not matrix");
    if (! isReal(a))
        error("argument not storage mode double");
    SEXP result;
    PROTECT(result = duplicate(a));
    int n = nrows(a);
    if (n != ncols(a))
        error("argument not square matrix");

    int info;
    F77_CALL(dpotrf)("L", &n, REAL(result), &n, &info FCONE);
    if (info != 0)
        error("Cholesky decomposition failed");
    F77_CALL(dpotri)("L", &n, REAL(result), &n, &info FCONE);
    if (info != 0)
        error("inverse computation failed");
    // but we still only have the lower triangle of result
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
            REAL(result)[j + n * i] = REAL(result)[i + n * j];
    UNPROTECT(1);
    return result;
}

/* arguments
 *       a: (input) square positive definite matrix
 *       n: (input) integer, zwz is n by n
 *  result: (output) [ w + Id ]^{- 1}
 *
 */

void mat_plus_identity_inverse(double *a, int n, double *result)
{
    memcpy(result, a, n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        result[i + n * i] += 1.0;

    int info;
    F77_CALL(dpotrf)("L", &n, result, &n, &info FCONE);
    if (info != 0)
        error("Cholesky decomposition failed");
    F77_CALL(dpotri)("L", &n, result, &n, &info FCONE);
    if (info != 0)
        error("inverse computation failed");
    // but we still only have the lower triangle of result
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
            result[j + n * i] = result[i + n * j];
}

/* arguments
 *     zwz: square positive definite matrix
 *     cee: vector of same dimension, standardized random effects
 *   nrand: vector of length number of variance components,
 *              which gives number of random effects for each variance component
 *   sigma: vector of square roots of variance components
 *       k: integer such that sigma[k] makes sense in C (0-origin indexing)
 *
 * returns: tr( [ A zwz A + Id ]^{- 1} [ E_k zwz A + A zwz E_k \bigr] )
 *              where A is a diagonal matrix whose diagonal is
 *              rep(sigma, times = nrand) and E_k is partial A / partial sigma_k
 *
 */

double log_det_mat_plus_identity_deriv(double *zwz, int n,
    double *sigma, int *nrand, int len_nu, int k)
{
    if (k < 0 || k >= len_nu)
        error("argument k out of range");
    for (int i = 0; i < len_nu; i++)
        if (nrand[i] <= 0)
            error("argument nrand must be positive-integer-valued");

    double *eek = R_Calloc(n, double);
    for (int i = 0, m = 0; i < len_nu; i++)
        for (int j = 0; j < nrand[i]; j++, m++)
            eek[m] = (i == k) ? 1.0 : 0.0;

    double *aaa = R_Calloc(n, double);
    for (int i = 0, m = 0; i < len_nu; i++)
        for (int j = 0; j < nrand[i]; j++, m++)
            aaa[m] = sigma[i];

    double *work = R_Calloc(n * n, double);
    memcpy(work, zwz, n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            work[i + n * j] *= aaa[i] * aaa[j];
    for (int i = 0; i < n; i++)
            work[i + n * i] += 1.0;

    int info;
    F77_CALL(dpotrf)("L", &n, work, &n, &info FCONE);
    if (info != 0)
        error("Cholesky decomposition failed");
    F77_CALL(dpotri)("L", &n, work, &n, &info FCONE);
    if (info != 0)
        error("inverse computation failed");
    // but we still only have the lower triangle of result
    for (int i = 0; i < n; i++)
        for (int j = 0; j < i; j++)
            work[j + n * i] = work[i + n * j];
    // now work is the inverse part of the result

    double *work_too = R_Calloc(n * n, double);
    memcpy(work_too, zwz, n * n * sizeof(double));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            work_too[i + n * j] *= aaa[i] * eek[j];

    double result = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result += work[i + n * j] *
                (work_too[j + n * i] + work_too[i + n * j]);

    R_Free(eek);
    R_Free(aaa);
    R_Free(work);
    R_Free(work_too);
    return result;
}

// calculates log(det(zwz %*% diag(dee) + diag(n)))
// where n = nin[0]
// zwz is n by n and symmetric positive semidefinite (latter not verified)
// dee is length(n)
double log_det_mat_times_diag_plus_identity_defunct(double *zwz, double *dee, int *nin)
{
    int n = nin[0];
    double *foo = (double *) R_alloc(n * n, sizeof(double));
    memcpy(foo, zwz, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            foo[i + n * j] *= dee[j];
        foo[i + n * i] += 1.0;
    }
    int *ipiv = (int *) R_alloc(n, sizeof(int));

    int info;
    F77_CALL(dgetrf)(&n, &n, foo, &n, ipiv, &info);
    if (info < 0)
        error("argument %d of LAPACK dgetrf had illegal value", -info);
    // if (info > 0)
    //     error("LU decomposition uninvertable");

    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += log(foo[i + n * i]);
    return result;
}

// calculates derivative of log(det(zwz %*% diag(dee) + diag(n))) with
//     respect to nu, where
// n = nin[0]
// zwz is n by n and symmetric positive semidefinite (latter not verified)
// dee = nu[dee.idx]
// dee_idx is integer vector of length n having values in 1:length(nu) - 1
//     zero-origin indexing for nu, like C not R
// result is length len_nu
void grad_log_det_mat_times_diag_plus_identity(double *zwz, double *nu,
   int *dee_idx, int *n_in, int *len_nu_in, double *result)
{
    int n = n_in[0];
    int len_nu = len_nu_in[0];

    double *foo = (double *) R_alloc(n * n, sizeof(double));
    memcpy(foo, zwz, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            foo[i + n * j] *= nu[dee_idx[j]];
        foo[i + n * i] += 1.0;
    }
    int *ipiv = (int *) R_alloc(n, sizeof(int));

    int info;
    F77_CALL(dgetrf)(&n, &n, foo, &n, ipiv, &info);
    if (info < 0)
        error("argument %d of LAPACK dgetrf had illegal value", -info);
    if (info > 0)
        error("LU decomposition uninvertable");

    int lwork = -1;
    double work_size;
    F77_CALL(dgetri)(&n, foo, &n, ipiv, &work_size, &lwork, &info);
    if (info != 0)
        error("LAPACK dgetri failed to return optimal workspace size");
    lwork = work_size;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgetri)(&n, foo, &n, ipiv, work, &lwork, &info);
    if (info < 0)
        error("argument %d of LAPACK dgetri had illegal value", -info);
    if (info > 0)
        error("LU decomposition uninvertable");
    // now foo is inverse of (zwz %*% D + identity)

    for (int i = 0; i < len_nu; i++)
        result[i] = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            result[dee_idx[i]] += foo[i + n * j] * zwz[j + n * i];
}

// complete redo of above, breaking into pieces
// step 1: LU decomposition of zwz %*% dee + identity
// step 2: log determinant from LU decomposition
// step 3: inverse of zwz %*% dee + identity from LU decomposition
// step 4: gradient of log determinant from inverse
// step 5: hessian of log determinant from inverse
//     point is that we only do the LU decomposition once

// step 1: LU decomposition of zwz %*% dee + identity
// input arguments
// n is nin[0]
// zwz is n by n and symmetric positive semidefinite (latter not verified)
// dee is length(n)
// output arguments
// lu is n by n matrix contains L and U
// ipiv is vector of length n, pivoting information
void LU_mat_times_diag_plus_identity(double *zwz, double *dee, int *nin,
    double *lu, int *ipiv)
{
    int n = nin[0];
    memcpy(lu, zwz, n * n * sizeof(double));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            lu[i + n * j] *= dee[j];
        lu[i + n * i] += 1.0;
    }

    int info;
    F77_CALL(dgetrf)(&n, &n, lu, &n, ipiv, &info);
    if (info < 0)
        error("argument %d of LAPACK dgetrf had illegal value", -info);
}

// step 2: log determinant from LU decomposition
// input arguments
// n is nin[0]
// lu is n by n matrix produced in step 1
// result is returned scalar
double log_det_mat_times_diag_plus_identity(double *lu, int *nin)
{
    int n = nin[0];
    double result = 0.0;
    for (int i = 0; i < n; i++)
        result += log(lu[i + n * i]);
    return result;
}

// step 3: inverse of zwz %*% dee + identity from LU decomposition
// input arguments
// n is nin[0]
// lu is n by n matrix produced in step 1
// ipiv is n vector produced in step 2
// input arguments
// lu is updated destructively, contains L and U going in, contains inverse
//     going out
void inverse_mat_times_diag_plus_identity(double *lu, int *ipiv, int *nin)
{
    int n = nin[0];
    int info;
    int lwork = -1;
    double work_size;
    F77_CALL(dgetri)(&n, lu, &n, ipiv, &work_size, &lwork, &info);
    if (info != 0)
        error("LAPACK dgetri failed to return optimal workspace size");
    lwork = work_size;
    double *work = (double *) R_alloc(lwork, sizeof(double));
    F77_CALL(dgetri)(&n, lu, &n, ipiv, work, &lwork, &info);
    if (info < 0)
        error("argument %d of LAPACK dgetri had illegal value", -info);
    if (info > 0)
        error("LU decomposition uninvertable");
}

// step 4: gradient of log determinant from LU decomposition
// input arguments
// n is n_in[0]
// len_nu is len_nu_in[0]
// lu is n by n matrix produced in step 3: (zwz %*% dee + identity)^-1
// zwz is n by n matrix
// dee is nu[dee_idx]
// dee_idx is integer vector of length n having values in 1:length(nu) - 1
//     zero-origin indexing for nu, like C not R,
//     so in C dee is nu[dee_idx]
// k is k_in[0]
// result is partial with respect to nu[k]
double partial_log_det_mat_times_diag_plus_identity(double *lu,
    int *dee_idx, double *zwz, int *n_in, int *k_in)
{
    int k = k_in[0];
    int n = n_in[0];

    double result = 0.0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (dee_idx[i] == k)
                result += lu[i + n * j] * zwz[j + n * i];
    return result;
}

// step 5: hessian of log determinant from LU decomposition
// input arguments
// n is n_in[0]
// len_nu is len_nu_in[0]
// lu is n by n matrix produced in step 3: (zwz %*% dee + identity)^-1
// zwz is n by n matrix
// dee is nu[dee_idx]
// dee_idx is integer vector of length n having values in 1:length(nu) - 1
//     zero-origin indexing for nu, like C not R,
//     so in C dee is nu[dee_idx]
// k1 is k1_in[0]
// k2 is k2_in[0]
// result is second partial with respect to nu[k1] and nu[k2]
double second_partial_log_det_mat_times_diag_plus_identity(double *lu,
    int *dee_idx, double *zwz, int *n_in, int *k1_in, int *k2_in)
{
    int k1 = k1_in[0];
    int k2 = k2_in[0];
    int n = n_in[0];

    double result = 0.0;
    // note no indentation (would be 7 deep)
    for (int i1 = 0; i1 < n; i1++)
    if (dee_idx[i1] == k2)
    for (int i2 = 0; i2 < n; i2++)
    for (int i3 = 0; i3 < n; i3++)
    if (dee_idx[i3] == k1)
    for (int i4 = 0; i4 < n; i4++)
    result += lu[i1 + n * i2] * zwz[i2 + n * i3] *
              lu[i3 + n * i4] * zwz[i4 + n * i1];
    return result;
}
