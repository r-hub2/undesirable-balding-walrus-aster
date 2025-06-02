
#include "aster.h"

/* does a %*% b where a is matrix b is vector */

void aster_mat_vec_mult(int *nrowin, int *ncolin, double *a, double *b,
    double *c)
{
    int nrow = nrowin[0];
    int ncol = ncolin[0];

    int i, j, k;

    for (i = 0; i < nrow; ++i)
        c[i] = 0.0;
    for (j = 0, k = 0; j < ncol; ++j)
        for (i = 0; i < nrow; ++i, ++k)
            c[i] += a[k] * b[j];
}

/* does b %*% a where a is matrix b is vector */

void aster_vec_mat_mult(int *nrowin, int *ncolin, double *a, double *b,
    double *c)
{
    int nrow = nrowin[0];
    int ncol = ncolin[0];

    int i, j, k;

    for (i = 0; i < ncol; ++i)
        c[i] = 0.0;

    for (j = 0, k = 0; j < ncol; ++j)
        for (i = 0; i < nrow; ++i, ++k)
            c[j] += a[k] * b[i];
}

/* does t(a) %*% diag(b) %*% a where a is matrix b is vector */

void aster_mat_vec_mat_mult(int *nrowin, int *ncolin, double *a, double *b,
    double *c)
{
    int nrow = nrowin[0];
    int ncol = ncolin[0];

    int i, j, k;

    for (i = 0; i < ncol * ncol; ++i)
        c[i] = 0.0;

    for (j = 0; j < ncol; ++j) {
        int jbase = nrow * j;
        for (k = 0; k < ncol; ++k) {
            int kbase = nrow * k;
            int m = j + k * ncol;
            for (i = 0; i < nrow; ++i)
                c[m] += a[jbase + i] * a[kbase + i] * b[i];
        }
    }
}

/* does diag[ a %*% b %*% t(a) ] where a and b are matrices */
/* nrow and ncol are dimensions of a -- hence b is square ncol x ncol */

void aster_diag_mat_mat_mat_mult(int *nrowin, int *ncolin, double *a,
    double *b, double *c)
{
    int nrow = nrowin[0];
    int ncol = ncolin[0];

    int i, j, k;

    for (i = 0; i < nrow; ++i) {
        c[i] = 0.0;
        for (j = 0; j < ncol; ++j) {
            int jbase = nrow * j;
            for (k = 0; k < ncol; ++k) {
                int kbase_for_a = nrow * k;
                int kbase_for_b = ncol * k;
                c[i] += a[jbase + i] * a[kbase_for_a + i] * b[kbase_for_b + j];
            }
        }
    }
}

