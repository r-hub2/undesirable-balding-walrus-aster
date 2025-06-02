
#include "aster.h"
#include "raster.h"

static double covxx(int i, int j, int j2, int nind, int nnode, int *pred,
    double *ctau, double *varx)
{
    int jbase = nind * (j - 1);

    /* totally one-origin indexing (including i) */
    if (i < 1 || i > nind)
        die("covxx: i = %d out of range\n", i);
    if (j < 1 || j > nnode)
        die("covxx: j = %d out of range\n", j);
    if (j2 < 1 || j2 > nnode)
        die("covxx: j2 = %d out of range\n", j2);

    if (j == j2) {
        return varx[jbase + (i - 1)];
    } else if (j < j2) {
        return covxx(i, j2, j, nind, nnode, pred, ctau, varx);
    } else /* j > j2 */ {
        int k = pred[j - 1];
        if (k > 0) {
            return ctau[jbase + (i - 1)] *
                covxx(i, k, j2, nind, nnode, pred, ctau, varx);
        } else {
            return 0.0;
        }
    }
}

#ifdef ASTER_OLD_STUFF
void aster_unco_hess(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *ctau, double *varx, double *modmat, double *hessian)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int ncoef = ncoefin[0];
    int ndata = nind * nnode;

    int i, j1, j2, m1, m2;

    aster_check_model(nindin, nnodein, pred, fam);

    /* i, j2, and j2 are 1-origin indexing */
    for (j1 = 1; j1 <= nnode; ++j1) {
        int j1base = nind * (j1 - 1);
        for (j2 = 1; j2 <= nnode; ++j2) {
            int j2base = nind * (j2 - 1);
            for (i = 1; i <= nind; ++i) {
                double foo = covxx(i, j1, j2, nind, nnode, pred, ctau, varx);
                for (m1 = 0; m1 < ncoef; ++m1) {
                    int m1base = ndata * m1;
                    for (m2 = 0; m2 < ncoef; ++m2) {
                        int m2base = ndata * m2;
                        hessian[m1 + ncoef * m2] +=
                            modmat[j1base + (i - 1) + m1base] *
                            modmat[j2base + (i - 1) + m2base] *
                            foo;
                    }
                }
            }
        }
    }
}
#endif /* ASTER_OLD_STUFF */

void aster_a_delsqpsi_m(int *nindin, int *nnodein, int *ncoefin,
    int *ncoef_for_ain, int *pred, int *fam, double *ctau, double *varx,
    double *a, double *m, double *result)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int ncoef = ncoefin[0];
    int ncoef_for_a = ncoef_for_ain[0];
    int ndata = nind * nnode;
    int ncoefsq = ncoef * ncoef;

    int i, j1, j2, m1, m2;

    aster_check_model(nindin, nnodein, pred, fam);

    for (i = 0; i < ncoefsq; ++i)
        result[i] = 0.0;

    /* i, j2, and j2 are 1-origin indexing */
    for (j1 = 1; j1 <= nnode; ++j1) {
        int j1base = nind * (j1 - 1);
        for (j2 = 1; j2 <= nnode; ++j2) {
            int j2base = nind * (j2 - 1);
            for (i = 1; i <= nind; ++i) {
                double foo = covxx(i, j1, j2, nind, nnode, pred, ctau, varx);
                for (m1 = 0; m1 < ncoef_for_a; ++m1) {
                    int m1base = ndata * m1;
                    for (m2 = 0; m2 < ncoef; ++m2) {
                        int m2base = ndata * m2;
                        result[m1 + ncoef * m2] +=
                            a[j1base + (i - 1) + m1base] *
                            m[j2base + (i - 1) + m2base] *
                            foo;
                    }
                }
            }
        }
    }
}

/* special case a is identity matrix */
void aster_id_delsqpsi_m(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *ctau, double *varx, double *m, double *result)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int ncoef = ncoefin[0];
    int ndata = nind * nnode;
    int nresult = nind * nnode * ncoef;

    int i, j1, j2, k;

    aster_check_model(nindin, nnodein, pred, fam);

    for (i = 0; i < nresult; ++i)
        result[i] = 0.0;

    /* i, j1, j2, and k are 1-origin indexing */
    for (j1 = 1; j1 <= nnode; ++j1) {
        int j1base = nind * (j1 - 1);
        for (j2 = 1; j2 <= nnode; ++j2) {
            int j2base = nind * (j2 - 1);
            for (i = 1; i <= nind; ++i) {
                double foo = covxx(i, j1, j2, nind, nnode, pred, ctau, varx);
                for (k = 1; k <= ncoef; ++k) {
                    int kbase = ndata * (k - 1);
                    result[kbase + j1base + (i - 1)] += foo *
                        m[kbase + j2base + (i - 1)];
                }
            }
        }
    }
}

