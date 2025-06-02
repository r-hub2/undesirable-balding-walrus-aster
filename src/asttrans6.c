
#include "aster.h"

void aster_xpred(int *nindin, int *nnodein, int *pred, int *fam,
    double *x, double *root, double *xpred)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, k, jbase, kbase;

    aster_check_model(nindin, nnodein, pred, fam);

    /* j and k are 1-origin indexing */
    for (j = nnode; j > 0; --j) {
        k = pred[j - 1];
        jbase = nind * (j - 1);
        kbase = nind * (k - 1);
        for (i = 0; i < nind; ++i)
            if (k > 0)
                xpred[jbase + i] = x[kbase + i];
            else
                xpred[jbase + i] = root[jbase + i];
    }
}

