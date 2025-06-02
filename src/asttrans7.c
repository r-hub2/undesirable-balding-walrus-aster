
#include "aster.h"

void aster_ctau2tau(int *nindin, int *nnodein, int *pred, int *fam,
    double *root, double *ctau, double *tau)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, k, jbase, kbase;

    aster_check_model_root(nindin, nnodein, pred, fam, root);

    /* j and k are 1-origin indexing */
    for (j = 1; j <= nnode; ++j) {
        k = pred[j - 1];
        jbase = nind * (j - 1);
        kbase = nind * (k - 1);
        for (i = 0; i < nind; ++i) {
            double foo = ctau[jbase + i];
            if (k > 0)
                foo *= tau[kbase + i];
            else
                foo *= root[jbase + i];
            tau[jbase + i] = foo;
        }
    }
}

