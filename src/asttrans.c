
#include "aster.h"

void aster_theta2phi(int *nindin, int *nnodein, int *pred, int *fam,
    double *theta, double *phi)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, k, jfam, jbase, kbase;

    aster_check_model(nindin, nnodein, pred, fam);

    for (i = 0; i < nind * nnode; ++i)
        phi[i] = theta[i];

    /* j and k are 1-origin indexing */
    for (j = nnode; j > 0; --j) {
        k = pred[j - 1];
        jfam = fam[j - 1];
        if (k > 0) {
            jbase = nind * (j - 1);
            kbase = nind * (k - 1);
            for (i = 0; i < nind; ++i) {
                int zero = 0;
                double foo;
                aster_family(&jfam, &zero, &theta[jbase + i], &foo);
                phi[kbase + i] -= foo;
            }
        }
    }
}

