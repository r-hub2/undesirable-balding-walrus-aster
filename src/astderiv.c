
#include "aster.h"
#include "raster.h"

void aster_D_beta2theta2phi(int *nindin, int *nnodein, int *ncoefin,
    int *pred, int *fam, double *theta, double *modmat, double *gradmat)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int ncoef = ncoefin[0];

    int i, j, k, m, jfam, jbase;

    aster_check_model(nindin, nnodein, pred, fam);

    for (i = 0; i < nind * nnode * ncoef; ++i)
        gradmat[i] = modmat[i];

    /* j and k are 1-origin indexing */
    for (j = nnode; j > 0; --j) {
        k = pred[j - 1];
        jfam = fam[j - 1];
        if (k > 0) {
            jbase = nind * (j - 1);
            for (i = 0; i < nind; ++i) {
                int one = 1;
                double foo;
                aster_family(&jfam, &one, &theta[jbase + i], &foo);
                for (m = 0; m < ncoef; ++m)
                    gradmat[nind * (nnode * m + (k - 1)) + i] -=
                    modmat[nind * (nnode * m + (j - 1)) + i] * foo;
            }
        }
    }
}

