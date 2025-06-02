
#include <stddef.h>
#include "aster.h"
#include "raster.h"

void aster_simulate_data(int *nindin, int *nnodein, int *pred, int *fam,
    double *theta, double *root, double *x)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, k, jfam, jbase, kbase;

    aster_check_model(nindin, nnodein, pred, fam);

    my_GetRNGstate();
    /* j and k are 1-origin indexing */
    for (j = 1; j <= nnode; ++j) {
        k = pred[j - 1];
        jfam = fam[j - 1];
        jbase = nind * (j - 1);
        kbase = nind * (k - 1);
        for (i = 0; i < nind; ++i) {
            double xpred;
            if (k > 0)
                xpred = x[kbase + i];
            else
                xpred = root[jbase + i];
            x[jbase + i] = aster_family_simulate(jfam, xpred, theta[jbase + i]);
        }
    }
    my_PutRNGstate();

    aster_check_model_data(nindin, nnodein, pred, fam, x, root);
}

