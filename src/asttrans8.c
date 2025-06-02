
#include "aster.h"
#include "raster.h"

void aster_theta2whatsis(int *nindin, int *nnodein, int *pred, int *fam,
    int *derivin, double *theta, double *result)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int deriv = derivin[0];

    int i, j, jfam, jbase;

    aster_check_model(nindin, nnodein, pred, fam);
    if (deriv < 0 || deriv > 2)
        die("deriv not 0, 1, or 2\n");

    /* j is 1-origin indexing */
    for (j = nnode; j > 0; --j) {
        jfam = fam[j - 1];
        jbase = nind * (j - 1);
        for (i = 0; i < nind; ++i) {
            double foo;
            aster_family(&jfam, &deriv, &theta[jbase + i], &foo);
            result[jbase + i] = foo;
        }
    }
}

