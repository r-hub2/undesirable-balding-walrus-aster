
#include "aster.h"

void aster_theta2ctau(int *nindin, int *nnodein, int *pred, int *fam,
    double *theta, double *ctau)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, jfam, jbase;

    aster_check_model(nindin, nnodein, pred, fam);

    /* j is 1-origin indexing */
    for (j = nnode; j > 0; --j) {
        jfam = fam[j - 1];
        jbase = nind * (j - 1);
        for (i = 0; i < nind; ++i) {
            int one = 1;
            double foo;
            aster_family(&jfam, &one, &theta[jbase + i], &foo);
            ctau[jbase + i] = foo;
        }
    }
}

