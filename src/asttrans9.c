
#include "aster.h"

void aster_tt2var(int *nindin, int *nnodein, int *pred, int *fam,
    double *x, double *root, double *theta, double *tau, double *var)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, k, jfam, jbase, kbase;

    aster_check_model_data(nindin, nnodein, pred, fam, x, root);

    /* j and k are 1-origin indexing */
    for (j = 1; j <= nnode; ++j) {
        k = pred[j - 1];
        jfam = fam[j - 1];
        jbase = nind * (j - 1);
        kbase = nind * (k - 1);
        for (i = 0; i < nind; ++i) {
            int one = 1;
            int two = 2;
            double psi_prime;
            double psi_double_prime;
            aster_family(&jfam, &one, &theta[jbase + i], &psi_prime);
            aster_family(&jfam, &two, &theta[jbase + i], &psi_double_prime);
            if (k > 0)
                var[jbase + i] = psi_double_prime * tau[kbase + i]
                    + psi_prime * psi_prime * var[kbase + i];
            else
                var[jbase + i] = psi_double_prime * root[jbase + i];
        }
    }
}

