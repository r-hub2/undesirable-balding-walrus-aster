
#include "aster.h"
#include "raster.h"

void aster_D_beta2theta2tau(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *beta, double *root, double *modmat, double *gradmat)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int ncoef = ncoefin[0];
    int ndata = nind * nnode;
    int one = 1;
    int two = 2;

    int i, j1, j2, k;
    double *theta, *psi_prime, *tau, *psi_prime_prime;

    aster_check_model(nindin, nnodein, pred, fam);

    for (i = 0; i < nind * nnode * ncoef; ++i)
        gradmat[i] = 0.0;

    theta = (double *) my_malloc(ndata * sizeof(double));
    psi_prime = (double *) my_malloc(ndata * sizeof(double));
    tau = (double *) my_malloc(ndata * sizeof(double));
    psi_prime_prime = (double *) my_malloc(ndata * sizeof(double));

    aster_mat_vec_mult(&ndata, &ncoef, modmat, beta, theta);
    aster_theta2whatsis(&nind, &nnode, pred, fam, &one, theta, psi_prime);
    aster_ctau2tau(&nind, &nnode, pred, fam, root, psi_prime, tau);
    aster_theta2whatsis(&nind, &nnode, pred, fam, &two, theta,
            psi_prime_prime);

    /* j1 and j2 are 1-origin indexing */
    for (j1 = nnode; j1 > 0; --j1) {
        j2 = j1;
        while (j2 > 0) {
            for (i = 0; i < nind; ++i) {
                double foo = tau[nind * (j1 - 1) + i] *
                    psi_prime_prime[nind * (j2 - 1) + i] /
                    psi_prime[nind * (j2 - 1) + i];
                for (k = 0; k < ncoef; ++k)
                    gradmat[nind * (nnode * k + (j1 - 1)) + i] +=
                        foo * modmat[nind * (nnode * k + (j2 - 1)) + i];
            }
            j2 = pred[j2 - 1];
        }
    }

    my_free(psi_prime_prime);
    my_free(tau);
    my_free(psi_prime);
    my_free(theta);
}

