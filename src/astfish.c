
#include "aster.h"
#include "raster.h"

void aster_fish_cond(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *beta, double *root, double *x,
    double *modmat, double *fish)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int ncoef = ncoefin[0];
    int ndata = nind * nnode;
    int two = 2;

    int i;
    double *theta, *ctau, *tau, *taupred, *psi_prime_prime, *hfulldiag;

    aster_check_model_data(nindin, nnodein, pred, fam, x, root);

    theta = (double *) my_malloc(ndata * sizeof(double));
    ctau = (double *) my_malloc(ndata * sizeof(double));
    tau = (double *) my_malloc(ndata * sizeof(double));
    taupred = (double *) my_malloc(ndata * sizeof(double));
    psi_prime_prime = (double *) my_malloc(ndata * sizeof(double));
    hfulldiag = (double *) my_malloc(ndata * sizeof(double));

    aster_mat_vec_mult(&ndata, &ncoef, modmat, beta, theta);
    aster_theta2ctau(&nind, &nnode, pred, fam, theta, ctau);
    aster_ctau2tau(&nind, &nnode, pred, fam, root, ctau, tau);
    aster_xpred(&nind, &nnode, pred, fam, tau, root, taupred);
    aster_theta2whatsis(&nind, &nnode, pred, fam, &two, theta,
        psi_prime_prime);
    for (i = 0; i < ndata; ++i)
        hfulldiag[i] = taupred[i] * psi_prime_prime[i];
    aster_mat_vec_mat_mult(&ndata, &ncoef, modmat, hfulldiag, fish);

    my_free(hfulldiag);
    my_free(psi_prime_prime);
    my_free(taupred);
    my_free(tau);
    my_free(ctau);
    my_free(theta);
}

