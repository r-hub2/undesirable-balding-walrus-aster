
#include "aster.h"
#include "raster.h"

void aster_D_beta2phi2tau(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *beta, double *root, double *origin, double *modmat,
    double *gradmat)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int ncoef = ncoefin[0];
    int ndata = nind * nnode;

    int i;
    double *phi, *theta, *psi_prime, *tau, *varx;

    aster_check_model(nindin, nnodein, pred, fam);

    for (i = 0; i < nind * nnode * ncoef; ++i)
        gradmat[i] = modmat[i];

    phi = (double *) my_malloc(ndata * sizeof(double));
    theta = (double *) my_malloc(ndata * sizeof(double));
    psi_prime = (double *) my_malloc(ndata * sizeof(double));
    tau = (double *) my_malloc(ndata * sizeof(double));
    varx = (double *) my_malloc(ndata * sizeof(double));

    aster_mat_vec_mult(&ndata, &ncoef, modmat, beta, phi);
    for (i = 0; i < ndata; ++i)
        phi[i] += origin[i];
    aster_phi2theta(&nind, &nnode, pred, fam, phi, theta);
    aster_theta2ctau(&nind, &nnode, pred, fam, theta, psi_prime);
    aster_ctau2tau(&nind, &nnode, pred, fam, root, psi_prime, tau);
    aster_tt2var(&nind, &nnode, pred, fam, root, root, theta, tau, varx);
    aster_id_delsqpsi_m(&nind, &nnode, &ncoef, pred, fam,
        psi_prime, varx, modmat, gradmat);

    my_free(varx);
    my_free(tau);
    my_free(psi_prime);
    my_free(theta);
    my_free(phi);
}

