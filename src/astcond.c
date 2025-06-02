
#include "aster.h"
#include "raster.h"

void aster_mlogl_cond(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, int *derivin, double *beta, double *root, double *x,
    double *origin,
    double *modmat, double *value, double *gradient, double *hessian)
{
    int nind = nindin[0];
    int nnode = nnodein[0];
    int ncoef = ncoefin[0];
    int deriv = derivin[0];
    int ndata = nind * nnode;
    int zero = 0;
    int one = 1;
    int two = 2;

    int i;
    double *theta, *xpred, *psi, *psi_prime, *gfull,
        *psi_prime_prime, *hfulldiag;

    aster_check_model_data(nindin, nnodein, pred, fam, x, root);

    theta = (double *) my_malloc(ndata * sizeof(double));
    xpred = (double *) my_malloc(ndata * sizeof(double));
    psi = (double *) my_malloc(ndata * sizeof(double));

    aster_mat_vec_mult(&ndata, &ncoef, modmat, beta, theta);
    for (i = 0; i < ndata; ++i)
        theta[i] += origin[i];
    aster_xpred(&nind, &nnode, pred, fam, x, root, xpred);
    aster_theta2whatsis(&nind, &nnode, pred, fam, &zero, theta, psi);

    value[0] = 0.0;
    for (i = 0; i < ndata; ++i)
        value[0] -= x[i] * theta[i] - xpred[i] * psi[i];

    if (my_is_na_or_nan(value[0]))
        value[0] = my_posinf();
    if (value[0] == my_neginf())
        die("calculated log likelihood + infinity, impossible");

    if (deriv >= 1 && value[0] < my_posinf()) {
        psi_prime = (double *) my_malloc(ndata * sizeof(double));
        gfull = (double *) my_malloc(ndata * sizeof(double));
        aster_theta2whatsis(&nind, &nnode, pred, fam, &one, theta, psi_prime);
        for (i = 0; i < ndata; ++i)
            gfull[i] = - (x[i] - xpred[i] * psi_prime[i]);
        aster_vec_mat_mult(&ndata, &ncoef, modmat, gfull, gradient);
        my_free(gfull);
        my_free(psi_prime);
    }

    if (deriv >= 2 && value[0] < my_posinf()) {
        psi_prime_prime = (double *) my_malloc(ndata * sizeof(double));
        hfulldiag = (double *) my_malloc(ndata * sizeof(double));
        aster_theta2whatsis(&nind, &nnode, pred, fam, &two, theta,
            psi_prime_prime);
        for (i = 0; i < ndata; ++i)
            hfulldiag[i] = xpred[i] * psi_prime_prime[i];
        aster_mat_vec_mat_mult(&ndata, &ncoef, modmat, hfulldiag, hessian);
        my_free(hfulldiag);
        my_free(psi_prime_prime);
    }

    my_free(psi);
    my_free(xpred);
    my_free(theta);
}

