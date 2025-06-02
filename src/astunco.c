
#include "aster.h"
#include "raster.h"

void aster_mlogl_unco(int *nindin, int *nnodein, int *ncoefin, int *pred,
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

    int i;
    double *phi, *theta, *xpred, *psi, *psi_prime, *tau, *gfull, *varx;

    aster_check_model_data(nindin, nnodein, pred, fam, x, root);

    phi = (double *) my_malloc(ndata * sizeof(double));
    theta = (double *) my_malloc(ndata * sizeof(double));
    xpred = (double *) my_malloc(ndata * sizeof(double));
    psi = (double *) my_malloc(ndata * sizeof(double));

    aster_mat_vec_mult(&ndata, &ncoef, modmat, beta, phi);
    for (i = 0; i < ndata; ++i)
        phi[i] += origin[i];
    aster_phi2theta(&nind, &nnode, pred, fam, phi, theta);
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
        tau = (double *) my_malloc(ndata * sizeof(double));
        gfull = (double *) my_malloc(ndata * sizeof(double));
        aster_theta2whatsis(&nind, &nnode, pred, fam, &one, theta, psi_prime);
        aster_ctau2tau(&nind, &nnode, pred, fam, root, psi_prime, tau);
        for (i = 0; i < ndata; ++i)
            gfull[i] = - (x[i] - tau[i]);
        aster_vec_mat_mult(&ndata, &ncoef, modmat, gfull, gradient);

        if (deriv >= 2) {
            varx = (double *) my_malloc(ndata * sizeof(double));
            aster_tt2var(&nind, &nnode, pred, fam, x, root, theta, tau, varx);
#ifdef ASTER_OLD_STUFF
            aster_unco_hess(&nind, &nnode, &ncoef, pred, fam,
                psi_prime, varx, modmat, hessian);
#else
            aster_a_delsqpsi_m(&nind, &nnode, &ncoef, &ncoef, pred, fam,
                psi_prime, varx, modmat, modmat, hessian);
#endif /* ASTER_OLD_STUFF */
            my_free(varx);
        }

        my_free(gfull);
        my_free(tau);
        my_free(psi_prime);
    }

    my_free(psi);
    my_free(xpred);
    my_free(theta);
    my_free(phi);
}

