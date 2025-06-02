

#include "mlogl.h"
#include "aster.h"
#include "raster.h"

double aster_mlogl_sat_unco(int nind, int nnode, int *pred, int *fam,
    double *phi, double *root, double *response, _Bool check)
{
    if (check)
        aster_check_model_data(&nind, &nnode, pred, fam, response, root);

    int ndata = nind * nnode;
    double *theta = (double *) my_malloc(ndata * sizeof(double));
    aster_phi2theta(&nind, &nnode, pred, fam, phi, theta);
    double value = aster_mlogl_sat_cond(nind, nnode, pred, fam,
        theta, root, response, 0);
    // there is memory leak here if aster_mlogl_sat_cond does not return
    my_free(theta);
    return value;
}

double aster_mlogl_sat_cond(int nind, int nnode, int *pred, int *fam,
    double *theta, double *root, double *response, _Bool check)
{
    if (check)
        aster_check_model_data(&nind, &nnode, pred, fam, response, root);

    int ndata = nind * nnode;
    double *xpred = (double *) my_malloc(ndata * sizeof(double));
    double *cumfun = (double *) my_malloc(ndata * sizeof(double));
    int zero = 0;

    aster_xpred(&nind, &nnode, pred, fam, response, root, xpred);
    aster_theta2whatsis(&nind, &nnode, pred, fam, &zero, theta, cumfun);

    double value = 0.0;
    for (int i = 0; i < ndata; ++i)
        value -= response[i] * theta[i] - xpred[i] * cumfun[i];

    my_free(cumfun);
    my_free(xpred);

    if (my_is_na_or_nan(value))
        value = my_posinf();
    if (value == my_neginf()) {
        die("calculated log likelihood is +infinity, impossible");
    }

    return value;
}

// We don't use these function pointers.  They are just here to make
// the compiler complain if the typedef is wrong.
// And we only need the typedef for the benefit of functions in other
// R packages calling these functions via the R_GetCCallable mechanism.

#include "mlogl-export.h"
#if defined(__GNUC__) || defined(__clang__)
static aster_mlogl_sat_either_funptr foo __attribute__ ((unused))
    = aster_mlogl_sat_unco;
static aster_mlogl_sat_either_funptr bar __attribute__ ((unused))
    = aster_mlogl_sat_cond;
#else
static aster_mlogl_sat_either_funptr foo = aster_mlogl_sat_unco;
static aster_mlogl_sat_either_funptr bar = aster_mlogl_sat_cond;
#endif /* defined(__GNUC__) || defined(__clang__) */

// Except.  We are going to use this typedef below in code
// that looks like the code that would be in another package
// calling this package.  The following function is not used
// anywhere in this package except in the file mlogl.R in the
// tests directory.

#include <stddef.h>
#include <R_ext/Rdynload.h>

void aster_export_exerciser(int *nind_in, int *nnode_in, int *pred, int *fam,
    double *parm, double *root, double *response, int *is_unco_in,
    double *value)
{
    int nind = nind_in[0];
    int nnode = nnode_in[0];
    int is_unco = is_unco_in[0];

    aster_mlogl_sat_either_funptr fun = NULL;
    if (is_unco) {
        fun = (aster_mlogl_sat_either_funptr)
            R_GetCCallable("aster", "aster_mlogl_sat_unco");
    } else {
        fun = (aster_mlogl_sat_either_funptr)
            R_GetCCallable("aster", "aster_mlogl_sat_cond");
    }
    value[0] = fun(nind, nnode, pred, fam, parm, root, response, 1);
}
