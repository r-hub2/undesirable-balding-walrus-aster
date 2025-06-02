
#include <stddef.h>
#include "aster.h"
#include "raster.h"

void aster_check_model(int *nindin, int *nnodein, int *pred, int *fam)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int j;

    if (nind <= 0)
        die("'nind' must be positive integer\n");
    if (nnode <= 0)
        die("'nnode' must be positive integer\n");

    /* pred values are 1-origin indexing */
    for (j = 0; j < nnode; ++j)
        if (pred[j] > j)
            /* error report as if using 1-origin indexing */
            die("pred[%d] = %d, invalid\n", j + 1, pred[j]);

    /* fam values are 1-origin indexing */
    for (j = 0; j < nnode; ++j)
        if (! aster_family_number_validate(fam[j]))
            die("family %d not valid", fam[j]);
}

void aster_check_model_data(int *nindin, int *nnodein, int *pred, int *fam,
    double *x, double *root)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, k, jfam, jbase, kbase;

    aster_check_model(nindin, nnodein, pred, fam);

    /* j and k are 1-origin indexing */
    for (j = nnode; j > 0; --j) {
        k = pred[j - 1];
        jfam = fam[j - 1];
        jbase = nind * (j - 1);
        kbase = nind * (k - 1);
        for (i = 0; i < nind; ++i) {
            double xnow = x[jbase + i];
            double xpred;
            if (k > 0)
                xpred = x[kbase + i];
            else
                xpred = root[jbase + i];
            if (! aster_family_validate(jfam, xnow, xpred))
                die("invalid data: family = %d, x = %f, xpred = %f\n",
                    jfam, xnow, xpred);
        }
    }
}

/* needed when predicting tau, which depends on x only at root nodes */

void aster_check_model_root(int *nindin, int *nnodein, int *pred, int *fam,
    double *root)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, k, jfam, jbase;

    aster_check_model(nindin, nnodein, pred, fam);

    /* j and k are 1-origin indexing */
    for (j = nnode; j > 0; --j) {
        k = pred[j - 1];
        jfam = fam[j - 1];
        jbase = nind * (j - 1);
        for (i = 0; i < nind; ++i) {
            if (k == 0) {
                double xpred = root[jbase + i];
                if (! aster_family_validate_parent(jfam, xpred))
                    die("invalid root data: family = %d, xpred = %f\n",
                        jfam, xpred);
            }
        }
    }
}

