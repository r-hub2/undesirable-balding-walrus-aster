
#include "aster.h"
#include "raster.h"

void aster_default_origin(int *nindin, int *nnodein, int *fam, double *result)
{
    int nind = nindin[0];
    int nnode = nnodein[0];

    int i, j, k;

    if (nind <= 0)
        die("'nind' must be positive integer\n");
    if (nnode <= 0)
        die("'nnode' must be positive integer\n");

    /* fam values are 1-origin indexing */
    for (j = 0, k = 0; j < nnode; ++j) {
        double foo = aster_family_origin(fam[j]);
        for (i = 0; i < nind; ++i, ++k)
            result[k] = foo;
    }
}

