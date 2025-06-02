
#ifndef ASTER_MLOGL_EXPORT_H
#define ASTER_MLOGL_EXPORT_H

// pointer to either of the functions
//     aster_mlogl_sat_unco
//     aster_mlogl_sat_cond
// defined in mlogl.h in the src directory

typedef double (*aster_mlogl_sat_either_funptr)(int nind, int nnode,
    int *pred, int *fam, double *phi, double *root, double *response,
    _Bool check);

#endif /* ASTER_MLOGL_EXPORT_H */

