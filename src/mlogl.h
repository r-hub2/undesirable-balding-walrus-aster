
#ifndef ASTER_MLOGL_H
#define ASTER_MLOGL_H

/* C function to evaluate saturated model aster log likelihoods
 * like
 *     aster_mlogl_unco
 *     aster_mlogl_cond
 * except those evaluate submodels specified by "origin" (offset) and model
 * matrix, and we don't want that.
 *
 * also, for now, we don't want to bother with derivatives (maybe add later)
 * because right now we are doing Bayesian.
 */

double aster_mlogl_sat_unco(int nind, int nnode, int *pred, int *fam,
    double *phi, double *root, double *response, _Bool check);

double aster_mlogl_sat_cond(int nind, int nnode, int *pred, int *fam,
    double *theta, double *root, double *response, _Bool check);

/* arguments:
 *     nind      number of individuals in the data
 *     nnode     number of nodes in the graph for one individual (all
 *                   individuals have same graph)
 *     pred      integer vector of length nnode, pred[i] - 1 is index of
 *                   predecessor of node i if pred[i] > 0.  Otherwise
 *     fam       integer vector of length nnode, fam[i] is the number of
 *                   the family for the arrow from the predecessor of node i
 *                   to node i
 *     phi       double vector of length nind * nnode, the saturated model
 *                   unconditional canonical parameter vector; the order is
 *                   first node of the graph for all individuals, second
 *                   node of the graph for all individuals (in the same order
 *                   as before), and so forth.
 *     theta     like phi except conditional rather than unconditional
 *                   canonical parameter vector
 *     root      like phi except data at root nodes, for individual i and node j
 *                   (zero-origin indexing) k = i + j * nind is the index
 *                   into root for this individual and node; if pred[j] == 0,
 *                   then root[k] is the data for the predecessor of this node
 *                   for this individual; if pred[j] > 0, then root[k] is
 *                   ignored
 *     response  like root except response data; if i, j, k are as above, then
 *                   the data at node j for individual i is response[k]; and
 *                   if pred[j] > 0, the predecessor node data is
 *                   response[i + (pred[j] - 1) * nind]
 *     check     if true check for validity, otherwise go for speed
 */

void aster_export_exerciser(int *nind_in, int *nnode_in, int *pred, int *fam,
    double *parm, double *root, double *response, int *is_unco_in,
    double *value);

#endif /* ASTER_MLOGL_H */

