
#ifndef ASTER_ASTER_H
#define ASTER_ASTER_H

void aster_family(int *famin, int *derivin, double *thetain, double *value);

int aster_family_validate(int fam, double x, double xpred);

int aster_family_validate_parent(int fam, double xpred);

int aster_family_number_validate(int fam);

double aster_family_origin(int fam);

double aster_family_simulate(int fam, double xpred, double theta);

#ifdef ASTER_OLD_STUFF
char *aster_family_name(int fam);
#endif /* ASTER_OLD_STUFF */

void aster_default_origin(int *nindin, int *nnodein, int *fam, double *result);

void aster_check_model(int *nindin, int *nnodein, int *pred, int *fam);

void aster_check_model_data(int *nindin, int *nnodein, int *pred, int *fam,
    double *x, double *root);

void aster_check_model_root(int *nindin, int *nnodein, int *pred, int *fam,
    double *root);

void aster_theta2phi(int *nindin, int *nnodein, int *pred, int *fam,
    double *theta, double *phi);

void aster_phi2theta(int *nindin, int *nnodein, int *pred, int *fam,
    double *phi, double *theta);

void aster_theta2whatsis(int *nindin, int *nnodein, int *pred, int *fam,
    int *derivin, double *theta, double *result);

void aster_theta2ctau(int *nindin, int *nnodein, int *pred, int *fam,
    double *theta, double *ctau);

void aster_xpred(int *nindin, int *nnodein, int *pred, int *fam,
    double *x, double *root, double *xpred);

void aster_ctau2tau(int *nindin, int *nnodein, int *pred, int *fam,
    double *root, double *ctau, double *tau);

void aster_tt2var(int *nindin, int *nnodein, int *pred, int *fam,
    double *x, double *root, double *theta, double *tau, double *var);

#ifdef ASTER_OLD_STUFF
void aster_unco_hess(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *ctau, double *varx, double *modmat, double *hessian);
#endif /* ASTER_OLD_STUFF */

void aster_a_delsqpsi_m(int *nindin, int *nnodein, int *ncoefin,
    int *ncoef_for_a, int *pred, int *fam, double *ctau, double *varx,
    double *a, double *m, double *result);

void aster_id_delsqpsi_m(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *ctau, double *varx, double *m, double *result);

#ifdef ASTER_OLD_STUFF
void aster_rnzp(int *nin, int *len_xpred_in, int *len_mu_in,
    double *xpred, double *mu, double *result);
#endif /* ASTER_OLD_STUFF */

void aster_rktp(int *nin, int *len_xpred_in, int *len_mu_in, int *len_k_in,
    double *xpred_in, double *mu_in, int *k_in, double *result);

void aster_rktnb(int *nin, int *len_xpred_in, int *len_mu_in, int *len_k_in,
    int *len_alpha_in, double *xpred_in, double *mu_in, int *k_in,
    double *alpha_in, double *result);

void aster_simulate_data(int *nindin, int *nnodein, int *pred, int *fam,
    double *theta, double *root, double *x);

void aster_mlogl_cond(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, int *derivin, double *beta, double *root, double *x,
    double *origin,
    double *modmat, double *value, double *gradient, double *hessian);

void aster_mlogl_unco(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, int *derivin, double *beta, double *root, double *x,
    double *origin,
    double *modmat, double *value, double *gradient, double *hessian);

void aster_fish_cond(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *beta, double *root, double *x,
    double *modmat, double *fish);

void aster_mat_vec_mult(int *nrowin, int *ncolin, double *a, double *b,
    double *c);

void aster_vec_mat_mult(int *nrowin, int *ncolin, double *a, double *b,
    double *c);

void aster_mat_vec_mat_mult(int *nrowin, int *ncolin, double *a, double *b,
    double *c);

void aster_diag_mat_mat_mat_mult(int *nrowin, int *ncolin, double *a,
    double *b, double *c);

void aster_D_beta2theta2phi(int *nindin, int *nnodein, int *ncoefin,
    int *pred, int *fam, double *theta, double *modmat, double *gradmat);

void aster_D_beta2phi2theta(int *nindin, int *nnodein, int *ncoefin,
    int *pred, int *fam, double *theta, double *modmat, double *gradmat);

void aster_D_beta2phi2tau(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *beta, double *root, double *origin, double *modmat,
    double *gradmat);

void aster_D_beta2theta2tau(int *nindin, int *nnodein, int *ncoefin, int *pred,
    int *fam, double *beta, double *root, double *modmat, double *gradmat);

void aster_clear_families(void);

void aster_add_family(char **name, double *hyper, int *nhyper);

void aster_get_family(int *famin, char **name, double *hyper, int *nhyper,
    char **hypername, double *origin);

void aster_get_superfamily(int *famin, char **name, int *nhyper,
    char **hypername);

void aster_byname_superfamily(char **name, int *nhyper, char **hypername);

#endif /* ASTER_ASTER_H */

