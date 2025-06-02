
// See Section 5.4 of Writing R Extensions
//
// Remove styles, which are deprecated in R-3.3.3
//
// See also Section 6.15 of Writing R Extensions

#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include "aster.h"
#ifdef ASTER_OLD_STUFF
#include "rraster.h"
#endif /* ASTER_OLD_STUFF */
#include "mlogl.h"
#include "objfun.h"
#include "matops2.h"

static R_NativePrimitiveArgType ast_fam_types[4] =
    {INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_chkmod_types[4] =
    {INTSXP, INTSXP, INTSXP, INTSXP};

static R_NativePrimitiveArgType ast_chkmodd_types[6] =
    {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_th2ph_types[6] =
    {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_whatsis_types[7] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_ctau2tau_types[7] =
    {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_tt2var_types[9] =
    {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_xpred_types[7] =
    {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};

#ifdef ASTER_OLD_STUFF
static R_NativePrimitiveArgType ast_rnzp_types[6] =
    {INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};
#endif /* ASTER_OLD_STUFF */

static R_NativePrimitiveArgType ast_rktp_types[8] =
    {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType ast_rktnb_types[10] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    INTSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_mlogl_types[14] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_fish_types[10] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP,
    REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_matops_types[5] =
    {INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_b2t2p_types[8] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_b2p2t_types[10] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
    REALSXP, REALSXP};

static R_NativePrimitiveArgType ast_b2t2t_types[9] =
    {INTSXP, INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP,
    REALSXP};

static R_NativePrimitiveArgType ast_addfam_types[3] =
    {STRSXP, REALSXP, INTSXP};

static R_NativePrimitiveArgType ast_getfam_types[6] =
    {INTSXP, STRSXP, REALSXP, INTSXP, STRSXP, REALSXP};

static R_NativePrimitiveArgType ast_getsup_types[4] =
    {INTSXP, STRSXP, INTSXP, STRSXP};

static R_NativePrimitiveArgType ast_bysup_types[3] = {STRSXP, INTSXP, STRSXP};

static R_NativePrimitiveArgType ast_origin_types[4] =
    {INTSXP, INTSXP, INTSXP, REALSXP};

static R_NativePrimitiveArgType ast_export_exerciser_types[9] =
    {INTSXP, INTSXP, INTSXP, INTSXP, REALSXP, REALSXP, REALSXP, LGLSXP,
    REALSXP};

static R_CMethodDef cMethods[] = {
    {"aster_family", (DL_FUNC) &aster_family, 4, ast_fam_types},
    {"aster_check_model", (DL_FUNC) &aster_check_model, 4, ast_chkmod_types},
    {"aster_check_model_data", (DL_FUNC) &aster_check_model_data,
        6, ast_chkmodd_types},
    {"aster_theta2phi", (DL_FUNC) &aster_theta2phi, 6, ast_th2ph_types},
    {"aster_phi2theta", (DL_FUNC) &aster_phi2theta, 6, ast_th2ph_types},
    {"aster_theta2whatsis", (DL_FUNC) &aster_theta2whatsis,
        7, ast_whatsis_types},
    {"aster_theta2ctau", (DL_FUNC) &aster_theta2ctau, 6, ast_th2ph_types},
    {"aster_xpred", (DL_FUNC) &aster_xpred, 7, ast_xpred_types},
    {"aster_ctau2tau", (DL_FUNC) &aster_ctau2tau, 7, ast_ctau2tau_types},
    {"aster_tt2var", (DL_FUNC) &aster_tt2var, 9, ast_tt2var_types},
#ifdef ASTER_OLD_STUFF
    {"aster_rnzp", (DL_FUNC) &aster_rnzp, 6, ast_rnzp_types},
#endif /* ASTER_OLD_STUFF */
    {"aster_rktp", (DL_FUNC) &aster_rktp, 8, ast_rktp_types},
    {"aster_rktnb", (DL_FUNC) &aster_rktnb, 10, ast_rktnb_types},
    {"aster_simulate_data", (DL_FUNC) &aster_simulate_data, 7, ast_xpred_types},
    {"aster_mlogl_cond", (DL_FUNC) &aster_mlogl_cond, 14, ast_mlogl_types},
    {"aster_mlogl_unco", (DL_FUNC) &aster_mlogl_unco, 14, ast_mlogl_types},
    {"aster_fish_cond", (DL_FUNC) &aster_fish_cond, 10, ast_fish_types},
    {"aster_mat_vec_mult", (DL_FUNC) &aster_mat_vec_mult, 5, ast_matops_types},
    {"aster_vec_mat_mult", (DL_FUNC) &aster_vec_mat_mult, 5, ast_matops_types},
    {"aster_mat_vec_mat_mult", (DL_FUNC) &aster_mat_vec_mat_mult,
        5, ast_matops_types},
    {"aster_diag_mat_mat_mat_mult", (DL_FUNC) &aster_diag_mat_mat_mat_mult,
        5, ast_matops_types},
    {"aster_D_beta2theta2phi", (DL_FUNC) &aster_D_beta2theta2phi,
        8, ast_b2t2p_types},
    {"aster_D_beta2phi2theta", (DL_FUNC) &aster_D_beta2phi2theta,
        8, ast_b2t2p_types},
    {"aster_D_beta2phi2tau", (DL_FUNC) &aster_D_beta2phi2tau,
        10, ast_b2p2t_types},
    {"aster_D_beta2theta2tau", (DL_FUNC) &aster_D_beta2theta2tau,
        9, ast_b2t2t_types},
    {"aster_clear_families", (DL_FUNC) &aster_clear_families, 0, NULL},
    {"aster_add_family", (DL_FUNC) &aster_add_family, 3, ast_addfam_types},
    {"aster_get_family", (DL_FUNC) &aster_get_family, 6, ast_getfam_types},
    {"aster_get_superfamily", (DL_FUNC) &aster_get_superfamily,
        4, ast_getsup_types},
    {"aster_byname_superfamily", (DL_FUNC) &aster_byname_superfamily,
        3, ast_bysup_types},
    {"aster_default_origin", (DL_FUNC) &aster_default_origin,
        4, ast_origin_types},
    {"aster_export_exerciser", (DL_FUNC) &aster_export_exerciser,
        9, ast_export_exerciser_types},
    {NULL, NULL, 0, NULL}
};
 
static R_CallMethodDef callMethods[]  = {
#ifdef ASTER_OLD_STUFF
    {"aster_families", (DL_FUNC) &aster_families, 0},
#endif /* ASTER_OLD_STUFF */
    {"objfun", (DL_FUNC) &objfun, 12},
    {"pos_def_mat_inv", (DL_FUNC) &pos_def_mat_inv, 1},
    {NULL, NULL, 0}
};

void attribute_visible R_init_aster(DllInfo *info)
{
    R_registerRoutines(info, cMethods, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
    R_RegisterCCallable("aster", "aster_mlogl_sat_unco",
        (DL_FUNC) aster_mlogl_sat_unco);
    R_RegisterCCallable("aster", "aster_mlogl_sat_cond",
        (DL_FUNC) aster_mlogl_sat_cond);
}

