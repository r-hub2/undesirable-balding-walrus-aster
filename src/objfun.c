#include <R.h>
#include <Rinternals.h>
#include "objfun.h"
#include "aster.h"
#include "matops2.h"
#include <string.h>
#include <math.h>
#include <stddef.h>

SEXP objfun(SEXP theta, SEXP modmat, SEXP nfixed, SEXP nrand,
    SEXP response, SEXP pred, SEXP fam, SEXP root, SEXP zwz,
    SEXP offset, SEXP standard_deviation, SEXP deriv)
{
    if (! isInteger(nfixed))
        error("argument nfixed must be storage mode integer");
    if (LENGTH(nfixed) != 1)
        error("argument nfixed must be length one");
    int len_alpha = INTEGER(nfixed)[0];
    if (len_alpha <= 0)
        error("argument nfixed must be positive-integer-valued");

    if (! isInteger(nrand))
        error("argument nrand must be storage mode integer");
    int *my_nrand = INTEGER(nrand);
    int len_nu = LENGTH(nrand);
    int len_bee = 0;
    for (int i = 0; i < len_nu; i++) {
        if (my_nrand[i] <= 0)
            error("argument nrand must be positive-integer-valued");
        len_bee += my_nrand[i];
    }

    if (! isMatrix(modmat))
        error("argument modmat must be matrix");
    int nfull = nrows(modmat); // number of nodes in full aster graph
    if (! isReal(modmat))
        error("argument modmat must be storage mode double");
    double *my_modmat = REAL(modmat);
    if (ncols(modmat) != len_alpha + len_bee)
        error("ncols(modmat) != nfixed + sum(nrand)");
    int my_ncols_modmat = len_alpha + len_bee;

    if (! isReal(theta))
        error("argument theta must be storage mode double");
    if (LENGTH(theta) != len_alpha + len_bee + len_nu)
        error("length(theta) != nfixed + sum(nrand) + length(nrand)");

    if (! isReal(response))
        error("argument response must be storage mode double");
    if (LENGTH(response) != nfull)
        error("length(response) != nrow(modmat)");
    double *my_response = REAL(response);

    if (! isReal(root))
        error("argument root must be storage mode double");
    if (LENGTH(root) != nfull)
        error("length(root) != nrow(modmat)");
    double *my_root = REAL(root);

    if (! isReal(offset))
        error("argument offset must be storage mode double");
    if (LENGTH(offset) != nfull)
        error("length(offset) != nrow(modmat)");
    double *my_offset = REAL(offset);

    if (! isInteger(pred))
        error("argument pred must be storage mode integer");
    int *my_pred = INTEGER(pred);
    int nnode = LENGTH(pred);
    int nind = nfull / nnode;
    if (nfull != nnode * nind)
        error("nrow(modmat) not divisible by length(pred)");

    if (! isInteger(fam))
        error("argument fam must be storage mode integer");
    if (LENGTH(fam) != nnode)
        error("length(fam) != length(pred)");
    int *my_fam = INTEGER(fam);

    if (! isLogical(standard_deviation))
        error("argument standard_deviation must be storage mode logical");
    if (LENGTH(standard_deviation) != 1)
        error("argument standard_deviation must have length 1");
    _Bool my_standard_deviation = LOGICAL(standard_deviation)[0];

    if (! isInteger(deriv))
        error("argument deriv must be storage mode integer");
    int my_deriv = INTEGER(deriv)[0];
    if (my_deriv < 0 || my_deriv > 2)
        error("argument deriv not 0, 1, or 2");

    int *dee_idx = (int *) R_alloc(len_bee, sizeof(int));
    for (int k = 0, m = 0; k < len_nu; k++)
        for (int i = 0; i < my_nrand[k]; i++, m++)
            dee_idx[m] = k;

    double *other_theta = (double *) R_alloc(len_alpha + len_bee + len_nu,
        sizeof(double));
    double *dee = (double *) R_alloc(len_bee, sizeof(double));
    double *alpha = REAL(theta);
    double *bee = NULL;
    double *alphabee = NULL;
    double *nu = NULL;
    double *cee = NULL;
    double *sigma = NULL;
    memcpy(other_theta, alpha, len_alpha * sizeof(double));
    if (my_standard_deviation) {
        cee = alpha + len_alpha;
        sigma = cee + len_bee;
        bee = other_theta + len_alpha;
        nu = bee + len_bee;
        for (int k = 0; k < len_nu; k++)
            nu[k] = sigma[k] * sigma[k];
        for (int i = 0; i < len_bee; i++) {
            int k = dee_idx[i];
            dee[i] = nu[k];
            bee[i] = sigma[k] * cee[i];
        }
        alphabee = other_theta;
    } else {
        bee = alpha + len_alpha;
        nu = bee + len_bee;
        cee = other_theta + len_alpha;
        sigma = cee + len_bee;
        for (int k = 0; k < len_nu; k++)
            sigma[k] = sqrt(nu[k]); // assigns NaN if nu[k] < 0.0
        for (int i = 0; i < len_bee; i++) {
            int k = dee_idx[i];
            dee[i] = nu[k];
            cee[i] = bee[i] / sigma[k]; // assigns Inf or NaN if nu[k] <= 0.0
            if (sigma[k] == 0.0)
                cee[i] = 0.0; // arbitrary
        }
        alphabee = REAL(theta);
    }

    if (! isMatrix(zwz))
        error("argument zwz must be matrix or NULL");
    if (nrows(zwz) != len_bee)
        error("nrow(zwz) != sum(nrand)");
    if (ncols(zwz) != len_bee)
        error("ncol(zwz) != sum(nrand)");
    // should check that matrix is symmetric positive definite
    // but takes too much time
    // actually maybe not positive definite or even semidefinite if subsampling
    double *my_zwz = REAL(zwz);

    double value;
    double *gradient = (double *) R_alloc(len_alpha + len_bee, sizeof(double));
    double *hessian = (double *) R_alloc((len_alpha + len_bee) *
        (len_alpha + len_bee), sizeof(double));

    aster_mlogl_unco(&nind, &nnode, &my_ncols_modmat, my_pred, my_fam,
        &my_deriv, alphabee, my_root, my_response, my_offset, my_modmat,
        &value, gradient, hessian);

    for (int i = 0; i < len_bee; i++) {
        double foo = cee[i] * cee[i] / 2.0;
        if (! R_FINITE(foo))
            value = R_PosInf;
        else
            value += foo;
        if (! my_standard_deviation)
            if (bee[i] != 0 && dee[i] == 0)
                value = R_PosInf;
    }

    double *lu = (double *) R_alloc(len_bee * len_bee, sizeof(double));
    int *ipiv = (int *) R_alloc(len_bee, sizeof(int));
    LU_mat_times_diag_plus_identity(my_zwz, dee, &len_bee, lu, ipiv);
    value += log_det_mat_times_diag_plus_identity(lu, &len_bee) / 2.0;

    // At this point have done value (for either parameterization)

    SEXP result, resultnames, val, grad, hess;
    PROTECT(result = allocVector(VECSXP, my_deriv + 1));
    PROTECT(resultnames = allocVector(STRSXP, my_deriv + 1));
    SET_STRING_ELT(resultnames, 0, mkChar("value"));
    if (my_deriv >= 1)
        SET_STRING_ELT(resultnames, 1, mkChar("gradient"));
    if (my_deriv >= 2)
        SET_STRING_ELT(resultnames, 2, mkChar("hessian"));
    namesgets(result, resultnames);

    // check that value is finite or +Inf
    // should have more careful analysis, but this is better than nothing
    if (! (R_finite(value) || value == R_PosInf))
        value = R_PosInf;

    PROTECT(val = ScalarReal(value));
    SET_VECTOR_ELT(result, 0, val);
    if (my_deriv == 0) {
        UNPROTECT(3);
        return result;
    }

    // At this point have returned value if that was all that was asked for
    // now only gradient and perhaps hessian left to do

    inverse_mat_times_diag_plus_identity(lu, ipiv, &len_bee);
    double *inv = lu;
    // now inv is (zwz %*% diag(dee) + diag(len_bee))^-1

    // C object gradient should have gradient for saturated aster model
    // with respect to (alpha, b)

    int len_alphabee = len_alpha + len_bee;
    int len_alphabeenu = len_alphabee + len_nu;
    PROTECT(grad = allocVector(REALSXP, len_alphabeenu));
    SET_VECTOR_ELT(result, 1, grad);
    double *my_grad = REAL(grad);

    // (14a) in design doc
    memset(my_grad, 0, len_alphabeenu * sizeof(double));
    // (15a) and (15b) in design doc
    memcpy(my_grad, gradient, len_alphabee * sizeof(double));
    double *my_grad_bee = my_grad + len_alpha;
    double *my_grad_nu = my_grad_bee + len_bee;

    // (15c) in design doc
    for (int k = 0; k < len_nu; k++)
         my_grad_nu[k] += partial_log_det_mat_times_diag_plus_identity(inv,
             dee_idx, my_zwz, &len_bee, &k) / 2.0;

    // At this point have gradient of p1 + p3 w. r. t. (alpha, b, nu)

    if (! my_standard_deviation) {

        // (17a) in design doc
        for (int i = 0; i < len_bee; i++)
            if (dee[i] > 0.0)
                my_grad_bee[i] += bee[i] / dee[i];

        // (17b) in design doc
        for (int i = 0; i < len_bee; i++) {
            int k = dee_idx[i];
            if (dee[i] > 0.0)
                my_grad_nu[k] -= (bee[i] * bee[i]) / (nu[k] * nu[k]) / 2.0;
        }

        // at this point have gradient of objective function with respect to
        // (alpha, b, nu) except for terms that are not differentiable, which
        // were skipped
        if (my_deriv == 1) {
            UNPROTECT(4);
            return result;
        }
    }

    // At this point have only Hessian to do if differentiating w. r. t
    // (alpha, b, nu), otherwise both gradient and Hessian.  Should say
    // Hessian if asked for.

    if (my_deriv == 2) {

        PROTECT(hess = allocMatrix(REALSXP, len_alphabeenu, len_alphabeenu));
        SET_VECTOR_ELT(result, 2, hess);
        double *my_hess = REAL(hess);
        // (14b) in design doc
        memset(my_hess, 0, len_alphabeenu * len_alphabeenu * sizeof(double));

#define NN len_alphabeenu
#define HESS_ALPHA_ALPHA(i, j) my_hess[(i) + NN * (j)]
#define HESS_ALPHA_BEE(i, j) my_hess[(i) + NN * (len_alpha + (j))]
#define HESS_BEE_ALPHA(i, j) my_hess[(len_alpha + (i)) + NN * (j)]
#define HESS_ALPHA_NU(i, j) my_hess[(i) + NN * (len_alphabee + (j))]
#define HESS_NU_ALPHA(i, j) my_hess[(len_alphabee + (i)) + NN * (j)]
#define HESS_BEE_BEE(i, j) my_hess[(len_alpha + (i)) + NN * (len_alpha + (j))]
#define HESS_BEE_NU(i, j) my_hess[(len_alpha + (i)) + NN * (len_alphabee + (j))]
#define HESS_NU_BEE(i, j) my_hess[(len_alphabee + (i)) + NN * (len_alpha + (j))]
#define HESS_NU_NU(i, j) my_hess[(len_alphabee + (i)) + NN * (len_alphabee + (j))]

        // (16a) and (16b) and (16c) in design doc
        for (int i = 0; i < len_alphabee; i++)
            for (int j = 0; j < len_alphabee; j++)
                my_hess[i + len_alphabeenu * j] = hessian[i + len_alphabee * j];

        // (16d) in design doc
        for (int k1 = 0; k1 < len_nu; k1++)
            for (int k2 = k1; k2 < len_nu; k2++) {
                double foo =
                    second_partial_log_det_mat_times_diag_plus_identity(
                    inv, dee_idx, my_zwz, &len_bee, &k1, &k2) / 2.0;
                HESS_NU_NU(k1, k2) -= foo;
                // force symmetry
                if (k1 != k2)
		    HESS_NU_NU(k2, k1) -= foo;
            }

        // at this point have derivatives of p1 + p3 w. r. t. (alpha, b, nu)

        if (! my_standard_deviation) {

            // (17c) in design doc
            for (int i = 0; i < len_bee; i++)
                if (dee[i] > 0.0)
                    HESS_BEE_BEE(i, i) += 1.0 / dee[i];

            // (17d) in design doc
            for (int i = 0; i < len_bee; i++) {
                int k = dee_idx[i];
                double foo = bee[i] / (nu[k] * nu[k]);
                HESS_BEE_NU(i, k) -= foo;
                HESS_NU_BEE(k, i) -= foo;
            }

            // (17e) in design doc
            for (int i = 0; i < len_bee; i++) {
                int k = dee_idx[i];
                double foo = (bee[i] * bee[i]) / (nu[k] * nu[k] * nu[k]);
                HESS_NU_NU(k, k) += foo;
            }

            // at this point have derivatives of p w. r. t. (alpha, b, nu)
            UNPROTECT(5);
            return result;
        }

        // for rest of function only doing derivatives w. r. t (alpha, c, sigma)

        // (18a) in design doc
        for (int i = 0; i < len_alpha; i++)
            for (int j = 0; j < len_bee; j++) {
                int k = dee_idx[j];
                double foo = HESS_ALPHA_BEE(i, j) * cee[j];
                HESS_ALPHA_NU(i, k) += foo;
                HESS_NU_ALPHA(k, i) += foo;
            }

        // (18b) in design doc
        for (int i = 0; i < len_alpha; i++)
            for (int j = 0; j < len_bee; j++) {
                int k = dee_idx[j];
                HESS_ALPHA_BEE(i, j) *= sigma[k];
                HESS_BEE_ALPHA(j, i) *= sigma[k];
            }

        // (18c) in design doc, middle term
        for (int k1 = 0; k1 < len_nu; k1++)
            for (int k2 = 0; k2 < len_nu; k2++)
                HESS_NU_NU(k1, k2) *= 4.0 * sigma[k1] * sigma[k2];
        // (18c) in design doc, third term
        for (int k = 0; k < len_nu; k++)
            HESS_NU_NU(k, k) += 2.0 * my_grad_nu[k];
        // (18c) in design doc, first term
        for (int i = 0; i < len_bee; i++) {
            int ki = dee_idx[i];
            for (int j = 0; j < len_bee; j++) {
                int kj = dee_idx[j];
                HESS_NU_NU(ki, kj) += HESS_BEE_BEE(i, j) * cee[i] * cee[j];
            }
        }

        // (18d) in design doc
        for (int i = 0; i < len_bee; i++) {
            int ki = dee_idx[i];
            for (int j = 0; j < len_bee; j++) {
                int kj = dee_idx[j];
                double foo = HESS_BEE_BEE(i, j) * cee[i] * sigma[kj];
                HESS_BEE_NU(j, ki) += foo;
                HESS_NU_BEE(ki, j) += foo;
            }
            double foo = my_grad_bee[i];
            HESS_BEE_NU(i, ki) += foo;
            HESS_NU_BEE(ki, i) += foo;
        }

        // (18e) in design doc
        for (int i = 0; i < len_bee; i++) {
            int ki = dee_idx[i];
            for (int j = 0; j < len_bee; j++) {
                int kj = dee_idx[j];
                HESS_BEE_BEE(i, j) *= sigma[ki] * sigma[kj];
            }
        }

        // (20b) in design doc
        for (int i = 0; i < len_bee; i++)
            HESS_BEE_BEE(i, i) += 1.0;
    }

    // Hessian done, back to gradient

    // (19a) in design doc, second term
    for (int k = 0; k < len_nu; k++)
        my_grad_nu[k] *= 2 * sigma[k];
    // (19a) in design doc, first term
    for (int j = 0; j < len_bee; j++) {
        int k = dee_idx[j];
        my_grad_nu[k] += cee[j] * my_grad_bee[j];
    }

    // (19b) in design doc
    for (int j = 0; j < len_bee; j++) {
        int k = dee_idx[j];
        my_grad_bee[j] *= sigma[k];
    }

    // (20a)
    for (int j = 0; j < len_bee; j++)
        my_grad_bee[j] += cee[j];

    UNPROTECT(3 + my_deriv);
    return result;
}

