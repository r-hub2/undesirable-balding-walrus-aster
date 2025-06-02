
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include "aster.h"
#include "raster.h"
#include "rraster.h"

#define BUFSIZE 1024

void die(const char *format, ...)
{
    char buf[BUFSIZE];
    va_list arg;

    va_start(arg, format);
    vsnprintf(buf, BUFSIZE, format, arg);
    va_end(arg);
    buf[BUFSIZE - 1] = '\0';
    error("%s", buf);
}

void my_warn(const char *format, ...)
{
    char buf[BUFSIZE];
    va_list arg;

    va_start(arg, format);
    vsnprintf(buf, BUFSIZE, format, arg);
    va_end(arg);
    buf[BUFSIZE - 1] = '\0';
    warning("%s", buf);
}

#ifdef ASTER_OLD_STUFF
SEXP aster_families(void)
{
    SEXP result;
    int nfam = 0;
    int i;

    for (i = 0; ; ++i) {
        /* family index is 1-origin */
        char *name = aster_family_name(i + 1);
        if (name == NULL)
            break;
        else
            ++nfam;
    }

    PROTECT(result = allocVector(STRSXP, nfam));
    for (i = 0; i < nfam; ++i)
        /* family index is 1-origin */
        SET_STRING_ELT(result, i, mkChar(aster_family_name(i + 1)));
    UNPROTECT(1);
    return result;
}
#endif /* ASTER_OLD_STUFF */

double my_expm1(double x)
{
    return expm1(x);
}

double my_log1p(double x)
{
    return log1p(x);
}

#ifdef ASTER_OLD_STUFF
double my_round(double x)
{
    return rint(x);
}
#endif /* ASTER_OLD_STUFF */

double my_rbinom(double n, double p)
{
    return rbinom(n, p);
}

double my_rpois(double mu)
{
    return rpois(mu);
}

double my_ppois(double x, double lambda, int lower_tail, int log_p)
{
    return ppois(x, lambda, lower_tail, log_p);
}

double my_dpois(double x, double lambda, int give_log)
{
    return dpois(x, lambda, give_log);
}

double my_rnbinom(double n /* size */, double p /* prob */)
{
    return rnbinom(n, p);
}

double my_pnbinom(double x, double n, double p, int lower_tail, int log_p)
{
    return pnbinom(x, n, p, lower_tail, log_p);
}

double my_dnbinom(double x, double n, double p, int give_log)
{
    return dnbinom(x, n, p, give_log);
}

double my_rnorm(double mu, double sigma)
{
    return rnorm(mu, sigma);
}

double my_nan(void)
{
    return R_NaN;
}

double my_is_finite(double foo)
{
    return R_finite(foo);
}

double my_is_na_or_nan(double foo)
{
    return R_IsNA(foo) || R_IsNaN(foo);
}

double my_posinf(void)
{
    return R_PosInf;
}

double my_neginf(void)
{
    return R_NegInf;
}

#ifdef ASTER_OLD_STUFF
double my_rnzp(double mu)
{
    if (mu <= 0.0)
        die("non-positive mu in non-zero-poisson simulator\n");
    if (mu >= 1.0) {
        /* case 1 in sim.Rnw write-up */
        for (;;) {
            double x = rpois(mu);
            if (x >= 1.0)
                return x;
        }
    } else /* mu < 1.0 */ {
        /* case 2 in sim.Rnw write-up */
        double tau;
        if (1.0 - mu == 1.0)
            return 1.0;
        tau = mu / (- my_expm1(- mu));
        for (;;) {
            double x = rpois(mu) + 1.0;
            double u = unif_rand();
            if (u < 1.0 / x)
                return x;
        }
    }
}

void aster_rnzp(int *nin, int *len_xpred_in, int *len_mu_in,
    double *xpred, double *mu, double *result)
{
    int n = nin[0];
    int len_xpred = len_xpred_in[0];
    int len_mu = len_mu_in[0];
    int i, j;

    GetRNGstate();
    for (i = 0; i < n; ++i) {
        double foo = xpred[i % len_xpred];
        double bar = 0.0;
        double qux = mu[i % len_mu];
        for (j = 0; j < foo; ++j)
            bar += my_rnzp(qux);
        result[i] = bar;
    }
    PutRNGstate();
}
#endif /* ASTER_OLD_STUFF */

double my_rktp(int k, double mu)
{
    int m;
    double mdoub;

    if (mu <= 0.0)
        die("non-positive mu in k-truncated-poisson simulator\n");
    if (k < 0)
        die("negative k in k-truncated-poisson simulator\n");

    mdoub = k + 1 - mu;
    if (mdoub < 0.0)
        mdoub = 0.0;
    m = mdoub;
    if (m < mdoub)
        m = m + 1;
    /* since mu > 0.0 we have 0.0 <= mdoub < k + 1 hence 0 <= m <= k + 1 */

    for (;;) {
        double x = rpois(mu) + m;
        if (m > 0) {
            double a = 1.0;
            int j;
            double u = unif_rand();
            for (j = 0; j < m; ++j)
                a *= (k + 1 - j) / (x - j);
            if (u < a && x > k)
                return x;
        } else {
            if (x > k)
                return x;
        }
    }
}

void aster_rktp(int *nin, int *len_xpred_in, int *len_mu_in, int *len_k_in,
    double *xpred_in, double *mu_in, int *k_in, double *result)
{
    int n = nin[0];
    int len_xpred = len_xpred_in[0];
    int len_mu = len_mu_in[0];
    int len_k = len_k_in[0];
    int i, j;

    GetRNGstate();
    for (i = 0; i < n; ++i) {
        double xpred = xpred_in[i % len_xpred];
        double mu = mu_in[i % len_mu];
        int k = k_in[i % len_k];
        double foo = 0.0;
        for (j = 0; j < xpred; ++j)
            foo += my_rktp(k, mu);
        result[i] = foo;
    }
    PutRNGstate();
}

double my_rktnb(double alpha, int k, double mu)
{
    int m;
    double mdoub;
    double p = alpha / (mu + alpha);
    double q = mu / (mu + alpha);

    if (alpha <= 0.0)
        die("non-positive size in k-truncated-neg-bin simulator\n");
    if (mu <= 0.0)
        die("non-positive mu in k-truncated-neg-bin simulator\n");
    if (k < 0)
        die("negative k in k-truncated-neg-bin simulator\n");

    mdoub = (k + 1.0) * p - alpha * q;
    if (mdoub < 0.0)
        mdoub = 0.0;
    m = mdoub;
    if (m < mdoub)
        m = m + 1;
    /* since p < 1.0 and q > 0.0 we have 0.0 <= mdoub < k + 1
       hence 0 <= m <= k + 1 */

    for (;;) {
        double x = rnbinom(alpha + m, p) + m;
        if (m > 0) {
            double a = 1.0;
            int j;
            double u = unif_rand();
            for (j = 0; j < m; ++j)
                a *= (k + 1 - j) / (x - j);
            if (u < a && x > k)
                return x;
        } else {
            if (x > k)
                return x;
        }
    }
}

void aster_rktnb(int *nin, int *len_xpred_in, int *len_mu_in, int *len_k_in,
    int *len_alpha_in, double *xpred_in, double *mu_in, int *k_in,
    double *alpha_in, double *result)
{
    int n = nin[0];
    int len_xpred = len_xpred_in[0];
    int len_mu = len_mu_in[0];
    int len_k = len_k_in[0];
    int len_alpha = len_alpha_in[0];
    int i, j;

    GetRNGstate();
    for (i = 0; i < n; ++i) {
        double xpred = xpred_in[i % len_xpred];
        double mu = mu_in[i % len_mu];
        int k = k_in[i % len_k];
        double alpha = alpha_in[i % len_alpha];
        double foo = 0.0;
        for (j = 0; j < xpred; ++j)
            foo += my_rktnb(alpha, k, mu);
        result[i] = foo;
    }
    PutRNGstate();
}

void my_GetRNGstate(void)
{
    GetRNGstate();
}

void my_PutRNGstate(void)
{
    PutRNGstate();
}

void *my_malloc(size_t size)
{
    void *foo = malloc(size);
    if (foo == NULL)
        die("malloc returned null\n");
    return foo;
}

void my_free(void *ptr)
{
    free(ptr);
}

