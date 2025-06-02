
#include <math.h>
#include <stddef.h>
#include <string.h>
#include "aster.h"
#include "raster.h"

/* Bernoulli */

#if defined(__GNUC__) || defined(__clang__)
static int bernoulli_parval(double xpred,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static int bernoulli_parval(double xpred, double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* must be nonnegative integer */
    int foo = 1;
    foo = foo && (xpred == ceil(xpred));
    foo = foo && (xpred >= 0.0);
    return foo;
}

static int bernoulli_validate(double x, double xpred, double hyper1,
    double hyper2)
{
    /* xpred must be valid */
    /* x must be integer between 0 and xpred */
    int foo = 1;
    foo = foo && bernoulli_parval(xpred, hyper1, hyper2);
    foo = foo && (x == ceil(x));
    foo = foo && (x >= 0.0);
    foo = foo && (x <= xpred);
    return foo;
}

#if defined(__GNUC__) || defined(__clang__)
static int bernoulli_hypval(double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static int bernoulli_hypval(double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    return 1;
}

#if defined(__GNUC__) || defined(__clang__)
static double bernoulli(int deriv, double theta,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static double bernoulli(int deriv, double theta, double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    double foo, bar;

    switch (deriv) {
    case 0:
        if (theta <= 0)
            return my_log1p(exp(theta));
        else
            return theta + my_log1p(exp(- theta));
    case 1:
        return 1.0 / (1.0 + exp(- theta));
    case 2:
        theta = fabs(theta);
        foo = exp(- theta);
        bar = 1.0 + foo;
        return (foo / bar) / bar;
    default:
        die("deriv %d not valid", deriv);
    }
}

#if defined(__GNUC__) || defined(__clang__)
static double bernoulli_simulate(double xpred, double theta,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static double bernoulli_simulate(double xpred, double theta, double hyper1,
    double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    if (xpred == 0.0) {
        return 0.0;
    } else {
        double p = 1.0 / (1.0 + exp(- theta));
        return my_rbinom(xpred, p);
    }
}

/* Poisson */

#if defined(__GNUC__) || defined(__clang__)
static int poisson_parval(double xpred,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static int poisson_parval(double xpred, double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* must be nonnegative real */
    return (xpred >= 0.0);
}

static int poisson_validate(double x, double xpred, double hyper1,
    double hyper2)
{
    /* xpred must be valid */
    /* x must be nonnegative integer */
    /* xpred == 0 implies x == 0 */
    int foo = 1;
    foo = foo && poisson_parval(xpred, hyper1, hyper2);
    foo = foo && (x == ceil(x));
    foo = foo && (x >= 0.0);
    foo = foo && (xpred > 0.0 || x == 0.0);
    return foo;
}

#if defined(__GNUC__) || defined(__clang__)
static int poisson_hypval(double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static int poisson_hypval(double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    return 1;
}

#if defined(__GNUC__) || defined(__clang__)
static double poisson(int deriv, double theta,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static double poisson(int deriv, double theta, double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    switch (deriv) {
    case 0:
    case 1:
    case 2:
        return exp(theta);
    default:
        die("deriv %d not valid", deriv);
    }
}

#if defined(__GNUC__) || defined(__clang__)
static double poisson_simulate(double xpred, double theta,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static double poisson_simulate(double xpred, double theta, double hyper1,
    double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    double mu = xpred * exp(theta);
    if (mu == 0.0)
        return 0.0;
    return my_rpois(mu);
}

#ifdef ASTER_OLD_STUFF
/* zero-truncated Poisson (now obsolete -- replace by k-truncated) */

static int non_zero_poisson_validate(double x, double xpred,
    double hyper1, double hyper2)
{
    int foo = 1;
    foo = foo && (xpred == my_round(xpred));
    foo = foo && (x == my_round(x));
    foo = foo && (xpred >= 0.0);
    foo = foo && (x >= 0.0);
    foo = foo && (xpred > 0.0 || x == 0.0);
    foo = foo && (xpred == 0.0 || x > 0.0);
    return foo;
}

static double non_zero_poisson(int deriv, double theta, double hyper1,
    double hyper2)
{
    double mu = exp(theta);
    double tau;
    if (1.0 - mu == 1.0)
        tau = 1.0;
    else
        tau = mu / (- my_expm1(- mu));
    switch (deriv) {
    case 0:
        return mu + my_log1p(- exp(- mu));
    case 1:
        return tau;
    case 2:
        return tau * (1.0 - tau * exp(- mu));
    default:
        die("deriv %d not valid", deriv);
    }
}

static double non_zero_poisson_simulate(double xpred, double theta,
    double hyper1, double hyper2)
{
    double mu = exp(theta);
    double result = 0.0;
    int i;

    for (i = 0; i < xpred; ++i)
        result += my_rnzp(mu);
    return result;
}

/* two-truncated Poisson (now obsolete -- replace by k-truncated) */

static int two_trunc_poisson_validate(double x, double xpred,
    double hyper1, double hyper2)
{
    int foo = 1;
    foo = foo && (xpred == my_round(xpred));
    foo = foo && (x == my_round(x));
    foo = foo && (xpred >= 0.0);
    foo = foo && (x >= 0.0);
    foo = foo && (xpred > 0.0 || x == 0.0);
    foo = foo && (xpred == 0.0 || x > 2.0);
    return foo;
}

static double two_trunc_poisson(int deriv, double theta, double hyper1,
    double hyper2)
{
    int k = 2;
    double mu = exp(theta);
    double foo, bar, baz, qux, alpha;

    switch (deriv) {
    case 0:
        return mu + my_ppois(k, mu, 0, 1);
    case 1:
        foo = my_ppois(k + 1, mu, 0, 0);
        if (foo == 0.0) {
            return mu + k + 1;
        } else {
            bar = my_dpois(k + 1, mu, 0);
            return mu + (k + 1) / (1.0 + foo / bar);
        }
    case 2:
        foo = my_ppois(k + 1, mu, 0, 0);
        if (foo == 0.0) {
            qux = 0.0;
        } else {
            bar = my_dpois(k + 1, mu, 0);
            qux = foo / bar;
        }
        alpha = (k + 1) / (1.0 + qux);
        if (qux < 1.0) {
            baz = qux / (1.0 + qux);
        } else {
            baz = 1.0 / (1.0 / qux + 1.0);
        }
        return mu * (1.0 - alpha * (1.0 - (k + 1) * baz / mu));
    default:
        die("deriv %d not valid", deriv);
    }
}

static double two_trunc_poisson_simulate(double xpred, double theta,
    double hyper1, double hyper2)
{
    double mu = exp(theta);
    double result = 0.0;
    int i;

    for (i = 0; i < xpred; ++i)
        result += my_rktp(2, mu);
    return result;
}
#endif /* ASTER_OLD_STUFF */

/* k-truncated Poisson */

#if defined(__GNUC__) || defined(__clang__)
static int trunc_poisson_parval(double xpred,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static int trunc_poisson_parval(double xpred, double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* must be nonnegative integer */
    int foo = 1;
    foo = foo && (xpred == ceil(xpred));
    foo = foo && (xpred >= 0.0);
    return foo;
}

static int trunc_poisson_validate(double x, double xpred,
    double hyper1, double hyper2)
{
    double k = hyper1;
    /* xpred must be valid */
    /* x must be integer */
    /* xpred == 0 implies x == 0 */
    /* xpred > 0 implies x >= xpred * (k + 1) */
    int foo = 1;
    foo = foo && trunc_poisson_parval(xpred, hyper1, hyper2);
    foo = foo && (x == ceil(x));
    foo = foo && (xpred > 0.0 || x == 0.0);
    foo = foo && (xpred == 0.0 || x >= xpred * (k + 1));
    return foo;
}

#if defined(__GNUC__) || defined(__clang__)
static int trunc_poisson_hypval(double hyper1,
    double hyper2 __attribute__ ((unused)))
#else
static int trunc_poisson_hypval(double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* hyper1 (truncation) must be nonnegative integer */
    int foo = 1;
    foo = foo && (hyper1 == ceil(hyper1));
    foo = foo && (hyper1 >= 0.0);
    return foo;
}

#if defined(__GNUC__) || defined(__clang__)
static double trunc_poisson(int deriv, double theta, double hyper1,
    double hyper2 __attribute__ ((unused)))
#else
static double trunc_poisson(int deriv, double theta,
    double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    int k = hyper1;
    double mu = exp(theta);
    double foo, bar, baz, qux, alpha;

    /* special case k == 0 backporting stuff from aster2 */
    /* also uses analysis in Stat 3701 notes and in the */
    /* devel/arithmetic/zero-truncated-poisson directory */
    /* of this package */
    if (k == 0) {
        if (deriv > 2 || deriv < 0)
            die("deriv %d not valid", deriv);
        double m = exp(theta);
        double mu;
        if (theta > (- 4.0)) {
            if (deriv == 0)
                return m + my_log1p(- exp(- m));
            mu = m / (- expm1(- m));
            if (deriv == 1)
                return mu;
            /* deriv == 2 */
            if (isinf(mu)) {
                return mu;
            } else {
                return mu * (1.0 - mu * exp(- m));
            }
        } else /* theta <= -4 */ {
            double bar = m / 2.0 * (1.0 + m / 3.0 * (1.0 + m / 4.0 *
                (1.0 + m / 5.0 * (1.0 + m / 6.0 * (1.0 + m / 7.0 *
                (1.0 + m / 8.0))))));
            if (deriv == 0)
                return theta + my_log1p(bar);
            mu = m + 1.0 / (1.0 + bar);
            if (deriv == 1)
                return mu;
            /* deriv == 2 */
            double msq = m * m;
            return m / 2 * (1 + m / 3 * (1 - msq / 30 *
                (1 - msq / 28 * (1 - msq / 30))));
        }
    }

    switch (deriv) {
    case 0:
        return mu + my_ppois(k, mu, 0, 1);
    case 1:
        foo = my_ppois(k + 1, mu, 0, 0);
        if (foo == 0.0) {
            return mu + k + 1;
        } else {
            bar = my_dpois(k + 1, mu, 0);
            return mu + (k + 1) / (1.0 + foo / bar);
        }
    case 2:
        foo = my_ppois(k + 1, mu, 0, 0);
        if (foo == 0.0) {
            qux = 0.0;
        } else {
            bar = my_dpois(k + 1, mu, 0);
            qux = foo / bar;
        }
        alpha = (k + 1) / (1.0 + qux);
        if (qux < 1.0) {
            baz = qux / (1.0 + qux);
        } else {
            baz = 1.0 / (1.0 / qux + 1.0);
        }
        return mu * (1.0 - alpha * (1.0 - (k + 1) * baz / mu));
    default:
        die("deriv %d not valid", deriv);
    }
}

#if defined(__GNUC__) || defined(__clang__)
static double trunc_poisson_simulate(double xpred, double theta,
    double hyper1, double hyper2 __attribute__ ((unused)))
#else
static double trunc_poisson_simulate(double xpred, double theta,
    double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    int k = hyper1;
    double mu = exp(theta);
    double result = 0.0;
    int i;

    for (i = 0; i < xpred; ++i)
        result += my_rktp(k, mu);
    return result;
}

/* (untruncated) negative binomial */

#if defined(__GNUC__) || defined(__clang__)
static int neg_bin_parval(double xpred,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static int neg_bin_parval(double xpred, double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* must be nonnegative real */
    return (xpred >= 0.0);
}

static int neg_bin_validate(double x, double xpred,
    double hyper1, double hyper2)
{
    /* xpred must be valid */
    /* x must be integer */
    /* xpred == 0 implies x == 0 */
    int foo = 1;
    foo = foo && neg_bin_parval(xpred, hyper1, hyper2);
    foo = foo && (x == ceil(x));
    foo = foo && (xpred > 0.0 || x == 0.0);
    return foo;
}

#if defined(__GNUC__) || defined(__clang__)
static int neg_bin_hypval(double hyper1,
    double hyper2 __attribute__ ((unused)))
#else
static int neg_bin_hypval(double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* hyper1 (size) must be nonnegative real */
    return (hyper1 >= 0.0);
}

#if defined(__GNUC__) || defined(__clang__)
static double neg_bin(int deriv, double theta, double hyper1,
    double hyper2 __attribute__ ((unused)))
#else
static double neg_bin(int deriv, double theta,
    double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    double alpha = hyper1;
    double p, q, foo;

    if (theta >= 0.0) {
        switch (deriv) {
        case 0:
            return my_posinf();
        case 1:
        case 2:
            return my_nan();
        default:
            die("deriv %d not valid", deriv);
        }
    }

    switch (deriv) {
    case 0:
        return alpha * my_log1p(1.0 / my_expm1(- theta));
    case 1:
    case 2:
        q = exp(theta);
        p = (- my_expm1(theta));
        foo = alpha * q / p;
        if (deriv == 1)
            return foo;
        else
            return foo / p;
    default:
        die("deriv %d not valid", deriv);
    }
}

#if defined(__GNUC__) || defined(__clang__)
static double neg_bin_simulate(double xpred, double theta,
    double hyper1, double hyper2 __attribute__ ((unused)))
#else
static double neg_bin_simulate(double xpred, double theta,
    double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    double alpha = hyper1;
    double p = - my_expm1(theta);
    return my_rnbinom(xpred * alpha, p);
}

/* truncated negative binomial */

#if defined(__GNUC__) || defined(__clang__)
static int trunc_neg_bin_parval(double xpred,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static int trunc_neg_bin_parval(double xpred, double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* must be nonnegative integer */
    int foo = 1;
    foo = foo && (xpred == ceil(xpred));
    foo = foo && (xpred >= 0.0);
    return foo;
}

static int trunc_neg_bin_validate(double x, double xpred,
    double hyper1, double hyper2)
{
    double k = hyper2;
    /* xpred must be valid */
    /* x must be integer */
    /* xpred == 0 implies x == 0 */
    /* xpred > 0 implies x >= xpred * (k + 1) */
    int foo = 1;
    foo = foo && trunc_neg_bin_parval(xpred, hyper1, hyper2);
    foo = foo && (x == ceil(x));
    foo = foo && (xpred > 0.0 || x == 0.0);
    foo = foo && (xpred == 0.0 || x >= xpred * (k + 1));
    return foo;
}

static int trunc_neg_bin_hypval(double hyper1, double hyper2)
{
    /* hyper1 (size) must be nonnegative real */
    /* hyper2 (truncation) must be nonnegative integer */
    int foo = 1;
    foo = foo && (hyper1 >= 0.0);
    foo = foo && (hyper2 == ceil(hyper2));
    foo = foo && (hyper2 >= 0.0);
    return foo;
}

static double trunc_neg_bin(int deriv, double theta,
    double hyper1, double hyper2)
{
    double alpha = hyper1;
    double k = hyper2;
    double p, q, foo, mu, numer, beta, beep;

    if (theta >= 0.0) {
        switch (deriv) {
        case 0:
            return my_posinf();
        case 1:
        case 2:
            return my_nan();
        default:
            die("deriv %d not valid", deriv);
        }
    }

    switch (deriv) {
    case 0:
        /* psi(theta) = - log(1 - exp(theta)) + log Pr(Y > k) */
        p = - my_expm1(theta);
        return alpha * my_log1p(1.0 / my_expm1(- theta)) +
            my_pnbinom(k, alpha, p, 0, 1);
    case 1:
    case 2:
        q = exp(theta);
        p = - my_expm1(theta);
        mu = alpha * q / p;
        numer = my_pnbinom(k + 1, alpha, p, 0, 0);
        if (numer == 0.0)
            beta = 0.0;
        else
            beta = numer / my_dnbinom(k + 1, alpha, p, 0);
        if (deriv == 1)
            return mu + (k + 1) / (1.0 + beta) / p;
        if (beta < 1.0)
            beep = beta / (1.0 + beta);
        else
            beep = 1.0 / (1.0 / beta + 1.0);
        foo = k + 1.0 + alpha;
        foo = (- q) + foo * q / (1.0 + beta) + (alpha - p * foo) * beep;
        foo *= (k + 1.0) / (1.0 + beta) / p;
        foo = mu - foo;
        foo /= p;
        return foo;
    default:
        die("deriv %d not valid", deriv);
    }
}

static double trunc_neg_bin_simulate(double xpred, double theta,
    double hyper1, double hyper2)
{
    double alpha = hyper1;
    double k = hyper2;
    double p = (- my_expm1(theta));
    double q = exp(theta);
    double mu = alpha * q / p;
    double result = 0.0;
    int i;

    for (i = 0; i < xpred; ++i)
        result += my_rktnb(alpha, k, mu);
    return result;
}

/* normal location (known scale) */

#if defined(__GNUC__) || defined(__clang__)
static int norm_loc_parval(double xpred,
    double hyper1 __attribute__ ((unused)),
    double hyper2 __attribute__ ((unused)))
#else
static int norm_loc_parval(double xpred, double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* must be nonnegative real */
    return (xpred >= 0.0);
}

static int norm_loc_validate(double x, double xpred,
    double hyper1, double hyper2)
{
    /* xpred must be valid */
    /* x must be real */
    /* xpred == 0 implies x == 0 */
    int foo = 1;
    foo = foo && norm_loc_parval(xpred, hyper1, hyper2);
    foo = foo && (xpred > 0.0 || x == 0.0);
    return foo;
}

#if defined(__GNUC__) || defined(__clang__)
static int norm_loc_hypval(double hyper1,
    double hyper2 __attribute__ ((unused)))
#else
static int norm_loc_hypval(double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    /* hyper1 (size) must be positive real */
    return (hyper1 > 0.0);
}

#if defined(__GNUC__) || defined(__clang__)
static double norm_loc(int deriv, double theta, double hyper1,
    double hyper2 __attribute__ ((unused)))
#else
static double norm_loc(int deriv, double theta,
    double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    double sigmasq = hyper1 * hyper1;

    switch (deriv) {
    case 0:
        return sigmasq * theta * theta / 2.0;
    case 1:
        return sigmasq * theta;
    case 2:
        return sigmasq;
    default:
        die("deriv %d not valid", deriv);
    }
}

#if defined(__GNUC__) || defined(__clang__)
static double norm_loc_simulate(double xpred, double theta,
    double hyper1, double hyper2 __attribute__ ((unused)))
#else
static double norm_loc_simulate(double xpred, double theta,
    double hyper1, double hyper2)
#endif /* defined(__GNUC__) || defined(__clang__) */
{
    double mu = theta * hyper1 * hyper1 * xpred;
    double sigma = hyper1 * sqrt(xpred);
    return my_rnorm(mu, sigma);
}

typedef double (*famfun_ptr)(int deriv, double theta, double hyper1,
    double hyper2);
typedef int (*famval_ptr)(double x, double xpred, double hyper1,
    double hyper2);
typedef int (*fampval_ptr)(double xpred, double hyper1, double hyper2);
typedef int (*famhval_ptr)(double hyper1, double hyper2);
typedef double (*famsim_ptr)(double xpred, double theta, double hyper1,
    double hyper2);

struct superfamtab {
    char * const name;
    famfun_ptr const psi;
    famval_ptr const validate;
    fampval_ptr const validparent;
    famhval_ptr const validhyper;
    famsim_ptr const simulate;
    int const mincard;
    int const maxcard;
    int const nhyper;
    char * const name_hyper1;
    char * const name_hyper2;
    double const origin;
};

static struct superfamtab const mysuperfamtab[] =
{
    {"bernoulli", bernoulli, bernoulli_validate, bernoulli_parval,
        bernoulli_hypval, bernoulli_simulate, 1, 1, 0, "", "", 0.0},
    {"poisson", poisson, poisson_validate, poisson_parval,
        poisson_hypval, poisson_simulate, 1, 1, 0, "", "", 0.0},
    {"truncated.poisson", trunc_poisson, trunc_poisson_validate,
        trunc_poisson_parval, trunc_poisson_hypval, trunc_poisson_simulate,
        1, 1, 1, "truncation", "", 0.0},
    {"negative.binomial", neg_bin, neg_bin_validate,
        neg_bin_parval, neg_bin_hypval, neg_bin_simulate,
        1, 1, 1, "size", "", -1.0},
    {"truncated.negative.binomial", trunc_neg_bin, trunc_neg_bin_validate,
        trunc_neg_bin_parval, trunc_neg_bin_hypval, trunc_neg_bin_simulate,
        1, 1, 2, "size", "truncation", -1.0},
    {"normal.location", norm_loc, norm_loc_validate,
        norm_loc_parval, norm_loc_hypval, norm_loc_simulate,
        1, 1, 1, "sd", "", 0.0},
    {NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0.0},
};

struct famtab {
    char *name;
    famfun_ptr psi;
    famval_ptr validate;
    fampval_ptr validparent;
    famsim_ptr simulate;
    int mincard;
    int maxcard;
    int nhyper;
    double hyper1;
    double hyper2;
    char *name_hyper1;
    char *name_hyper2;
    double origin;
};

static struct famtab myfamtab[MAXNFAM];
static int nfam = 0;

void aster_clear_families(void)
{
    nfam = 0;
}

void aster_add_family(char **name, double *hyper, int *nhyper)
{
    int ifam;
    double hyper1 = 0.0;
    double hyper2 = 0.0;

    if (nfam == MAXNFAM)
        die("number of families exceeds family table size");
    for (ifam = 0; ; ++ifam) {
        if (mysuperfamtab[ifam].name == NULL)
            die("family \"%s\" not found", name[0]);
        if (strcmp(mysuperfamtab[ifam].name, name[0]) == 0)
            break;
    }
    if (mysuperfamtab[ifam].nhyper != nhyper[0])
        die("family \"%s\" has %d hyperparameters, %d specified",
            name[0], mysuperfamtab[ifam].nhyper, nhyper[0]);
    if (nhyper[0] >= 1)
        hyper1 = hyper[0];
    if (nhyper[0] >= 2)
        hyper2 = hyper[1];
    if (! mysuperfamtab[ifam].validhyper(hyper1, hyper2))
        die("family \"%s\" specified with invalid hyperparameters", name[0]);

    myfamtab[nfam].name = mysuperfamtab[ifam].name;
    myfamtab[nfam].psi = mysuperfamtab[ifam].psi;
    myfamtab[nfam].validate = mysuperfamtab[ifam].validate;
    myfamtab[nfam].validparent = mysuperfamtab[ifam].validparent;
    myfamtab[nfam].simulate = mysuperfamtab[ifam].simulate;
    myfamtab[nfam].mincard = mysuperfamtab[ifam].mincard;
    myfamtab[nfam].maxcard = mysuperfamtab[ifam].maxcard;
    myfamtab[nfam].nhyper = mysuperfamtab[ifam].nhyper;
    myfamtab[nfam].hyper1 = hyper1;
    myfamtab[nfam].hyper2 = hyper2;
    myfamtab[nfam].name_hyper1 = mysuperfamtab[ifam].name_hyper1;
    myfamtab[nfam].name_hyper2 = mysuperfamtab[ifam].name_hyper2;
    myfamtab[nfam].origin = mysuperfamtab[ifam].origin;
    ++nfam;
}

void aster_get_family(int *famin, char **name, double *hyper, int *nhyper,
    char **hypername, double *origin)
{
    int fam = famin[0];
    if (fam <= 0 || fam > nfam) {
        name[0] = "";
        return;
    }

    int ifam = fam - 1; /* fam is 1-origin index */

    name[0] = myfamtab[ifam].name;
    nhyper[0] = myfamtab[ifam].nhyper;
    origin[0] = myfamtab[ifam].origin;
    if (nhyper[0] >= 1) {
        hyper[0] = myfamtab[ifam].hyper1;
        hypername[0] = myfamtab[ifam].name_hyper1;
    }
    if (nhyper[0] >= 2) {
        hyper[1] = myfamtab[ifam].hyper2;
        hypername[1] = myfamtab[ifam].name_hyper2;
    }
    return;
}

void aster_get_superfamily(int *famin, char **name, int *nhyper,
    char **hypername)
{
    int fam = famin[0];
    if (fam <= 0) {
        name[0] = "";
        return;
    }

    int ifam = fam - 1; /* fam is 1-origin index */

    for (int i = 0; i <= ifam; ++i)
        if (mysuperfamtab[ifam].name == NULL) {
           name[0] = "";
           return;
        }

    name[0] = mysuperfamtab[ifam].name;
    nhyper[0] = mysuperfamtab[ifam].nhyper;
    if (nhyper[0] >= 1)
        hypername[0] = mysuperfamtab[ifam].name_hyper1;
    if (nhyper[0] >= 2)
        hypername[1] = mysuperfamtab[ifam].name_hyper2;
    return;
}

void aster_byname_superfamily(char **name, int *nhyper, char **hypername)
{
    int ifam;

    for (ifam = 0; ; ++ifam) {
        if (mysuperfamtab[ifam].name == NULL)
            die("family \"%s\" not found", name[0]);
        if (strcmp(mysuperfamtab[ifam].name, name[0]) == 0)
            break;
    }

    nhyper[0] = mysuperfamtab[ifam].nhyper;
    if (nhyper[0] >= 1)
        hypername[0] = mysuperfamtab[ifam].name_hyper1;
    if (nhyper[0] >= 2)
        hypername[1] = mysuperfamtab[ifam].name_hyper2;
    return;
}

void aster_family(int *famin, int *derivin, double *thetain, double *value)
{
    int fam = famin[0];
    int deriv = derivin[0];
    double theta = thetain[0];

    if (fam <= 0 || fam > nfam)
         die("family %d not valid", fam);

    int ifam = fam - 1; /* fam is 1-origin index */
    double hyper1 = myfamtab[ifam].hyper1;
    double hyper2 = myfamtab[ifam].hyper2;
    value[0] = myfamtab[ifam].psi(deriv, theta, hyper1, hyper2);
}

int aster_family_validate(int fam, double x, double xpred)
{
    if (fam <= 0 || fam > nfam)
         die("family %d not valid", fam);

    int ifam = fam - 1; /* fam is 1-origin index */
    double hyper1 = myfamtab[ifam].hyper1;
    double hyper2 = myfamtab[ifam].hyper2;

    return myfamtab[ifam].validate(x, xpred, hyper1, hyper2);
}

int aster_family_validate_parent(int fam, double xpred)
{
    if (fam <= 0 || fam > nfam)
         die("family %d not valid", fam);

    int ifam = fam - 1; /* fam is 1-origin index */
    double hyper1 = myfamtab[ifam].hyper1;
    double hyper2 = myfamtab[ifam].hyper2;

    return myfamtab[ifam].validparent(xpred, hyper1, hyper2);
}

int aster_family_number_validate(int fam)
{
    return fam > 0 && fam <= nfam;
}

double aster_family_origin(int fam)
{
    if (fam <= 0 || fam > nfam)
         die("family %d not valid", fam);

    int ifam = fam - 1; /* fam is 1-origin index */
    return myfamtab[ifam].origin;
}

double aster_family_simulate(int fam, double xpred, double theta)
{
    if (fam <= 0 || fam > nfam)
         die("family %d not valid", fam);

    int ifam = fam - 1; /* fam is 1-origin index */
    double hyper1 = myfamtab[ifam].hyper1;
    double hyper2 = myfamtab[ifam].hyper2;

    return myfamtab[ifam].simulate(xpred, theta, hyper1, hyper2);
}

#ifdef ASTER_OLD_STUFF
struct funtab {
    char *name;
    famfun_ptr fun;
    famval_ptr validate;
    famsim_ptr simulate;
};

static struct funtab myfuntab[] =
{
    {"bernoulli", bernoulli, bernoulli_validate, bernoulli_simulate},
    {"poisson", poisson, poisson_validate, poisson_simulate},
    {"non.zero.poisson", non_zero_poisson, non_zero_poisson_validate,
        non_zero_poisson_simulate},
    {"two.trunc.poisson", two_trunc_poisson, two_trunc_poisson_validate,
        two_trunc_poisson_simulate},
    {NULL, NULL, NULL, NULL},
};

void aster_family(int *famin, int *derivin, double *thetain, double *value)
{
    int fam = famin[0];
    int deriv = derivin[0];
    double theta = thetain[0];

    int i;
    for (i = 0; myfuntab[i].fun != NULL; ++i)
        if (i == (fam - 1)) /* fam is 1-origin index */ {
            value[0] = myfuntab[i].fun(deriv, theta, 0.0, 0.0);
            return;
        }
    die("family %d not valid", fam);
}

int aster_family_validate(int fam, double x, double xpred)
{
    int i;
    for (i = 0; myfuntab[i].validate != NULL; ++i)
        if (i == (fam - 1)) /* fam is 1-origin index */
            return myfuntab[i].validate(x, xpred, 0.0, 0.0);
    die("family %d not valid", fam);
}

char *aster_family_name(int fam)
{
    int i;
    for (i = 0; myfuntab[i].name != NULL; ++i)
        if (i == (fam - 1)) /* fam is 1-origin index */
            return myfuntab[i].name;
    return NULL;
}

double aster_family_simulate(int fam, double xpred, double theta)
{
    int i;
    for (i = 0; myfuntab[i].simulate != NULL; ++i)
        if (i == (fam - 1)) /* fam is 1-origin index */
            return myfuntab[i].simulate(xpred, theta, 0.0, 0.0);
    die("family %d not valid", fam);
}
#endif /* ASTER_OLD_STUFF */

