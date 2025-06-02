
#ifndef ASTER_R_ASTER_H
#define ASTER_R_ASTER_H

#include <stdarg.h>
#include <stdlib.h>

#ifndef __GNUC__
void die(const char *format, ...);
#else
void die(const char *format, ...) __attribute__ ((__noreturn__));
#endif /* __GNUC__ */

void my_warn(const char *format, ...);

double my_expm1(double x);

double my_log1p(double x);

#ifdef ASTER_OLD_STUFF
double my_round(double x);
#endif /* ASTER_OLD_STUFF */

double my_rbinom(double n, double p);

double my_rpois(double mu);

double my_ppois(double x, double lambda, int lower_tail, int log_p);

double my_dpois(double x, double lambda, int give_log);

double my_rnbinom(double n /* size */, double p /* prob */);

double my_pnbinom(double x, double n, double p, int lower_tail, int log_p);

double my_dnbinom(double x, double n, double p, int give_log);

double my_rnorm(double mu, double sigma);

double my_nan(void);

double my_is_finite(double foo);

double my_is_na_or_nan(double foo);

double my_posinf(void);

double my_neginf(void);

#ifdef ASTER_OLD_STUFF
double my_rnzp(double mu);
#endif /* ASTER_OLD_STUFF */

double my_rktp(int k, double mu);

double my_rktnb(double alpha, int k, double mu);

void my_GetRNGstate(void);

void my_PutRNGstate(void);

void *my_malloc(size_t size);

void my_free(void *ptr);

#endif /* ASTER_R_ASTER_H */

