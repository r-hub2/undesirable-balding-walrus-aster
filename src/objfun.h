#ifndef ASTER_OBJFUN_H
#define ASTER_OBJFUN_H

SEXP objfun(SEXP theta, SEXP modmat, SEXP nfixed, SEXP nrand,
    SEXP response, SEXP pred, SEXP fam, SEXP root, SEXP zwz,
    SEXP offset, SEXP standard_deviation, SEXP deriv);

#endif /* ASTER_OBJFUN_H */

