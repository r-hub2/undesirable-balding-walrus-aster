\name{newpickle}
\alias{newpickle}
\title{Penalized Quasi-Likelihood for Aster Models}
\usage{
newpickle(alphaceesigma, fixed, random, obj, y, origin, zwz, deriv = 0)
}
\note{
  Not intended for use by naive users.  Use \code{\link{reaster}}.
  Actually no longer used by other functions in this package.
}
\description{
Evaluates the objective function for approximate maximum likelihood for
an aster model with random effects.  Uses Laplace approximation to
integrate out the random effects analytically.
The \dQuote{quasi} in the title is a misnomer in the context of aster
models but the acronym PQL for this procedure is well-established in
the generalized linear mixed models literature.
}
\arguments{
  \item{alphaceesigma}{the parameter value where the function is evaluated,
    a numeric vector, see details.}
  \item{fixed}{the model matrix for fixed effects.  The number of rows
    is \code{nrow(obj$data)}.  The number of columns is the number of fixed
    effects.}
  \item{random}{the model matrix or matrices for random effects.
    The number of rows is \code{nrow(obj$data)}.  The number of columns
    is the number of random effects in a group.  Either a matrix or a list
    each element of which is a matrix.}
  \item{obj}{aster model object, the result of a call to \code{\link{aster}}.}
  \item{y}{response vector.  May be omitted, in which case \code{obj$x}
    is used.  If supplied, must be a matrix of the same dimensions as
    \code{obj$x}.}
  \item{origin}{origin of aster model.  May be omitted, in which case
    default origin (see \code{\link{aster}}) is used.  If supplied, must be
    a matrix of the same dimensions \code{obj$x}.}
  \item{zwz}{A possible value of \eqn{Z^T W Z}{t(Z) W Z}, where \eqn{Z} is the
    model matrix for all random effects and \eqn{W} is the variance matrix of
    the response.  May be missing, in which case it is calculated from
    \code{alphaceesigma}.  See details.}
  \item{deriv}{Number of derivatives wanted, either zero or one.
    Must be zero if \code{zwz} is missing.}
}
\details{
Define
\deqn{p(\alpha, c, \sigma) = m(a + M \alpha + Z A c) + c^T c / 2 + \log \det[A Z^T W(a + M \alpha + Z A c) Z A + I]}{p(alpha, c, sigma) = m(a + M alpha + Z A c) + t(c) c / 2 + log det[A t(Z) W(a + M alpha + Z A c) Z A + I]}
where \eqn{m} is minus the log likelihood function of a saturated aster model,
where \eqn{W} is the Hessian matrix of \eqn{m},
where \eqn{a} is a known vector (the \emph{offset vector} in the terminology
of \code{\link{glm}} but the \emph{origin} in the terminology
of \code{\link{aster}}), where \eqn{M} is a known matrix, the model matrix for
fixed effects (the argument \code{fixed} of this function),
\eqn{Z} is a known matrix, the model matrix for random effects
(either the argument \code{random} of this functions if it is a matrix or
\code{Reduce(cbind, random)} if \code{random} is a list of matrices),
where \eqn{A} is a diagonal matrix whose diagonal is the vector
\code{rep(sigma, times = nrand)}
where \code{nrand} is \code{sapply(random, ncol)}
when \code{random} is a list of
matrices and \code{ncol(random)} when \code{random} is a matrix,
and where \eqn{I} is the identity matrix.
This function evaluates \eqn{p(\alpha, c, \sigma)}{p(alpha, c, sigma)}
when \code{zwz} is missing.
Otherwise it evaluates the same thing except that
\deqn{Z^T W(a + M \alpha + Z A c) Z}{t(Z) W(a + M alpha + Z A c) Z}
is replaced by \code{zwz}.
Note that \eqn{A} is a function of \eqn{\sigma}{sigma} although the
notation does not explicitly indicate this.
}
\value{
  a list with components \code{value} and \code{gradient},
  the latter missing if \code{deriv == 0}.
}
\examples{
data(radish)

pred <- c(0,1,2)
fam <- c(1,3,2)

### need object of type aster to supply to penmlogl and pickle

aout <- aster(resp ~ varb + fit : (Site * Region + Block + Pop),
    pred, fam, varb, id, root, data = radish)

### model matrices for fixed and random effects

modmat.fix <- model.matrix(resp ~ varb + fit : (Site * Region),
    data = radish)
modmat.blk <- model.matrix(resp ~ 0 + fit:Block, data = radish)
modmat.pop <- model.matrix(resp ~ 0 + fit:Pop, data = radish)

rownames(modmat.fix) <- NULL
rownames(modmat.blk) <- NULL
rownames(modmat.pop) <- NULL

idrop <- match(aout$dropped, colnames(modmat.fix))
idrop <- idrop[! is.na(idrop)]
modmat.fix <- modmat.fix[ , - idrop]

nfix <- ncol(modmat.fix)
nblk <- ncol(modmat.blk)
npop <- ncol(modmat.pop)

alpha.start <- aout$coefficients[match(colnames(modmat.fix),
    names(aout$coefficients))]
cee.start <- rep(0, nblk + npop)
sigma.start <- rep(1, 2)
alphaceesigma.start <- c(alpha.start, cee.start, sigma.start)

foo <- newpickle(alphaceesigma.start, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)
}
\keyword{misc}
