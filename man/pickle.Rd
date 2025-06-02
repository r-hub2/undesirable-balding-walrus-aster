\name{pickle}
\alias{pickle}
\alias{pickle1}
\alias{pickle2}
\alias{pickle3}
\alias{makezwz}
\title{Penalized Quasi-Likelihood for Aster Models}
\usage{
pickle(sigma, parm, fixed, random, obj, y, origin, cache, \dots)
makezwz(sigma, parm, fixed, random, obj, y, origin)
pickle1(sigma, parm, fixed, random, obj, y, origin, cache, zwz,
    deriv = 0, \dots)
pickle2(alphasigma, parm, fixed, random, obj, y, origin, cache, zwz,
    deriv = 0, \dots)
pickle3(alphaceesigma, fixed, random, obj, y, origin, zwz, deriv = 0)
}
\note{
  Not intended for use by naive users.  Use \code{\link{reaster}},
  which calls them.
}
\description{
Evaluates an approximation to minus the log likelihood for an
aster model with random effects.  Uses Laplace approximation to
integrate out the random effects analytically.
The \dQuote{quasi} in the title is a misnomer in the context of aster
models but the acronym PQL for this procedure is well-established in
the generalized linear mixed models literature.
}
\arguments{
  \item{sigma}{vector of square roots of variance components,
    one component for each group of random effects.  Negative values are
    allowed; the vector of variance components is \code{sigma^2}.}
  \item{parm}{starting value for inner optimization.  Ignored if
    \code{cache$parm} exists, in which case the latter is used.
    For \code{pickle} and \code{pickle1} and \code{makezwz}, length
    is number of effects (fixed and random).
    For \code{pickle2}, length
    is number of random effects.   For all, random effects are rescaled,
    divided by the corresponding component of \code{sigma} if that is
    nonzero and equal to zero otherwise.}
  \item{alphasigma}{the concatenation of the vector of fixed effects
    and the vector of square roots of variance components.}
  \item{alphaceesigma}{the concatenation of the vector of fixed effects,
    the vector of rescaled random effects,
    and the vector of square roots of variance components.}
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
  \item{cache}{If not missing, an environment in which to cache the value
    of \code{parm} found during previous evaluations.  If supplied \code{parm}
    is taken from \code{cache}.}
  \item{zwz}{A possible value of \eqn{Z^T W Z}{t(Z) W Z}, where \eqn{Z} is the
    model matrix for all random effects and \eqn{W} is the variance matrix of
    the response.}
  \item{deriv}{Number of derivatives wanted.  For \code{pickle1}
    or \code{pickle2}, either zero or one.  For \code{pickle3},
    zero, one or two.}
  \item{\dots}{additional arguments passed to \code{\link[trust]{trust}}, which
    is used to maximize the penalized log likelihood.}
}
\details{
Define
\deqn{p(\alpha, c, \sigma) = m(a + M \alpha + Z A c) + c^T c / 2 + \log \det[A Z^T \widehat{W} Z A + I] / 2}{p(alpha, c, sigma) = m(a + M alpha + Z A c) + t(c) c / 2 + log det[A t(Z) What Z A + I] / 2}
where \eqn{m} is minus the log likelihood function of a saturated aster model,
\eqn{a} is a known vector (the \emph{offset vector} in the terminology
of \code{\link{glm}} but the \emph{origin} in the terminology
of \code{\link{aster}}), \eqn{M} is a known matrix, the model matrix for
fixed effects (the argument \code{fixed} of these functions),
\eqn{Z} is a known matrix, the model matrix for random effects
(either the argument \code{random} of these functions if it is a matrix or
\code{Reduce(cbind, random)} if \code{random} is a list of matrices),
\eqn{A} is a diagonal matrix whose diagonal is the vector
\code{rep(sigma, times = nrand)}
where \code{nrand} is \code{sapply(random, ncol)}
when \code{random} is a list of
matrices and \code{ncol(random)} when \code{random} is a matrix,
\eqn{\widehat{W}}{What} is any symmetric positive semidefinite matrix
(more on this below),
and \eqn{I} is the identity matrix.
Note that \eqn{A} is a function of \eqn{\sigma}{sigma} although the
notation does not explicitly indicate this.

Let \eqn{c^*}{cstar} denote the minimizer of
\eqn{p(\alpha, c, \sigma)} considered as a function of \eqn{c}
for fixed \eqn{\alpha}{alpha} and \eqn{\sigma}{sigma},
and let \eqn{\tilde{\alpha}}{alphatwiddle} and \eqn{\tilde{c}}{ctwiddle}
denote the (joint) minimizers of \eqn{p(\alpha, c, \sigma)} considered as
a function of \eqn{\alpha}{alpha} and \eqn{c} for fixed \eqn{\sigma}{sigma}.
Note that \eqn{c^*}{cstar} is a function of \eqn{\alpha}{alpha}
and \eqn{\sigma}{sigma} although the notation does not explicitly
indicate this.
Note that \eqn{\tilde{\alpha}}{alphatwiddle} and \eqn{\tilde{c}}{ctwiddle}
are functions of \eqn{\sigma}{sigma} (only) although the notation
does not explicitly indicate this.
Now define
\deqn{q(\alpha, \sigma) = p(\alpha, c^*, \sigma)}{q(alpha, sigma) = p(alpha, cstar, sigma)}
and
\deqn{r(\sigma) = p(\tilde{\alpha}, \tilde{c}, \sigma)}{r(sigma) = p(alphatwiddle, cstar, sigma)}
Then \code{pickle1} evaluates \eqn{r(\sigma)}{r(sigma)},
\code{pickle2} evaluates \eqn{q(\alpha, \sigma)}{q(alpha, sigma)}, and
\code{pickle3} evaluates \eqn{p(\alpha, c, \sigma)}{p(alpha, c, sigma)},
where \eqn{Z^T \widehat{W} Z}{t(Z) What Z} in the formulas above is
specified by the argument \code{zwz} of these functions.
All of these functions supply derivative (gradient) vectors if
\code{deriv = 1} is specified, and \code{pickle3} supplies the
second derivative (Hessian) matrix if \code{deriv = 2} is specified.

Let \eqn{W} denote the second derivative function of \eqn{m}, that is,
\eqn{W(\varphi)}{W(phi)} is the second derivative matrix of the function
\eqn{m} evaluated at the point \eqn{\varphi}{phi}.  The idea is that
\eqn{\widehat{W}}{What} should be approximately the value of
\eqn{W(a + M \hat{\alpha} + Z \widehat{A} \hat{c})}{W(a + M alphahat + Z Ahat chat)},
where \eqn{\hat{\alpha}}{alphahat}, \eqn{\hat{c}}{chat},
and \eqn{\hat{\sigma}}{sigmahat} are the (joint) minimizers of \eqn{p(\alpha, c, \sigma)}{p(alpha, c, sigma)}
and \eqn{\widehat{A} = A(\hat{\sigma})}{Ahat = A(sigmahat)}.
In aid of this, the function \code{makezwz}
evaluates
\eqn{Z^T W(a + M \alpha + Z A c) Z}{t(Z) W(a + M alpha + Z A c) Z}
for any \eqn{\alpha}{alpha}, \eqn{c}, and \eqn{\sigma}{sigma}.

\code{pickle} evaluates the function
\deqn{s(\sigma) = m(a + M \tilde{\alpha} + Z A \tilde{c}) + \tilde{c}^T \tilde{c} / 2 + \log \det[A Z^T W(a + M \tilde{\alpha} + Z A \tilde{c}) Z A + I]}{s(sigma) = m(a + M alphatwiddle + Z A ctwiddle) + t(ctwiddle) ctwiddle / 2 + log det[A t(Z) W(a + M alphatwiddle + Z A ctwiddle) Z A + I]}
no derivatives can be computed because no derivatives of the function
\eqn{W} are computed for aster models.

The general idea is the one uses \code{pickle} with a no-derivative
optimizer, such as the \code{"Nelder-Mead"} method of the \code{optim}
function to get a crude estimate of \eqn{\hat{\sigma}}{sigmahat}.
Then one uses \code{\link[trust]{trust}} with objective
function \code{\link{penmlogl}} to estimate the corresponding
\eqn{\hat{\alpha}}{alphahat} and \eqn{\hat{c}}{chat} (example below).
Then one use \code{makezwz} to produce the corresponding \code{zwz}
(example below).
These estimates can be improved using \code{\link[trust]{trust}} with objective
function \code{pickle3} using this \code{zwz} (example below), and this
step may be iterated until convergence.
Finally, \code{\link{optim}} is used with objective function \code{pickle2}
to estimate the Hessian matrix of \eqn{q(\alpha, \sigma)}{q(alpha, sigma)},
which is approximate observed information because
\eqn{q(\alpha, \sigma)}{q(alpha, sigma)} is approximate minus log likelihood.
}
\value{
  For \code{pickle}, a scalar, minus the
  (PQL approximation of) the log likelihood.
  For \code{pickle1} and \code{pickle2}, a list having components \code{value}
  and \code{gradient} (present only when \code{deriv = 1}).
  For \code{pickle3}, a list having components \code{value},
  \code{gradient} (present only when \code{deriv >= 1}),
  and \code{hessian} (present only when \code{deriv = 2}).
}
\examples{
data(radish)
library(trust)

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

### try penmlogl

sigma.start <- c(1, 1)

alpha.start <- aout$coefficients[match(colnames(modmat.fix),
    names(aout$coefficients))]
parm.start <- c(alpha.start, rep(0, nblk + npop))

tout <- trust(objfun = penmlogl, parm.start, rinit = 1, rmax = 10,
    sigma = sigma.start, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout)
tout$converged

### crude estimate of variance components

eff.blk <- tout$argument[seq(nfix + 1, nfix + nblk)]
eff.pop <- tout$argument[seq(nfix + nblk + 1, nfix + nblk + npop)]
sigma.crude <- sqrt(c(var(eff.blk), var(eff.pop)))

### try optim and pickle

cache <- new.env(parent = emptyenv())
oout <- optim(sigma.crude, pickle, parm = tout$argument,
    fixed = modmat.fix, random = list(modmat.blk, modmat.pop),
    obj = aout, cache = cache)
oout$convergence == 0
### estimated variance components
oout$par^2

### get estimates of fixed and random effects

tout <- trust(objfun = penmlogl, tout$argument, rinit = 1, rmax = 10,
    sigma = oout$par, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout, fterm = 0)
tout$converged

sigma.better <- oout$par
alpha.better <- tout$argument[1:nfix]
c.better <- tout$argument[- (1:nfix)]
zwz.better <- makezwz(sigma.better, parm = c(alpha.better, c.better),
    fixed = modmat.fix, random = list(modmat.blk, modmat.pop), obj = aout)

### get better estimates

objfun <- function(alphaceesigma, zwz)
    pickle3(alphaceesigma, fixed = modmat.fix,
    random = list(modmat.blk, modmat.pop), obj = aout, zwz = zwz, deriv = 2)
tout <- trust(objfun, c(alpha.better, c.better, sigma.better),
    rinit = 1, rmax = 10, zwz = zwz.better)
tout$converged
alpha.mle <- tout$argument[1:nfix]
c.mle <- tout$argument[nfix + 1:(nblk + npop)]
sigma.mle <- tout$argument[nfix + nblk + npop + 1:2]
zwz.mle <- makezwz(sigma.mle, parm = c(alpha.mle, c.mle),
    fixed = modmat.fix, random = list(modmat.blk, modmat.pop), obj = aout)
### estimated variance components
sigma.mle^2

### preceding step can be iterated "until convergence"

### get (approximate) Fisher information

objfun <- function(alphasigma) pickle2(alphasigma, parm = c.mle,
    fixed = modmat.fix, random = list(modmat.blk, modmat.pop),
    obj = aout, zwz = zwz.mle)$value
gradfun <- function(alphasigma) pickle2(alphasigma, parm = c.mle,
    fixed = modmat.fix, random = list(modmat.blk, modmat.pop),
    obj = aout, zwz = zwz.mle, deriv = 1)$gradient
oout <- optim(c(alpha.mle, sigma.mle), objfun, gradfun, method = "BFGS",
    hessian = TRUE)
oout$convergence == 0
fish <- oout$hessian
}
\keyword{misc}
