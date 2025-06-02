\name{quickle}
\alias{quickle}
\title{Penalized Quasi-Likelihood for Aster Models}
\usage{
quickle(alphanu, bee, fixed, random, obj, y, origin, zwz, deriv = 0)
}
\note{
  Not intended for use by naive users.  Use \code{\link{summary.reaster}},
  which calls it.
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
  \item{alphanu}{the parameter vector value where the function is evaluated,
    a numeric vector, see details.}
  \item{bee}{the random effects vector that is used as the starting point
    for the inner optimization, which maximizes the penalized log likelihood
    to find the optimal random effects vector matching \code{alphanu}.}
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
    the response.  See details.  Typically constructed by
    the function \code{\link{makezwz}}.}
  \item{deriv}{Number of derivatives wanted, zero, one, or two.}
}
\details{
Define
\deqn{p(\alpha, b, \nu) = m(a + M \alpha + Z b) + {\textstyle \frac{1}{2}} b^T D^{- 1} b + {\textstyle \frac{1}{2}} \log \det[Z^T W Z D + I]}{p(alpha, b, nu) = m(a + M alpha + Z b) + t(b) D^(- 1) b / 2 + log det[t(Z) W Z D + I] / 2}
where \eqn{m} is minus the log likelihood function of a saturated aster model,
where \eqn{a} is a known vector (the \emph{offset vector} in the terminology
of \code{\link{glm}} but the \emph{origin} in the terminology
of \code{\link{aster}}),
where \eqn{M} is a known matrix, the model matrix for fixed effects
(the argument \code{fixed} of this function),
where \eqn{Z} is a known matrix, the model matrix for random effects
(either the argument \code{random} of this function if it is a matrix or
\code{Reduce(cbind, random)} if \code{random} is a list of matrices),
where \eqn{D} is a diagonal matrix whose diagonal is the vector
\code{rep(nu, times = nrand)}
where \code{nrand} is \code{sapply(random, ncol)}
when \code{random} is a list of
matrices and \code{ncol(random)} when \code{random} is a matrix,
where \eqn{W} is an arbitrary symmetric positive semidefinite matrix
(\eqn{Z^T W Z}{t(Z) W Z} is the argument \code{zwz} of this function),
and where \eqn{I} is the identity matrix.
Note that \eqn{D} is a function of \eqn{\nu}{nu}
although the notation does not explicitly indicate this.

The argument \code{alphanu} of this function is the concatenation
of the parameter vectors \eqn{\alpha}{alpha} and \eqn{\nu}{\nu}.
The argument \code{bee} of this function is a possible value of \eqn{b}.
The length of \eqn{\alpha}{alpha} is the column dimension of \eqn{M}.
The length of \eqn{b} is the column dimension of \eqn{Z}.
The length of \eqn{\nu} is the length of the argument \code{random}
of this function if it is a list and is one otherwise.

Let \eqn{b^*}{bstar} denote the minimizer
of \eqn{p(\alpha, b, \nu)}{p(alpha, b, nu)} considered as a function of
\eqn{b} for fixed \eqn{\alpha}{alpha} and \eqn{\nu}{nu}, so \eqn{b^*}{bstar}
is a function of \eqn{\alpha}{alpha} and \eqn{\nu}{nu}.
This function evaluates
\deqn{q(\alpha, \nu) = p(\alpha, b^*, \nu)}{q(alpha, nu) = p(alpha, bstar, nu)}
and its gradient vector and Hessian matrix (if requested).
Note that \eqn{b^*}{bstar} is a function of \eqn{\alpha}{alpha}
and \eqn{\nu}{nu} although the notation does not explicitly indicate this.
}
\value{
  a list with some of the following components: \code{value}, \code{gradient},
  \code{hessian}, \code{alpha}, \code{bee}, \code{nu}.  The first three are
  the requested derivatives.  The second three are the corresponding parameter
  values: \code{alpha} and \code{nu} are the corresponding parts of the
  argument \code{alphanu}, the value of \code{bee} is the result of the inner
  optimization (\eqn{b^*}{bstar} in the notation in details),
  not the argument \code{bee} of this function.
}
\examples{
data(radish)

pred <- c(0,1,2)
fam <- c(1,3,2)

rout <- reaster(resp ~ varb + fit : (Site * Region),
    list(block = ~ 0 + fit : Block, pop = ~ 0 + fit : Pop),
    pred, fam, varb, id, root, data = radish)

alpha.mle <- rout$alpha
bee.mle <- rout$b
nu.mle <- rout$sigma^2
zwz.mle <- rout$zwz
obj <- rout$obj
fixed <- rout$fixed
random <- rout$random
alphanu.mle <- c(alpha.mle, nu.mle)

qout <- quickle(alphanu.mle, bee.mle, fixed, random, obj,
    zwz = zwz.mle, deriv = 2)
}
\keyword{misc}
