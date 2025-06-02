\name{penmlogl}
\alias{penmlogl}
\alias{penmlogl2}
\title{Penalized Minus Log Likelihood for Aster Models}
\usage{
penmlogl(parm, sigma, fixed, random, obj, y, origin, deriv = 2)
penmlogl2(parm, alpha, sigma, fixed, random, obj, y, origin)
}
\note{
  Not intended for use by naive users.  Use \code{\link{reaster}},
  which calls them.
}
\description{
Penalized minus log likelihood for an aster model, and its first and second
derivative.  The penalization allows for (approximate) random effects.
These functions are called inside \code{\link{pickle}},
\code{\link{pickle1}}, \code{\link{pickle2}}, \code{\link{pickle3}},
and \code{\link{reaster}}.
}
\arguments{
  \item{parm}{for \code{penmlogl}, parameter value (vector of regression
    coefficients and rescaled random effects) at which we evaluate the
    penalized log likelihood.  For \code{penmlogl2} the vector of rescaled
    random effects only (see next item).}
  \item{alpha}{the vector of fixed effects.  For \code{penmlogl2}, the
    concatenation \code{c(alpha, parm)} is the same as \code{parm} that
    is supplied to \code{pemnmlogl}.}
  \item{sigma}{vector of square roots of variance components, one component
    for each group of random effects.}
  \item{fixed}{the model matrix for fixed effects.  The number of rows
    is \code{nrow(obj$data)}.
    The number of columns is the number of fixed effects.}
  \item{random}{the model matrix or matrices for random effects.
    Each has the same number of rows as \code{fixed}.  The number of columns
    is the number of random effects in a group.  Either a matrix or a list
    of matrices.}
  \item{obj}{aster model object, the result of a call to \code{\link{aster}}.}
  \item{y}{response vector.  May be omitted, in which case \code{obj$x}
    is used.  If supplied, must be a matrix of the same dimensions as
    \code{obj$x}.}
  \item{origin}{origin of aster model.  May be omitted, in which case
    default origin (see \code{\link{aster}}) is used.  If supplied, must be
    a matrix of the same dimensions \code{obj$x}.}
  \item{deriv}{number of derivatives wanted.  Allowed values are 0, 1, or 2.}
}
\details{
  Consider an aster model with random effects and canonical parameter vector
  of the form
  \deqn{M \alpha + Z_1 b_1 + \cdots + Z_k b_k}{M alpha + Z[1] b[1] + \ldots +
  Z[k] b[k]}
  where \eqn{M} and each \eqn{Z_j}{Z[j]} are known matrices having the same
  row dimension, where \eqn{\alpha}{alpha} is a vector of unknown parameters
  (the fixed effects) and each \eqn{b_j}{b[j]} is a vector of random effects
  that are supposed to be (marginally) independent and identically distributed
  mean-zero normal with variance \code{sigma[j]^2}.

  These functions evaluate minus the \dQuote{penalized log likelihood}
  for this model, which considers the random effects as parameters but
  adds a penalization term
  \deqn{b_1^2 / (2 \sigma_1^2) + \cdots + b_k^2 / (2 \sigma_k^2)}{b[1]^2 / (2 * sigma[1]^2) + \ldots + b[k]^2 / (2 sigma[k]^2)}
  to minus the log likelihood.

  To properly deal with random effects that are zero, random effects
  are rescaled by their standard deviation.
  The rescaled random effects are
  \eqn{c_i = b_i / \sigma_i}{c[i] = b[i] / sigma[i]}.
  If \eqn{\sigma_i = 0}{sigma[i] = 0}, then the corresponding rescaled
  random effects \eqn{c_i}{c[i]} are also zero.
}
\value{
  a list containing some of the following components:
  \item{value}{minus the penalized log likelihood.}
  \item{gradient}{minus the first derivative vector of the penalized
    log likelihood.}
  \item{hessian}{minus the second derivative matrix of the penalized
    log likelihood.}
  \item{argument}{the value of the \code{parm} argument for this function.}
  \item{scale}{the vector by which \code{parm} must be scaled to obtain
    the true random effects.}
  \item{mlogl.gradient}{gradient for evaluation of log likelihood;
    \code{gradient} is this plus gradient of penalty.}
  \item{mlogl.hessian}{hessian for evaluation of log likelihood;
    \code{hessian} is this plus hessian of penalty.}
}
\seealso{For an example using this function see the example
  for \code{\link{pickle}}.}
\keyword{misc}
