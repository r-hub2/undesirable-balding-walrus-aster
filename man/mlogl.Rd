\name{mlogl}
\alias{mlogl}
\title{Minus Log Likelihood for Aster Models}
\usage{
mlogl(parm, pred, fam, x, root, modmat, deriv = 0,
    type = c("unconditional", "conditional"), famlist = fam.default(),
    origin, origin.type = c("model.type", "unconditional", "conditional"))
}
\description{
Minus the Log Likelihood for an Aster model, and its first and second
derivative.  This function is called inside R function \code{\link{aster}}.
Users generally do not need to call it directly.
}
\arguments{
  \item{parm}{parameter value (vector of regression coefficients)
    where we evaluate the log likelihood, etc.
    We also refer to \code{length(parm)} as \code{ncoef}.}
  \item{pred}{integer vector determining the graph.
    \code{pred[j]} is the index of the predecessor of
    the node with index \code{j} unless the predecessor is a root
    node, in which case \code{pred[j] == 0}.
    We also refer to \code{length(pred)} as \code{nnode}.
    This argument is required to satisfy \code{pred[j] < j} for all \code{j}.}
  \item{fam}{an integer vector of length \code{nnode} determining
    the exponential family structure of the aster model.  Each element
    is an index into the vector of family specifications given by
    the argument \code{famlist}.}
  \item{x}{the response.  If a matrix, rows are individuals, and columns are
    variables (nodes of graphical model).  So \code{ncol(x) == nnode} and
    we also refer to \code{nrow(x)} as \code{nind}.  If not a matrix, then
    \code{x} must be as if it were such a matrix and then dimension
    information removed by \code{x = as.numeric(x)}.}
  \item{root}{A matrix or vector like \code{x}.
    Data \code{root[i, j]} is the data for a root node that is
    the predecessor of the response \code{x[i, j]}
    and is ignored when \code{pred[j] > 0}.}
  \item{modmat}{a three-dimensional array, \code{nind} by \code{nnode} by
    \code{ncoef}, the model matrix.  Or a matrix, \code{nind * nnode} by
    \code{ncoef} (when \code{x} and \code{root} are one-dimensional
    of length \code{nind * nnode}).}
  \item{deriv}{derivatives wanted: 0, 1, or 2.}
  \item{type}{type of model.  The value of this argument can be abbreviated.}
  \item{famlist}{a list of family specifications (see \code{\link{families}}).}
  \item{origin}{Distinguished point in parameter space.  May be missing,
    in which case an unspecified default is provided.  See \code{\link{aster}}
    for further explanation.}
  \item{origin.type}{Parameter space in which specified distinguished point
    is located.  If \code{"conditional"} then argument \code{"origin"} is
    a conditional canonical parameter value.
    If \code{"unconditional"} then argument \code{"origin"} is
    an unconditional canonical parameter value.
    If \code{"model.type"} then the type is taken from argument \code{"type"}.
    The value of this argument can be abbreviated.}
}
\value{
  a list containing some of the following components:
  \item{value}{minus the log likelihood.}
  \item{gradient}{minus the first derivative vector of the log likelihood
    (minus the score).}
  \item{hessian}{minus the second derivative matrix of the log likelihood
    (observed Fisher information).}
}
\keyword{misc}
