\name{objfun}
\alias{objfun}
\alias{objfun.factory}
\title{Function Factory that Makes Objective Function for Aster Models
   with Random Effects}
\concept{regression}
\concept{exponential family}
\concept{graphical model}
\concept{mixed model}
\concept{random effects}
\description{
  A function factory that returns an R function that evaluates the
  objective function for fitting aster models with random effects
}
\usage{
objfun.factory(fixed, random, response, pred, fam, root, zwz,
    famlist = fam.default(), offset, standard.deviation = FALSE,
    deriv = 0)
}
\arguments{

  \item{fixed}{the model matrix for fixed effects.}

  \item{random}{either a model matrix or list of model matrices.
    Each model matrix operates on the random effects for one
    variance component. See Details section for more.}

  \item{response}{a vector specifying the response.  Length
    is the row dimension of all model matrices.}

  \item{pred}{an integer vector of length \code{nnode} determining
    the dependence graph of the aster model.  \code{pred[j]} is
    the index of the predecessor of
    the node with index \code{j} unless the predecessor is a root
    node, in which case \code{pred[j] == 0}.  See details section
    of \code{\link{aster}} for further requirements.}

  \item{fam}{an integer vector of length \code{nnode} determining
    the exponential family structure of the aster model.  Each element
    is an index into the vector of family specifications given by
    the argument \code{famlist}.}

  \item{root}{a vector whose length is the row dimension of all model
    matrices.  For nodes whose predecessors are root nodes specifies the
    value of the constant at that root node.  Typically the vector having
    all components equal to one.  Always positive integer valued.  If
    components greater than one, then graph for that "individual" in scare
    quotes is actually for data for multiple individuals.}

  \item{zwz}{a positive definite symmetric matrix (neither symmetry nor
    positive definiteness is checked, but incorrect answers with no error
    or warning will result if condition is not met)
    whose dimensions are the number of random effects.  Should be
    \eqn{Z^T W(\varphi) Z}{t(Z) \%*\% W(phi) \%*\% Z} where \eqn{Z} is
    the model matrix for random effects and \eqn{W(\,\cdot\,)}{W()} is
    minus the second derivative of the saturated model aster log likelihood
    with respect to its unconditional canonical parameter.  But it does
    not have to be.  Use R function \code{\link{makezwz}} to make this.}

  \item{famlist}{a list of family specifications (see \code{\link{families}}).}

  \item{offset}{a vector whose length is the row dimension of all model
    matrices.  Distinguished point in parameter space.  May be missing,
    in which case an unspecified default is provided.  See details of
    \code{\link{aster}} for further explanation, where it is called
    \dQuote{origin} rather than \dQuote{offset}.}

  \item{standard.deviation}{scalar logical.  If \code{FALSE}, use the
    parameterization \code{alpha}, \code{b}, \code{nu}, fixed effects,
    random effects, and variance components, respectively.  If \code{TRUE},
    use the parameterization \code{alpha}, \code{c}, \code{sigma}.
    See Details for more about these parameterizations.}
   
  \item{deriv}{0, 1, or 2.  Number of derivatives.  Zero means only the
    function value, one means gradient, two means Hessian.}
}
\details{
  See the help page for the function \code{\link{aster}} for specification
  of aster models.  This function only fits unconditional aster models
  (those with default values of the \code{aster} function arguments
  \code{type} and \code{origin.type}.

  The only difference between this function and functions \code{aster}
  and \code{mlogl} is
  that some effects (also called coefficients, also called canonical affine
  submodel canonical parameters) are treated as random.
  The unconditional canonical
  parameter vector of the (saturated) aster model is treated as
  an affine function of fixed and random effects
  \deqn{\varphi = a + M \alpha + \sum_{i = 1}^k Z_i b_i}{phi = a + M alpha
  + Z[1] b[1] + \dots + Z[k] b[k]}
  where \eqn{M} and the \eqn{Z_i}{Z[i]} are model matrices specified by
  the arguments \code{fixed} and \code{random}, where \eqn{\alpha}{alpha}
  is a vector of
  fixed effects and each \eqn{b_i}{b[i]} is a vector of random
  effects that are assumed to be (marginally) normally distributed with
  mean vector zero and variance matrix \eqn{\nu_i}{nu[i]} times
  the identity matrix.  The random effects vectors \eqn{b_i}{b[i]} are
  assumed independent.

  In what follows, \eqn{Z} is the single matrix produced
  by \code{cbind}'ing the \eqn{Z_i}{Z[i]},
  and \eqn{b} is the single vector produced by \code{c}'ing
  the \eqn{b_i}{b[i]},
  and \eqn{\nu} is the single vector produced by \code{c}'ing the
  \eqn{\nu_i}{nu[i]}, and \eqn{D} is the variance-covariance matrix
  of \eqn{b}, which is diagonal with all diagonal components equal to
  some component of \eqn{\nu}{nu}.

  We can write
  \deqn{D = \sum_k \nu_k E_k}{D = nu[1] E[1] + ... + nu[k] E[k]}
  where the nus are unknown parameters to estimate and the Es are known
  diagonal matrices whose components are zero-or-one-valued having the
  property that they sum to the identity matrix and multiply to the zero
  matrix (so each component of \eqn{b} is exactly one component of
  \eqn{\nu}{nu}).

  The R function returned by this function evaluates
  \deqn{p_W(\alpha, b, \nu) =
  - l(\varphi) + \tfrac{1}{2} b^T D^{-1} b + \tfrac{1}{2}
  \log\mathop{\rm det}\left(Z^T W Z D + \text{Id} \right)}{p_W(alpha, b, nu)
  = - l(phi)
  + b \%*\% solve(D) \%*\% b / 2
  + log(det(t(Z) \%*\% W \%*\% Z \%*\% D + Id)) / 2}
  where \eqn{l} is the log likelihood for the saturated aster model specified
  by arguments \code{pred}, \code{fam}, and \code{famlist},
  where \eqn{\varphi}{phi} is given by the displayed equation above,
  and where \eqn{Z^T W Z}{t(Z) \%*\% W \%*\% Z} is argument \code{zwz}.
  This is equation (7) in Geyer, et al. (2013).

  The displayed equation above only defines \eqn{p_W}{p_W} when all of the
  variance components (components of \eqn{\nu}{nu} and diagonal components
  of \eqn{D}) are strictly positive because otherwise \eqn{D} is not
  invertible.  In order to be lower semicontinuous, the R function returned
  by this function returns \code{Inf} whenever any component of \eqn{\nu}{nu}
  is nonpositive, except it returns zero if \eqn{b_i = 0}{b[i] = 0} whenever
  the corresponding diagonal element of \eqn{D} is zero.

  When argument \code{standard.deviation} is \code{TRUE}, the arguments
  to the objective function are replaced by the change-of-variables
  \deqn{\nu = \sigma^2}{nu = sigma^2} and \deqn{b = A c}{b = A \%*\% c}
  where the first equation works componentwise (as in R vector arithmetic),
  components of the vector \eqn{\sigma}{sigma} are allowed to be
  negative, and \eqn{D = A^2} so
  \deqn{A = \sum_k \sigma_k E_k}{D = sigma[1] E[1] + ... + sigma[k] E[k]}

  When the \eqn{(\alpha. c, \sigma)}{(alpha, c, sigma)} variables are used,
  \eqn{p_W}{p_W} is continuous and infinitely differentiable for all values
  of these variables.
 
  When the \eqn{(\alpha. b, \nu)}{(alpha, b, nu)} variables are used,
  \eqn{p_W}{p_W} is continuous and infinitely differentiable for all values
  of \eqn{\alpha}{alpha} and \eqn{b} but only strictly positive values
  \eqn{\nu}{nu}.  When some components of \eqn{\nu}{nu} are zero,
  then subderivatives exist, but not derivatives (Geyer, et al., 2013).
}

\value{
  returns an R function having one argument which is a numeric vector
  of all the parameters and random effects in the order
  \eqn{(\alpha, b, \nu)}{(alpha, b, nu)} for one set of variables or
  \eqn{(\alpha, c, \sigma)}{(alpha, c, \sigma)} for other set.

  The R function returned by this function evaluates the function \eqn{p_W}
  defined in the Details section and returns a list with components named
  \code{value}, \code{gradient}, and \code{hessian} in case the \code{deriv}
  argument is 2.  The returned list has only the first two of these in case
  the \code{deriv} argument is 1.
  The returned list has only the first one of these in case
  the \code{deriv} argument is 0.

  In case the \eqn{(\alpha, b, \nu)}{(alpha, b, nu)} parameterization is
  being used and some components of \eqn{\nu}{nu} are zero, the derivatives
  returned are for the terms of the objective function that are differentiable.
}

\section{Warning about Negative Binomial}{
The negative binomial and truncated negative binomial families are fundamentally
incompatible with random effects.  The reason is that the canonical parameter
space for a one-parameter negative binomial or truncated negative binomial
is the negative half line.  Thus the conditional canonical parameter
\eqn{\theta}{theta} for such a node must be negative valued.  The aster
transform is so complicated that it is unclear what the corresponding
constraint on the unconditional canonical parameter \eqn{\varphi}{phi} is,
but there is a constraint: its parameter space is not the whole real line.
A normal random effect, in contrast, does have support the whole real line.
It wants to make parameters that are constrained to have any real number.
The code only warns about this situation, because if the random effects do
not influence any negative binomial or truncated negative binomial nodes
of the graph, then there would be no problem.
}
\section{Warning about Individual Random Effects}{
The Breslow-Clayton approximation assumes the complete data log likelihood
is approximately quadratic considered as a function of random effects only.
This will be the case by the law of large numbers if the number of individuals
is much larger than the number of random effects.  Thus Geyer, et al. (2013)
warn against trying to put a random effect for each individual in the model.
If you do that, the code will try to fit the model, but it will take forever
and no theory says the results will make any sense.
}
\references{
Breslow, N. E., and Clayton, D. G. (1993).
Approximate Inference in Generalized Linear Mixed Models.
\emph{Journal of the American Statistical Association}, \bold{88}, 9--25.
\doi{10.1080/01621459.1993.10594284}.

Geyer, C. J., Ridley, C. E., Latta, R. G., Etterson, J. R.,
    and Shaw, R. G. (2013)
Local Adaptation and Genetic Effects on Fitness: Calculations
    for Exponential Family Models with Random Effects.
\emph{Annals of Applied Statistics}, \bold{7}, 1778--1795.
\doi{10.1214/13-AOAS653}.
}
\examples{
library(aster)
data(radish)
pred <- c(0,1,2)
fam <- c(1,3,2)
rout <- reaster(resp ~ varb + fit : (Site * Region),
    list(block = ~ 0 + fit : Block, pop = ~ 0 + fit : Pop),
    pred, fam, varb, id, root, data = radish)
objfun <- with(rout, objfun.factory(fixed, random, response,
    obj$pred, obj$fam, as.vector(obj$root), zwz))
theta.hat <- with(rout, c(alpha, b, nu))
all.equal(objfun(theta.hat)$value, rout$deviance / 2)
}
\keyword{models}
\keyword{regression}
