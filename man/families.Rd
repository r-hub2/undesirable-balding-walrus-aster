\name{families}
\alias{families}
\alias{fam.bernoulli}
\alias{fam.poisson}
\alias{fam.truncated.poisson}
\alias{fam.negative.binomial}
\alias{fam.truncated.negative.binomial}
\alias{fam.normal.location}
\alias{fam.default}
\alias{famfun}
\title{Families for Aster Models}
\usage{
fam.bernoulli()
fam.poisson()
fam.truncated.poisson(truncation)
fam.negative.binomial(size)
fam.truncated.negative.binomial(size, truncation)
fam.normal.location(sd)
fam.default()
famfun(fam, deriv, theta)
}
\description{
Families (response models) known to the package.
These functions construct simple family specifications used
in specifying models for \code{\link{aster}} and \code{\link{mlogl}}.
They are mostly for convenience, since the specifications are easy
to construct by hand.
}
\arguments{
  \item{truncation}{the truncation point, called \eqn{k} in the details
  section below.}
  \item{size}{the sample size.  May be non-integer.}
  \item{sd}{the standard deviation.  May be non-integer.}
  \item{fam}{a family specification, which is a list of class \code{"astfam"}
      containing, at least one element named \code{"name"} and perhaps other
      elements specifying hyperparameters.}
  \item{deriv}{derivative wanted: 0, 1, or 2.}
  \item{theta}{value of the canonical parameter.}
}
\details{
Currently implemented families are
  \describe{
    \item{\code{"bernoulli"}}{Bernoulli.  The mean value parameter
      \eqn{\mu}{mu} is the success probability.  The canonical parameter is
      \eqn{\theta = \log(\mu) - \log(1 - \mu)}{theta = log(mu) - log(1 - mu)},
      also called logit of \eqn{\mu}{mu}.  The first derivative of the
      cumulant function has the value \eqn{\mu}{mu} and the second derivative
      of the cumulant function has the value \eqn{\mu (1 - \mu)}{mu (1 - mu)}.}
    \item{\code{"poisson"}}{Poisson.  The mean value parameter
      \eqn{\mu}{mu} is the mean of the Poisson distribution.
      The canonical parameter is \eqn{\theta = \log(\mu)}{theta = log(mu)}.
      The first and second derivatives of the cumulant function both have
      the value \eqn{\mu}{mu}.}
    \item{\code{"truncated.poisson"}}{Poisson conditioned on being
      strictly greater than \eqn{k}, specified by
      the argument \code{truncation}.
      Let \eqn{\mu}{mu} be the mean of the corresponding untruncated
      Poisson distribution.  Then the canonical parameters for both
      truncated and untruncated distributions are the same
      \eqn{\theta = \log(\mu)}{theta = log(mu)}.
      Let \eqn{Y} be a Poisson random variable having the same mean parameter
      as this distribution, and define
      \deqn{\beta = \frac{\Pr\{Y > k + 1\}}{\Pr\{Y = k + 1\}}}{beta = Pr(Y > k + 1) / Pr(Y = k + 1)}
      Then the mean value parameter and first derivative of the
      cumulant function of this distribution has the value
      \deqn{\tau = \mu + \frac{k + 1}{1 + \beta}}{tau = mu + (k + 1) / (1 + beta)}
      and the second derivative of the cumulant function has the value
      \deqn{\mu \left[ 1 - \frac{k + 1}{1 + \beta} \left( 1 - \frac{k + 1}{\mu} \cdot \frac{\beta}{1 + \beta} \right) \right]}{mu [ 1 - (k + 1) / (1 + beta) ( 1 - (k + 1) / mu * beta / (1 + beta) ) ]}.}
    \item{\code{"negative.binomial"}}{Negative binomial.  The size parameter
      \eqn{\alpha}{alpha} may be noninteger, meaning the cumulant function
      is \eqn{\alpha}{alpha} times the cumulant function of the geometric
      distribution.  The mean value parameter \eqn{\mu}{mu} is the mean of
      the negative binomial distribution.  The success probability parameter
      is \deqn{p = \frac{\alpha}{\mu + \alpha}.}{p = alpha / (mu + alpha).}
      The canonical parameter
      is \eqn{\theta = \log(1 - p)}{theta = log(1 - p)}.
      Since \eqn{1 - p < 1}, the canonical parameter space is restricted,
      the set of \eqn{\theta}{theta} such that \eqn{\theta < 0}{theta < 0}.
      This is, however, a
      regular exponential family (the log likelihood goes to minus infinity
      as \eqn{\theta}{theta} converges to the boundary of the parameter space,
      so the constraint \eqn{\theta < 0}{theta < 0} plays no role in maximum
      likelihood estimation so long as the optimization software is not too
      stupid.  There will be no problems so long as the default optimizer
      (\code{\link[trust]{trust}}) is used.  Since zero is not in the canonical
      parameter space a negative default origin is used.  The first derivative
      of the cumulant function has the value
      \deqn{\mu = \alpha \frac{1 - p}{p}}{mu = alpha (1 - p) / p}
      and the second derivative has the value
      \deqn{\alpha \frac{1 - p}{p^2}.}{alpha (1 - p) / p^2.}}
    \item{\code{"truncated.negative.binomial"}}{Negative binomial conditioned
      on being strictly greater than \eqn{k}, specified by
      the argument \code{truncation}.
      Let \eqn{p} be the success probability parameter of the corresponding
      untruncated negative binomial distribution.  Then the canonical
      parameters for both
      truncated and untruncated distributions are the same
      \eqn{\theta = \log(1 - p)}{theta = log(1 - p)}, and consequently
      the canonical parameter spaces are the same,
      the set of \eqn{\theta}{theta} such that \eqn{\theta < 0}{theta < 0},
      and both models are regular exponential families.
      Let \eqn{Y} be an untruncated negative binomial random variable having
      the same size and success probability parameters as this distribution.
      and define
      \deqn{\beta = \frac{\Pr\{Y > k + 1\}}{\Pr\{Y = k + 1\}}}{beta = Pr(Y > k + 1) / Pr(Y = k + 1)}
      Then the mean value parameter and first derivative of the
      cumulant function of this distribution has the value
      \deqn{\tau = \mu + \frac{k + 1}{p (1 + \beta)}}{tau = mu + (k + 1) / (p (1 + beta))}
      and the second derivative is too complicated to write here (the
      formula can be found in the vignette \code{trunc.pdf}.}
    \item{\code{"normal.location"}}{Normal, unknown mean, known variance.
      The sd (standard deviation) parameter
      \eqn{\sigma}{sigma} may be noninteger, meaning the cumulant function
      is \eqn{\sigma^2}{sigma^2} times the cumulant function of the standard
      normal distribution.  The mean value parameter \eqn{\mu}{mu} is the
      mean of the normal distribution.
      The canonical parameter
      is \eqn{\theta = \mu / \sigma^2}{theta = mu / sigma^2}.
      The first derivative of the cumulant function has the value
      \deqn{\mu = \sigma^2 \theta}{mu = sigma^2 theta}
      and the second derivative has the value
      \deqn{\sigma^2.}{sigma^2.}}
  }
}
\value{
    For all but \code{fam.default},
    a list of class \code{"astfam"} giving name and values of any
    hyperparameters.
    For \code{fam.default},
    a list each element of which is of class \code{"astfam"}.
    The list of families which were hard coded in earlier versions of the
    package.
}
\seealso{\code{\link{aster}} and \code{\link{mlogl}}}
\examples{
### mean of poisson with mean 0.2
famfun(fam.poisson(), 1, log(0.2))
### variance of poisson with mean 0.2
famfun(fam.poisson(), 2, log(0.2))
### mean of poisson with mean 0.2 conditioned on being nonzero
famfun(fam.truncated.poisson(trunc = 0), 1, log(0.2))
### variance of poisson with mean 0.2 conditioned on being nonzero
famfun(fam.truncated.poisson(trunc = 0), 2, log(0.2))
}
\keyword{misc}
