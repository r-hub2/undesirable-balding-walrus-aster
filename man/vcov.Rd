\name{vcov.aster}
\alias{vcov.aster}
\alias{vcov.reaster}
\title{Calculate Asymptotic Variance-Covariance Matrix for a Fitted Model
    Object}
\usage{
\method{vcov}{aster}(object, info = c("expected", "observed"), \dots)
\method{vcov}{reaster}(object, standard.deviation = TRUE, re.too = FALSE,
  complete = TRUE, \dots)
}
\arguments{
  \item{object}{an object of class \code{"aster"} or class \code{"reaster"},
    usually the object returned by a call to R function \code{\link{aster}}
    or R function \code{\link{reaster}}.}
  \item{info}{the type of Fisher information used to compute the
    variance-covariance matrix.}
  \item{standard.deviation}{if \code{FALSE}, then some of the parameters are
     variance components, if \code{TRUE}, then the square roots of these
     variance components.}
  \item{re.too}{if \code{TRUE}, treat the estimated random effects vector,
     the mode of the conditional distribution of random effects given
     observed data, as part of the estimates vector described by the
     returned variance-covariance matrix.  If also \code{standard.deviation}
     is \code{TRUE}, then these estimates are standardized.}
  \item{complete}{If \code{TRUE}, the full variance-covariance matrix,
     including entries for parameters estimated to be zero (for which the
     entries will also be zero).}
  \item{\dots}{required by the generic function.  Ignored by these methods.}
}
\description{
  \code{vcov()} is a generic function and functions with names beginning
  in \code{vcov.} will be methods for this function.
}
\value{
  A matrix of the estimated approximate, large sample, covariances between
  the parameter estimates.

  Entries for parameters estimated to be zero are zero, not because there
  is no sampling variability but rather because the sampling distribution
  is highly non-normal and one-sided so covariances are useless.

  This returned matrix has dimension attributes, like any R matrix.
  For the function \code{vcov.reaster}, it also
  has attributes, which are some of \code{is.alpha}, \code{is.b},
  \code{is.nu}, \code{is.c}, and \code{is.sigma}.  These are all logical
  vectors that can serve as index vectors for the matrix.
  \itemize{
    \item \code{is.alpha} extracts elements of the variance-covariance
      matrix for estimates of fixed effects.
    \item \code{is.nu} extracts elements of the variance-covariance
      matrix for estimates of variance components
      (if \code{standard.deviation == FALSE} was specified).
    \item \code{is.b} extracts elements of the variance-covariance
      matrix for estimates of random effects
      (if \code{standard.deviation == FALSE & re.too = TRUE} was specified).
    \item \code{is.nu} extracts elements of the variance-covariance
      matrix for estimates of square roots of variance components
      (if \code{standard.deviation == TRUE} was specified).
    \item \code{is.c} extracts elements of the variance-covariance
      matrix for standardized estimates of random effects
      (if \code{standard.deviation == TRUE & re.too = TRUE} was specified).
  }
  See \code{\link{reaster}} for more about these parameterizations.
}
\examples{
data(radish2)

pred <- c(0,1,2)
fam <- c(1,3,2)

rout <- reaster(resp ~ varb + fit : (Site * Region),
    list(block = ~ 0 + fit : Block, pop = ~ 0 + fit : Pop),
    pred, fam, varb, id, root, data = radish2)

vout <- vcov(rout, standard.deviation = FALSE, re.too = TRUE)
attributes(vout)
is.alpha <- attr(vout, "is.alpha")
vout[is.alpha, is.alpha]
}
\seealso{
  \code{\link{aster}}, \code{\link{reaster}}, \code{\link{vcov}}
}
\keyword{models}
