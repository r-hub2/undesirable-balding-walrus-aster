\name{summary.reaster}
\alias{summary.reaster}
\alias{print.summary.reaster}
\title{Summarizing Aster Model with Random Effects Fits}
\usage{
\method{summary}{reaster}(object, standard.deviation = TRUE, \dots)

\method{print}{summary.reaster}(x, digits = max(3, getOption("digits") - 3),
      signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"reaster"}, usually, a result of a
    call to \code{\link{reaster}}.}
  \item{standard.deviation}{if \code{TRUE}, treat the parameters described
    in the \dQuote{variance components} section of the printout are square
    roots of variance components (that is, standard deviations) rather
    than the variance components themselves.  Warning: if \code{FALSE}
    so actual variance components are described, (asymptotic, approximate)
    standard errors are zero when they the variance components are zero
    (see details section below).}
  \item{x}{an object of class \code{"summary.reaster"}, usually, a result of a
    call to \code{summary.reaster}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, \dQuote{significance stars}
    are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  These functions are all \code{\link{methods}} for class \code{reaster} or
  \code{summary.reaster} objects.
}
\details{
The \code{\link{reaster}} function only does approximate maximum likelihood.
Even if it did actual maximum likelihood, standard errors would be only
approximate.  Standard errors for variance components are derived via
the delta method from standard errors for square roots of variance
components (standard deviations).  Hence P-values for variance components
and square roots of variance components do not agree exactly (although
they do asymptotically).
}
\value{
  \code{summary.reaster} returns an object of class \code{"summary.reaster"}.
}
\seealso{
  \code{\link{reaster}}, \code{\link{summary}}.
}
\keyword{models}
