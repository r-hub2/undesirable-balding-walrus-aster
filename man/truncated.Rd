\name{truncated}
\alias{truncated}
\alias{rktp}
\alias{rktnb}
\alias{rnzp}
\title{K-Truncated Distributions}
\description{
  Random generation for the \eqn{k}-truncated Poisson distribution
  or for the \eqn{k}-truncated negative binomial distribution, where
  \dQuote{\eqn{k}-truncated} means conditioned on being strictly greater
  than \eqn{k}.  If \code{xpred} is not one, then the random variate is
  the sum of \code{xpred} such random variates.
}
\usage{
rktp(n, k, mu, xpred = 1)
rktnb(n, size, k, mu, xpred = 1)
rnzp(n, mu, xpred = 1)
}
\arguments{
  \item{n}{number of random values to return.  If \code{length(n) > 1},
    the length is taken to be the number required.}
  \item{size}{the size parameter for the negative binomial distribution.}
  \item{k}{truncation limit.}
  \item{xpred}{number of trials.}
  \item{mu}{vector of positive means.}
}
\details{
  \code{rktp} simulates \eqn{k}-truncated Poisson random variates.
  \code{rktnb} simulates \eqn{k}-truncated negative binomial random variates.
  \code{rnzp} simulates zero-truncated Poisson random variates
  (maintained only for backward compatibility, it now calls \code{rktp}).
}
\value{
  a vector of random deviates.
}
\seealso{\code{\link{families}}}
\examples{
rktp(10, 2, 0.75)
rktnb(10, 2.222, 2, 0.75)
}
\keyword{distribution}
