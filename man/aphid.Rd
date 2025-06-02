\name{aphid}
\docType{data}
\alias{aphid}
\title{Life History Data on Uroleucon rudbeckiae}
\description{
  Data on life history traits for the brown ambrosia aphid
  \emph{Uroleucon rudbeckiae}
}
\usage{data(aphid)}
\format{
  A data frame with records for 18 insects.
  Data are already in \dQuote{long} format; no need to reshape.
   \describe{
    \item{resp}{Response vector.}
    \item{varb}{Categorical.  Gives node of graphical model corresponding
      to each component of \code{resp}.  See details below.}
    \item{root}{All ones.  Root variables for graphical model.}
    \item{id}{Categorical.  Indicates individual plants.}
  }
}
\details{
The levels of \code{varb} indicate nodes of the graphical model to which
the corresponding elements of the response vector \code{resp} belong.
This is the typical \dQuote{long} format produced by the R \code{reshape}
function.  For each individual, there are several response variables.
All response variables are combined in one vector \code{resp}.
The variable \code{varb} indicates which \dQuote{original} variable
the number was for.  The variable \code{id} indicates which individual
the number was for.  The levels of \code{varb}, which are the names
of the \dQuote{original} variables are the following.
\code{S1} through \code{S13} are Bernoulli: one if alive, zero if dead.
\code{B2} through \code{B9} are conditionally Poisson: the number of
offspring in the corresponding time period.  Some variables in the
original data that were zero have been deleted.
}
\references{
These data were published in the following, where they were analyzed
by non-aster methods.

Lenski, R. E. and Service, P. M. (1982).
The statistical analysis of population growth rates calculated
    from schedules of survivorship and fecunidity.
\emph{Ecology}, \bold{63}, 655-662.
\doi{10.2307/1936785}.

These data are reanalyzed by aster methods in the following.

Shaw, R. G., Geyer, C. J., Wagenius, S., Hangelbroek, H. H.,
   and Etterson, J. R. (2008)
Unifying life history analyses for inference of fitness and population growth.
\emph{American Naturalist}, \bold{172}, E35-E47.
\doi{10.1086/588063}.
}
\examples{
data(aphid)
### wide version
aphidw <- reshape(aphid, direction = "wide", timevar = "varb",
    v.names = "resp", varying = list(levels(aphid$varb)))
}
\keyword{datasets}

