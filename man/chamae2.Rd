\name{chamae2}
\docType{data}
\alias{chamae2}
\title{Life History Data on Chamaecrista fasciculata}
\description{
  Data on life history traits for the partridge pea
  \emph{Chamaecrista fasciculata}
}
\usage{data(chamae2)}
\format{
  A data frame with records for 2239 plants.
  Data are already in \dQuote{long} format; no need to reshape.
   \describe{
    \item{resp}{Response vector.}
    \item{varb}{Categorical.  Gives node of graphical model corresponding
      to each component of \code{resp}.  See details below.}
    \item{root}{All ones.  Root variables for graphical model.}
    \item{id}{Categorical.  Indicates individual plants.}
    \item{STG1N}{Numerical.  Reproductive stage.  Integer with only 3 values
      in this dataset.}
    \item{LOGLVS}{Numerical.  Log leaf number.}
    \item{LOGSLA}{Numerical.  Log leaf thickness.}
    \item{BLK}{Categorical.  Block within experiment.}
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
of the \dQuote{original} variables are
\describe{
  \item{fecund}{Fecundity.  Bernoulli, One if any fruit, zero if no fruit.}
  \item{fruit}{Integer.  Number of fruits observed.}
}
}
\source{
Julie Etterson
\url{https://sites.google.com/d.umn.edu/dr-julie-r-etterson/home}
}
\references{
These data are a subset of data previously analyzed by non-aster methods in the following.

Etterson, J. R. (2004).
Evolutionary potential of \emph{Chamaecrista fasciculata} in
    relation to climate change.  I. Clinal patterns of selection along
    an environmental gradient in the great plains.
\emph{Evolution}, \bold{58}, 1446-1458.
\doi{10.1111/j.0014-3820.2004.tb01726.x}.

Etterson, J. R., and Shaw, R. G. (2001).
Constraint to adaptive evolution in response to global warming.
\emph{Science}, \bold{294}, 151-154.
\doi{10.1126/science.1063656}.

These data are reanalyzed by aster methods in the following.

Shaw, R. G., Geyer, C. J., Wagenius, S., Hangelbroek, H. H.,
   and Etterson, J. R. (2008)
Unifying life history analyses for inference of fitness and population growth.
\emph{American Naturalist}, \bold{172}, E35-E47.
\doi{10.1086/588063}.
}
\examples{
data(chamae2)
### wide version
chamae2w <- reshape(chamae2, direction = "wide", timevar = "varb",
    v.names = "resp", varying = list(levels(chamae2$varb)))
}
\keyword{datasets}

