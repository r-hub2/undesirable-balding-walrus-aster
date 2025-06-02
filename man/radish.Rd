\name{radish}
\docType{data}
\alias{radish}
\alias{radish2}
\title{Life History Data on Raphanus sativus}
\description{
  Data on life history traits for the 
  invasive California wild radish \emph{Raphanus sativus}
}
\usage{
  data(radish)
  data(radish2)
}
\format{
  A data frame with records for 286 plants.
  Data are already in \dQuote{long} format; no need to reshape.
   \describe{
    \item{resp}{Response vector.}
    \item{varb}{Categorical.  Gives node of graphical model corresponding
      to each component of \code{resp}.  See details below.}
    \item{root}{All ones.  Root variables for graphical model.}
    \item{id}{Categorical.  Indicates individual plants.}
    \item{Site}{Categorical.  Experimental site where plant was grown.
      Two sites in this dataset.}
    \item{Block}{Categorical.  Block nested within site.}
    \item{Region}{Categorical.  Region from which individuals
      were obtained: northern, coastal California (N) or southern, inland
      California (S).}
    \item{Pop}{Categorical.  Wild population nested within region.}
    \item{varbFlowering}{Indicator (zero or one).  Shorthand
      for \code{as.numeric(radish$varb == "Flowering")}.}
    \item{varbFlowers}{Indicator (zero or one).  Shorthand
      for \code{as.numeric(radish$varb == "Flowers")}.}
    \item{fit}{Indicator (zero or one).  Shorthand
      for \code{as.numeric(radish$varb == "Fruits")}.  So-called because the
      components of \code{outcome} indicated are the best surrogate of
      Darwinian fitness in these data.}
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
  \item{Flowering}{Indicator (zero or one).  Bernoulli, One if individual
    survived to produce flowers.}
  \item{Flowers}{Integer.  Zero-truncated Poisson, number of flowers observed.}
  \item{Fruits}{Integer.  Poisson, number of fruits observed.}
}
Graphical model is
\deqn{1 \longrightarrow \mbox{Flowering} \longrightarrow \mbox{Flowers} \longrightarrow \mbox{Fruits}}{1 -> Flowering -> Flowers -> Fruits.}

Dataset \code{radish2} has exactly the same structure but the estimated value
of the \code{pop} variance component is zero.  These data were simulated by
R function \code{\link{raster}} in this package.
}
\source{
Caroline Ridley
}
\references{
These data are a subset of data previously analyzed using fixed effect
aster methods (R function \code{aster}) in the following.

Ridley, C. E. and Ellstrand, N. C. (2010).
Rapid evolution of morphology and adaptive life history in
    the invasive California wild radish (\emph{Raphanus sativus}) and
    the implications for management.
\emph{Evolutionary Applications}, \bold{3}, 64--76.

These data are a subset of data previously analyzed using random effect
aster methods (R function \code{reaster}) in the following.

Geyer, C. J., Ridley, C. E., Latta, R. G., Etterson, J. R.,
    and Shaw, R. G. (2013)
Local Adaptation and Genetic Effects on Fitness: Calculations
    for Exponential Family Models with Random Effects.
\emph{Annals of Applied Statistics}, \bold{7}, 1778--1795.
\doi{10.1214/13-AOAS653}.
}
\seealso{
\code{\link{pickle}}
}
\examples{
data(radish)
}
\keyword{datasets}
