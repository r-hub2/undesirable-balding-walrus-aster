\name{oats}
\docType{data}
\alias{oats}
\title{Life History Data on Avena barbata}
\description{
  Data on life history traits for the 
  invasive California wild oat \emph{Avena barbata}
}
\usage{data(oats)}
\format{
  A data frame with records for 821 plants.
  Data are already in \dQuote{long} format; no need to reshape.
   \describe{
    \item{resp}{Response vector.}
    \item{varb}{Categorical.  Gives node of graphical model corresponding
      to each component of \code{resp}.  See details below.}
    \item{root}{All ones.  Root variables for graphical model.}
    \item{id}{Categorical.  Indicates individual plants.}
    \item{Plant.id}{Categorical.  Another indicator of individual plants.}
    \item{Env}{Categorical.  Environment in which plant was grown, a
      combination of experimental site and year.}
    \item{Gen}{Categorical.  Ecotype of plant: mesic (M) or xeric (X).}
    \item{Fam}{Categorical.  Accession, nested within ecotype.}
    \item{Site}{Categorical.  Experiment site.  Two sites in these data.}
    \item{Year}{Categorical.  Year in which data were collected.
        Four years in these data.}
    \item{fit}{Indicator (zero or one).  Shorthand
      for \code{as.numeric(oats$varb == "Spike")}.  So-called because the
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
  \item{Surv}{Indicator (zero or one).  Bernoulli, One if individual
    survived to produce flowers.}
  \item{Spike}{Integer.  Zero-truncated Poisson, number of spikelets
    (compound floral structures) observed.}
}
Graphical model is
\deqn{1 \longrightarrow \mbox{Surv} \longrightarrow \mbox{Spike}}{1 -> Surv -> Spike.}
}
\source{
Robert Latta
\url{https://www.dal.ca/faculty/science/biology/faculty-staff/our-faculty/robert-latta.html}
}
\references{
These data are a subset of data previously analyzed using non-aster methods
in the following.

Latta, R. G. (2009).
Testing for local adaptation in \emph{Avena barbata},
a classic example of ecotypic divergence.
\emph{Molecular Ecology}, \bold{18}, 3781--3791.
\doi{https://doi.org/10.1111/j.1365-294X.2009.04302.x}.

These data are reanalyzed by aster methods in the following.

Geyer, C. J., Ridley, C. E., Latta, R. G., Etterson, J. R.,
    and Shaw, R. G. (2013)
Local Adaptation and Genetic Effects on Fitness: Calculations
    for Exponential Family Models with Random Effects.
\emph{Annals of Applied Statistics}, \bold{7}, 1778--1795.
\doi{10.1214/13-AOAS653}.
}
\examples{
data(oats)
}
\keyword{datasets}

