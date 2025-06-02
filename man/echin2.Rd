
\name{echin2}
\docType{data}
\alias{echin2}
\title{Life History Data on Echinacea angustifolia}
\description{
  Data on life history traits for the narrow-leaved purple coneflower
  \emph{Echinacea angustifolia}
}
\usage{data(echin2)}
\format{
  A data frame with records for 557 plants observed over five years.
  Data are already in \dQuote{long} format; no need to reshape.
   \describe{
    \item{resp}{Response vector.}
    \item{varb}{Categorical.  Gives node of graphical model corresponding
        to each component of \code{resp}.  See details below.}
    \item{root}{All ones.  Root variables for graphical model.}
    \item{id}{Categorical.  Indicates individual plants.}
    \item{flat}{Categorical.  Position in growth chamber.}
    \item{row}{Categorical.  Row in the field.}
    \item{posi}{Numerical.  Position within row in the field.}
    \item{crosstype}{Categorical.  See details.}
    \item{yearcross}{Categorical.  Year in which cross was done.}
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
  \item{lds1}{Survival for the first month in the growth
    chamber.}
  \item{lds2}{Ditto for 2nd month in the growth chamber.}
  \item{lds3}{Ditto for 3rd month in the growth chamber.}
  \item{ld01}{Survival for first year in the field.}
  \item{ld02}{Ditto for 2nd year in the field.}
  \item{ld03}{Ditto for 3rd year in the field.}
  \item{ld04}{Ditto for 4th year in the field.}
  \item{ld05}{Ditto for 5th year in the field.}
  \item{roct2003}{Rosette count, measure of size and vigor,
    recorded for 3rd year in the field.}
  \item{roct2004}{Ditto for 4th year in the field.}
  \item{roct2005}{Ditto for 5th year in the field.}
}

These data are complicated by the experiment being done in two parts.
Plants start their life indoors in a growth chamber.  The predictor
variable \code{flat} only makes sense during this time in which three
response variables \code{lds1}, \code{lds2}, and \code{lds3} are observed.
After three months in the growth chamber, the plants (if they
survived, i. e., if \code{lds3 == 1}) were planted in an experimental
field plot outdoors.  The variables \code{row} and \code{posi} only make
sense during this time in which all of the rest of the response variables
are observed.  Because of certain predictor variables only making sense
with respect to certain components of the response vector, the R formula
mini-language is unable to cope, and model matrices must be constructed
"by hand".

\emph{Echinacea angustifolia} is native to North American tallgrass prairie,
which was once extensive but now exists only in isolated remnants.
To evaluate the effects of different mating regimes on the fitness of
resulting progeny, crosses were conducted to produce progeny of (a) mates
from different remnants, (b) mates chosen at random from the same remnant,
and (c) mates known to share maternal parent.  These three categories are
the three levels of \code{crosstype}.
}
\source{
Stuart Wagenius,
\url{https://www.chicagobotanic.org/research/staff/wagenius}
}
\references{
These data are analyzed in the following.                                     

Shaw, R. G., Geyer, C. J., Wagenius, S., Hangelbroek, H. H.,
   and Etterson, J. R. (2008)
Unifying life history analyses for inference of fitness and population growth.
\emph{American Naturalist}, \bold{172}, E35-E47.
\doi{10.1086/588063}.
}
\examples{
data(echin2)
}
\keyword{datasets}

