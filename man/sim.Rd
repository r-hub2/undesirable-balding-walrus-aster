\name{sim}
\docType{data}
\alias{sim}
\alias{beta.true}
\alias{mu.true}
\alias{phi.true}
\alias{theta.true}
\alias{fam}
\alias{pred}
\alias{vars}
\alias{ladata}
\alias{redata}
\title{Simulated Life History Data}
\description{
  Data on life history traits for four years and five fitness components
}
\usage{data(sim)}
\format{
  Loads nine objects.
  The objects \code{beta.true},  \code{mu.true}, \code{phi.true}, and
  \code{theta.true} are the simulation truth parameter values in
  different parametrizations.
    \describe{
      \item{beta.true}{Regression coefficient vector for model
      \code{resp ~ varb + 0 + z1 + z2 + I(z1^2) + I(z1*z2) + I(z2^2)}.}
      \item{mu.true}{Unconditional mean value parameter vector for same
      model.} 
      \item{phi.true}{Unconditional canonical value parameter vector for
      same model.} 
      \item{theta.true}{Conditional canonical value parameter vector for
      same model.}
    }
  The objects \code{fam},  \code{pred}, and \code{vars}
  specify the aster model graphical and probabilistic structure.
    \describe{
      \item{fam}{Integer vector giving the families of the variables in
      the graph.}
      \item{pred}{Integer vector giving the predecessors of the variables in
      the graph.}
      \item{vars}{Character vector giving the names of the variables in
      the graph.}
    }
  The objects \code{ladata} and \code{redata} are the simulated data
  in two forms \code{"wide"} and \code{"long"} in the terminology
  of the \code{reshape} function.
    \describe{
      \item{ladata}{Data frame with variables \code{y}, \code{z1},
      \code{z2} used for Lande-Arnold type estimation of fitness landscape.
      \code{y} is the response, fitness, and \code{z1} and \code{z1} are
      predictor variables, phenotypes.}
      \item{redata}{Data frame with variables \code{resp}, \code{z1},
      \code{z2}, \code{varb}, \code{id}, \code{root}
      used for aster type estimation of fitness landscape.
      \code{resp} is the response, containing all components of fitness,
      and \code{z1} and \code{z1} are predictor variables, phenotypes.
      \code{varb} is a factor whose levels are are elements of \code{vars}
      indicating which elements of \code{resp} go with which nodes of the
      aster model graphical structure.  The variables \code{z1} and \code{z2}
      have been set equal to zero except when \code{grep("nseed", varb)} is
      \code{TRUE}.  For the rationale see Section 3.2 of TR 669 referenced
      below.
    }
  }
}
\source{
Geyer, C. J and Shaw, R. G. (2008)
Supporting Data Analysis for a talk to be given at Evolution 2008.
Technical Report No. 669.  School of Statistics, University of Minnesota.
\url{https://hdl.handle.net/11299/56204}.
}
\references{
Geyer, C. J and Shaw, R. G. (2009)
Hypothesis Tests and Confidence Intervals
    Involving Fitness Landscapes fit by Aster Models.
Technical Report No. 671.  School of Statistics, University of Minnesota.
\url{https://hdl.handle.net/11299/56219}.
}
\examples{
data(sim)
\dontrun{
### CRAN policy says examples must take < 5 sec.  This doesn't.
out6 <- aster(resp ~ varb + 0 + z1 + z2 + I(z1^2) + I(z1*z2) + I(z2^2),
    pred, fam, varb, id, root, data = redata)
summary(out6)
}
lout <- lm(y ~ z1 + z2 + I(z1^2) + I(z1*z2) + I(z2^2), data = ladata)
summary(lout)
}
\keyword{datasets}

