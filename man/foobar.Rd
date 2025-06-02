\name{foobar}
\docType{data}
\alias{foobar}
\concept{direction of recession}
\title{Toy Life History Data having Directions of Recession}
\description{
    Toy life history data created to exhibit the phenomenon of directions
    of recession.   It was analyzed in a special topics course on aster models
    Fall 2018.
}
\usage{data(foobar)}
\format{
  \code{data(foobar)} loads four R objects.
  \describe{
    \item{fam}{a vector of family indices.}
    \item{pred}{a vector of predecessor indices.}
    \item{vars}{a vector of variable names associated with the nodes
      of the aster graph.}
    \item{redata}{a data frame having data on 300 individuals, each of
      which has \code{length(vars) == 4} components of fitness,
      so the aster graph for one individual has four nodes.  This data
      frame is already in long format; no need to \code{\link{reshape}}.
      The variables in this data frame are:
      \describe{
        \item{resp}{the response vector.}
        \item{varb}{Categorical.  Gives node of graphical model corresponding
          to each component of \code{resp}.  See details below.}
        \item{root}{All ones.  Root variables for graphical model.}
        \item{id}{Categorical.  Individual ID numbers.}
        \item{trt}{Categorical.  Treatment.}
        \item{blk}{Categorical.  Block.}
        \item{fit}{Bernoulli (zero-or-one-valued).  Indicator variable
          of the fitness nodes
          of the graph; in these data there is just one node for fitness.}
      }
    }
  }
}
\details{
The levels of \code{varb} indicate nodes of the graphical model to which
the corresponding elements of the response vector \code{resp} belong.
This is the typical \dQuote{long} format produced by the R \code{reshape}
function.  For each individual, there are several response variables.
All response variables are combined in one vector \code{resp}.
The variable \code{varb} indicates which \dQuote{original} variable
the the corresponding component of the response vector was.
The variable \code{id} indicates which individual the corresponding
component of the response vector was.
}
\source{
Charles J. Geyer
\url{https://www.stat.umn.edu/geyer/8931aster/foobar.rda}
}
\references{
These data are analyzed in deck 9 of the slides for a special topics
course on aster models taught fall semester 2018
(\url{https://www.stat.umn.edu/geyer/8931aster/slides/s9.pdf}).
}
\examples{
data(foobar)
library(aster)
aout <- aster(resp ~ varb + fit : (trt * blk), pred,
    fam, varb, id, root, data = redata)
foo <- try(summary(aout))
# gives an error about directions of recession
# get directions of recession
dor <- attr(foo, "condition")$dor
dor
# found one apparent direction of recession
# from regular pattern
# it looks like a true direction of recession
dor <- dor / max(abs(dor))
dor
# but what does it do?  For that need to map to saturated model
# parameter space
modmat <- aout$modmat
dim(modmat)
# oof!  modmat is three-dimensional.  Need an actual matrix here.
modmat <- matrix(as.vector(modmat), ncol = length(dor))
dor.phi <- drop(modmat \%*\% dor)
names(dor.phi) <- with(redata, paste(id, as.character(varb), sep = "."))
dor.phi[dor.phi != 0]
fam.default()[fam[vars == "seeds"]]
# since the support of the Poisson distribution is bounded above,
# actually this must be minus the DOR (if it is a DOR at all).
# check that all components of response vector for which dor.phi == 1 are zero
# (lower bound of Poisson range)
all(redata$resp[dor.phi == 1] == 0)
# so minus dor.phi is a true direction of recession in the saturated model
# canonical parameter space, and minus dor is a true direction of recession
# in the submodel canonical parameter space
#
# try to get more info
trt.blk <- with(redata,
    paste(as.character(trt), as.character(blk), sep = "."))
unique(trt.blk[dor.phi == 1])
# the reason for the direction of recession is that every individual getting
# treatment a in block A had zero seeds.
#
# the reason the submodel DOR, R object dor, was so hard to interpret was
# because fit:trta:blkA is not in the model.  So let's force it in
redata <- transform(redata, trt = relevel(trt, ref = "b"),
    blk = relevel(blk, ref = "B"))
# Note: following code is copied exactly from above.  Only difference
# is releveling in the immediately preceding statement
aout <- aster(resp ~ varb + fit : (trt * blk), pred,
    fam, varb, id, root, data = redata)
foo <- try(summary(aout))
dor <- attr(foo, "condition")$dor
dor <- dor / max(abs(dor))
dor
# now it is obvious from looking at this dor that all individuals in trt a
# and blk A are at the lower end (zero) of the Poisson range.
# maybe the other dor we had previously would be "obvious" to someone
# sufficiently skilled in understanding the meaning of regression coefficients
# but not "obvious" to your humble author
#
# as for what to do about this, see the course slides cited in the reference
# section.  There is no single Right Thing to do.
}
\keyword{datasets}

