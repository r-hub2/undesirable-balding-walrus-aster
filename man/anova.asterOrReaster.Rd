\name{anova.asterOrReaster}
\alias{anova.asterOrReaster}
\alias{anovaAsterOrReasterList}
\alias{anova.aster}
\alias{anova.reaster}
\title{Analysis of Deviance for Reaster Model Fits}
\usage{
\method{anova}{asterOrReaster}(object, \dots,
    tolerance = .Machine$double.eps ^ 0.75)
anovaAsterOrReasterList(objectlist, tolerance = .Machine$double.eps ^ 0.75)
}
\description{
  Compute an analysis of deviance table for two or more aster model
  fits with or without random effects.
}
\arguments{
  \item{object, \dots}{objects of class \code{"asterOrReaster"}, typically
    the result of a call to \code{\link{aster}} or \code{\link{reaster}}.}
  \item{objectlist}{list of objects of class \code{"asterOrReaster"}.}
  \item{tolerance}{tolerance for comparing nesting of model matrices.}
}
\details{
  Constructs a table having a row for the
  degrees of freedom and deviance for each model.
  For all but the first model, the change in degrees of freedom and
  deviance is also given, as is the corresponding asymptotic \eqn{P}-value.

  For objects of class \code{"reaster"}, the
  quantity called deviance is only approximate.  See references on
  help for \code{\link{reaster}}.

  When objects of class \code{"reaster"} are among those supplied,
  degrees of freedom for fixed effects and degrees of freedom for
  variance components are reported separately, because tests for fixed
  effects are effectively two-tailed and tests for variance components
  are effectively one-tailed.

  In case models being compared differ by one variance component, the
  reference distribution is half a chi-square with the fixed effect
  degrees of freedom (difference of number of fixed effects in the two
  models) and half a chi-square with one more degrees of freedom.

  In case models being compared differ by two or more variance components,
  we do not know how to how to do the test.  The reference distribution is
  a mixture of chi-squares but the mixing weights are difficult to calculate.
  An error is given in this case.
}
\value{
  An object of class \code{"anova"} inheriting from class \code{"data.frame"}.
}
\section{Warning}{
  The comparison between two or more models by \code{anova} or
  \code{anovaAsterOrReasterList} will only be valid if they
  are (1) fitted to the same dataset, 
  (2) models are nested,
  (3) have the same
  dependence graph and exponential families.
  Some of this is currently checked.  Some warnings are given.
}
\seealso{
  \code{\link{aster}}, \code{\link{reaster}}, \code{\link{anova}}.
}
\examples{
### see package vignette for explanation ###
library(aster)
data(echinacea)
vars <- c("ld02", "ld03", "ld04", "fl02", "fl03", "fl04",
    "hdct02", "hdct03", "hdct04")
redata <- reshape(echinacea, varying = list(vars), direction = "long",
    timevar = "varb", times = as.factor(vars), v.names = "resp")
redata <- data.frame(redata, root = 1)
pred <- c(0, 1, 2, 1, 2, 3, 4, 5, 6)
fam <- c(1, 1, 1, 1, 1, 1, 3, 3, 3)
hdct <- grepl("hdct", as.character(redata$varb))
redata <- data.frame(redata, hdct = as.integer(hdct))
level <- gsub("[0-9]", "", as.character(redata$varb))
redata <- data.frame(redata, level = as.factor(level))
aout1 <- aster(resp ~ varb + hdct : (nsloc + ewloc + pop),
    pred, fam, varb, id, root, data = redata)
aout2 <- aster(resp ~ varb + level : (nsloc + ewloc) + hdct : pop,
    pred, fam, varb, id, root, data = redata)
aout3 <- aster(resp ~ varb + level : (nsloc + ewloc + pop),
    pred, fam, varb, id, root, data = redata)
anova(aout1, aout2, aout3)

# now random effects models and models without random effects mixed
\dontrun{
### CRAN policy says examples must take < 5 sec.
### This doesn't (on their computers).
data(radish)
pred <- c(0,1,2)
fam <- c(1,3,2)
rout2 <- reaster(resp ~ varb + fit : (Site * Region),
    list(block = ~ 0 + fit : Block, pop = ~ 0 + fit : Pop),
    pred, fam, varb, id, root, data = radish)
rout1 <- reaster(resp ~ varb + fit : (Site * Region),
    list(block = ~ 0 + fit : Block),
    pred, fam, varb, id, root, data = radish)
rout0 <- aster(resp ~ varb + fit : (Site * Region),
    pred, fam, varb, id, root, data = radish)
anova(rout0, rout1, rout2)
}
}
\keyword{models}
\keyword{regression}

