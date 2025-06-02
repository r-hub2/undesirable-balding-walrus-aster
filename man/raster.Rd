\name{raster}
\alias{raster}
\title{Aster Model Simulation}
\description{
  Random generation of data for Aster models.
}
\usage{
raster(theta, pred, fam, root, famlist = fam.default())
}
\arguments{
  \item{theta}{canonical parameter of the conditional model.
      A matrix, rows represent individuals and columns represent
      nodes in the graphical model.}
  \item{pred}{integer vector of length \code{ncol(theta)} determining
      the graph.  \code{pred[j]} is the index of the predecessor of
      the node with index \code{j} unless the predecessor is a root
      node, in which case \code{pred[j] == 0}.}
  \item{fam}{integer vector of length \code{ncol(theta)} determining
      the exponential family structure of the aster model.  Each element
      is an index into the vector of family specifications given by
      the argument \code{famlist}.}
  \item{root}{A matrix of the same dimensions as \code{theta}.
      Data \code{root[i, j]} is the data for the founder that is
      the predecessor of the \code{[i, j]} node.}
  \item{famlist}{a list of family specifications (see \code{\link{families}}).}
}
\value{
  A matrix of the same dimensions as \code{theta}.  The random data
  for an aster model with the specified graph, parameters, and root
  data.
}
\examples{
### see package vignette for explanation ###
data(echinacea)
vars <- c("ld02", "ld03", "ld04", "fl02", "fl03", "fl04",
    "hdct02", "hdct03", "hdct04")
redata <- reshape(echinacea, varying = list(vars),
     direction = "long", timevar = "varb", times = as.factor(vars),
     v.names = "resp")
redata <- data.frame(redata, root = 1)
pred <- c(0, 1, 2, 1, 2, 3, 4, 5, 6)
fam <- c(1, 1, 1, 1, 1, 1, 3, 3, 3)
hdct <- grep("hdct", as.character(redata$varb))
hdct <- is.element(seq(along = redata$varb), hdct)
redata <- data.frame(redata, hdct = as.integer(hdct))
aout4 <- aster(resp ~ varb + nsloc + ewloc + pop * hdct - pop,
    pred, fam, varb, id, root, data = redata)
newdata <- data.frame(pop = levels(echinacea$pop))
for (v in vars)
    newdata[[v]] <- 1
newdata$root <- 1
newdata$ewloc <- 0
newdata$nsloc <- 0
renewdata <- reshape(newdata, varying = list(vars),
    direction = "long", timevar = "varb", times = as.factor(vars),
    v.names = "resp")
hdct <- grep("hdct", as.character(renewdata$varb))
hdct <- is.element(seq(along = renewdata$varb), hdct)
renewdata <- data.frame(renewdata, hdct = as.integer(hdct))
beta.hat <- aout4$coef
theta.hat <- predict(aout4, model.type = "cond", parm.type = "canon")
theta.hat <- matrix(theta.hat, nrow = nrow(aout4$x), ncol = ncol(aout4$x))
xstar <- raster(theta.hat, pred, fam, aout4$root)
aout4star <- aster(xstar, aout4$root, pred, fam, aout4$modmat, beta.hat)
beta.star <- aout4star$coef
print(cbind(beta.hat, beta.star))
}
\seealso{\code{\link{aster}}}
\keyword{distribution}
