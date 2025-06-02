\name{predict.aster}
\alias{predict.aster}
\alias{predict.aster.formula}
\concept{regression}
\concept{exponential family}
\concept{graphical model}
\title{Predict Method for Aster Model Fits}
\usage{
\method{predict}{aster}(object, x, root, modmat, amat,
    parm.type = c("mean.value", "canonical"),
    model.type = c("unconditional", "conditional"),
    is.always.parameter = FALSE,
    se.fit = FALSE, info = c("expected", "observed"),
    info.tol = sqrt(.Machine$double.eps), newcoef = NULL,
    gradient = se.fit, \ldots)

\method{predict}{aster.formula}(object, newdata, varvar, idvar, root, amat,
    parm.type = c("mean.value", "canonical"),
    model.type = c("unconditional", "conditional"),
    is.always.parameter = FALSE,
    se.fit = FALSE, info = c("expected", "observed"),
    info.tol = sqrt(.Machine$double.eps), newcoef = NULL,
    gradient = se.fit, \ldots)
}
\arguments{
  \item{object}{a fitted object of class inheriting from \code{"aster"}
    or \code{"aster.formula"}.}

  \item{modmat}{a model matrix to use instead of \code{object$modmat}.
    Must have the same structure (three-dimensional array, first index runs
    over individuals, second over nodes of the graphical model, third over
    covariates).  Must have the same second and third dimensions as
    \code{object$modmat}.  The second and third components of
    \code{dimnames(modmat)} and \code{dimnames(object$modmat)} must also be
    the same.

    May be missing, in which case \code{object$modmat} is used.

    \code{predict.aster.formula} constructs such a \code{modmat} from
    \code{object$formula}, the data frame \code{newdata}, and variables
    in the environment of the formula.  When \code{newdata} is missing,
    \code{object$modmat} is used.}

  \item{x}{response.  Ignored and may be missing unless
    \code{parm.type = "mean.value"} and \code{model.type = "conditional"}.
    Even then may be missing when \code{modmat} is missing,
    in which case \code{object$x} is used.  A matrix whose first and
    second dimensions and the corresponding dimnames agrees with
    those of \code{modmat} and \code{object$modmat}.

    \code{predict.aster.formula} constructs such an \code{x} from
    the response variable name in \code{object$formula},
    the data frame \code{newdata},
    and the variables in the environment of the formula.  When \code{newdata}
    is missing, \code{object$x} is used.}

  \item{root}{root data.  Ignored and may be missing unless
    \code{parm.type == "mean.value"}.
    Even then may be missing when \code{modmat} is missing,
    in which case \code{object$root} is used.  A matrix of the
    same form as \code{x}.

    \code{predict.aster.formula} looks up the variable supplied as
    the argument \code{root} in the data frame \code{newdata} or in
    the variables in the environment of the formula and makes it a matrix
    of the same form as \code{x}.  When \code{newdata}
    is missing, \code{object$root} is used.}

  \item{amat}{if \code{zeta} is the requested prediction (mean value
    or canonical, unconditional or conditional, depending on \code{parm.type}
    and \code{model.type}), then we predict the linear function
    \code{t(amat) \%*\% zeta}.  May be missing, in which case the identity
    linear function is used.

    For \code{predict.aster}, a three-dimensional
    array with \code{dim(amat)[1:2] == dim(modmat)[1:2]}.

    For \code{predict.aster.formula}, a three-dimensional array
    of the same dimensions as required for \code{predict.aster}
    (even though \code{modmat} is not provided).  First dimension
    is number of individuals in \code{newdata}, if provided, otherwise
    number of individuals in \code{object$data}.  Second dimension
    is number of variables (\code{length(object$pred)}).
  }

  \item{parm.type}{the type of parameter to predict.  The default is
    mean value parameters (the opposite of the default
    for \code{\link{predict.glm}}), the expected value of a linear function
    of the response under the MLE probability model (also called the
    MLE of the mean value parameter).  The expectation is unconditional
    or conditional depending on \code{parm.type}.

    The alternative \code{"canonical"} is the value of a linear function
    of the MLE of canonical parameters under the MLE probability model.
    The canonical parameter is unconditional
    or conditional depending on \code{parm.type}.

    The value of this argument can be abbreviated.
  }

  \item{model.type}{the type of model in which to predict.  The default is
    \code{"unconditional"} in which case the parameters (either mean value
    or canonical, depending on the value of \code{parm.type}) are those
    of an unconditional model.
    The alternative is \code{"conditional"} in which case the parameters
    are those of a conditional model.

    The value of this argument can be abbreviated.
  }

  \item{is.always.parameter}{logical, default \code{FALSE}.
      Only affects the result when \code{parm.type = "mean.value"} and
      \code{model.type = "conditional"}.  \code{TRUE} means the conditional
      mean value parameter is produced.  \code{FALSE} means the conditional
      mean values themselves are produced (which depend on data so are not
      parameters).  See Conditional Mean Values Section below for further
      explanation.}

  \item{se.fit}{logical switch indicating if standard errors are required.}

  \item{info}{the type of Fisher information to use to compute standard errors.}

  \item{info.tol}{tolerance for eigenvalues of Fisher information.
    If \code{eval} is the vector of eigenvalues of the information matrix,
    then \code{eval < cond.tol * max(eval)} are considered zero.  Hence the
    corresponding eigenvectors are directions of constancy or recession of
    the log likelihood.}

  \item{newdata}{optionally, a data frame in which to look for variables with
    which to predict.  If omitted, see \code{modmat} above.  See also details
    section below.}

  \item{varvar}{a variable of length \code{nrow(newdata)}, typically a
    variable in \code{newdata}
    that is a factor whose levels are character strings
    treated as variable names.  The number of variable names is \code{nnode}.
    Must be of the form \code{rep(vars, each = nind)} where \code{vars} is
    a vector of variable names.  Not used if \code{newdata} is missing.}

  \item{idvar}{a variable of length \code{nrow(newdata)}, typically a
    variable in \code{newdata}
    that indexes individuals.  The number
    of individuals is \code{nind}.
    Must be of the form \code{rep(inds, times = nnode)} where \code{inds} is
    a vector of labels for individuals.  Not used if \code{newdata} is missing.}

  \item{newcoef}{if not \code{NULL},
    a variable of length \code{object$coefficients} and used
    in its place when one wants predictions at other than the fitted
    coefficient values.}

  \item{gradient}{if \code{TRUE} return the gradient (Jacobian of the
      transformation) matrix.  This matrix has number of rows equal to the
      length of the fitted values and number of columns equal to the number
      of regression coefficients.  It is the derivative matrix (matrix of
      partial derivatives) of the mapping from regression coefficients to
      whatever the predicted values are, which depends on what the
      arguments \code{newdata}, \code{amat}, \code{parm.type}, and
      \code{model.type} are.}

  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  Obtains predictions (parameter estimates) and optionally estimates
  standard errors of those predictions (parameter estimates) from
  a fitted Aster model object.
}
\section{Conditional Mean Values}{
Both the original aster paper (Geyer, et al., 2007) and this package are
weird about conditional mean values.  Equation (10) of that paper defines
(using different notation from what we use here)
\deqn{\xi_j = E(x_j | x_{p(j)})}{xi[j] = E(x[j] | x[p(j)])}
where \eqn{x_j}{x[j]} are components of the response vector and \eqn{p(j)}
denotes denotes the predecessor of node \eqn{j}.  That paper explicitly says
that this is is not a parameter because it depends on the data.  In fact
\deqn{E(x_j | x_{p(j)}) = x_{p(j)} E(x_j | x_{p(j)} = 1)}{E(x[j] | x[p(j)]) = x[p(j)] E(x[j] | x[p(j)] = 1)}
(this is equation (3) of that paper in different notation).
Thus it is weird to use a Greek letter to denote this.

There should be a conditional mean value parameter, and
Geyer (2010, equation (11b)) defines it as
\deqn{\xi_j = E(y_j | y_{p(j)} = 1)}{xi[j] = E(y[j] | y[p(j)] = 1)}
(This equation only makes sense when the conditioning event
\eqn{x_{p(j)} = 1}{x[p(j)] = 1} is possible, which it is not for
\eqn{k}-truncated arrows for \eqn{k > 0}.  Then a more complicated definition
must be used.  By definition \eqn{x_j}{x[j]} is the sum
of \eqn{x_{p(j)}}{x[p(j)]} independent and identically distributed random
variables, and \eqn{\xi_j}{xi[j]} is always the mean of one of those
random variables.)
This gives us the important relationship between conditional and unconditional
mean value parameters
\deqn{\mu_j = \xi_j \mu_{p(j)}}{mu[j] = xi[j] mu[p(j)]}
which holds for all successor nodes \eqn{j}.
All later writings of this author use this definition of \eqn{\xi}{xi}
as does the R package \code{aster2} (Geyer, 2017).
This is one of six important parameterizations of an unconditional aster model
(Geyer, 2010, Sections 2.7 and 2.8).  The R package \code{aster2} uses all
of them.

This function (as of version 1.0 of this package) has an argument
\code{is.always.parameter} to switch between these two definitions
in case \code{parm.type = "mean.value"} and \code{model.type = "conditional"}
are specified.  Then \code{is.always.parameter = TRUE} specifies that the latter
definition of \eqn{\xi}{xi} is produced (which is a parameter, with all other
options for \code{parm.type} and \code{model.type}).  The option
\code{is.always.parameter = FALSE} specifies that the former
definition of \eqn{\xi}{xi} is produced (which is a conditional expectation
but not a parameter) and is what this function produced in versions of this
package before 1.0.
}
\details{
Note that \code{model.type} need have nothing to do with the type
of the fitted aster model, which is \code{object$type}.

Whether the
fitted model is conditional or unconditional, one typically wants
\emph{unconditional} mean value parameters, because conditional mean
value parameters for hypothetical individuals depend on the hypothetical
data \code{x}, which usually makes no scientific sense.

If one asks for \emph{conditional} mean value parameters, one gets
them only if \code{is.always.parameter = TRUE} is specified.
Otherwise, conditional expectations that are not parameters (because
they depend on data) are produced.
See Conditional Mean Values Section for more about this.

Similarly, if \code{object$type == "conditional"}, then the conditional
canonical parameters are a linear function of the regression coefficients
\eqn{\theta = M \beta}{theta = M beta}, where \eqn{M} is the model matrix,
but one can predict either \eqn{\theta}{theta} or the unconditional
canonical parameters \eqn{\varphi}{phi},
as selected by \code{model.type}.
Similarly, if \code{object$type == "unconditional"},
so \eqn{\varphi = M \beta}{phi = M beta}, one can predict either
\eqn{\theta}{theta} or \eqn{\varphi}{phi}
as selected by \code{model.type}.

The specification of the prediction model is confusing because there
are so many possibilities.  First the \dQuote{usual} case.
The fit was done using a formula, found in \code{object$formula}.
A data frame \code{newdata} that has the same variables as \code{object$data},
the data frame used in the fit, but may have different rows (representing
hypothetical individuals) is supplied.
But \code{newdata} must specify \emph{all nodes}
of the graphical model for each (hypothetical, new) individual,
just like \code{object$data} did for real observed individuals.
Hence \code{newdata} is typically constructed using \code{\link{reshape}}.
See also the details section of \code{\link{aster}}.

In this \dQuote{usual} case we need \code{varvar} and \code{idvar} to
tell us what rows of \code{newdata} correspond to which individuals and
nodes (the same role they played in the original fit by \code{\link{aster}}).
If we are predicting canonical parameters, then we do not need \code{root} or
\code{x}.
If we are predicting unconditional mean value parameters, then
we also need \code{root} but not \code{x}.
If we are predicting conditional mean value parameters, then
we also need both \code{root} and \code{x}.
In the \dQuote{usual} case, these are found in \code{newdata} and
not supplied as arguments to \code{predict}.  Moreover, \code{x}
is not named \code{"x"} but is the response in \code{out$formula}.

The next case, \code{predict(object)} with no other arguments,
is often used with linear models (\code{\link{predict.lm}}),
but we expect will be little used for aster models.  As for linear
models, this \dQuote{predicts} the observed data.  In this case
\code{modmat}, \code{x}, and \code{root} are found in \code{object}
and nothing is supplied as an argument to \code{predict.aster}, except
perhaps \code{amat} if one wants a function of predictions for the observed
data.

The final case, also perhaps little used, is a fail-safe mode for problems
in which the R formula language just cannot be bludgeoned into doing what
you want.  This is the same reason \code{\link{aster.default}} exists.
Then a model matrix can be constructed \dQuote{by hand}, and the function
\code{predict.aster} is used instead of \code{predict.aster.formula}.

Note that it is possible to use a \dQuote{constructed by hand}
model matrix even if \code{object} was produced by \code{\link{aster.formula}}.
Simply explicitly call \code{predict.aster} rather than \code{predict}
to override the R method dispatch (which would
call \code{predict.aster.formula} in this case).
}
\value{
  If \code{se.fit = FALSE} and \code{gradient = FALSE}, a vector of predictions.
  If \code{se.fit = TRUE}, a list with components
  \item{fit}{Predictions}
  \item{se.fit}{Estimated standard errors}
  \item{gradient}{The gradient of the transformation from
    regression coefficients to predictions}
  If \code{gradient = TRUE}, a list with components
  \item{fit}{Predictions}
  \item{gradient}{The gradient of the transformation from
    regression coefficients to predictions}
}
\references{
Geyer, C. J. (2010)
A Philosophical Look at Aster Models.
Technical Report No. 676.  School of Statistics, University of Minnesota.
\url{https://purl.umn.edu/57163}.

Geyer, C.~J. (2017).
R package \code{aster2} (Aster Models), version 0.3.
\url{https://cran.r-project.org/package=aster2}.

Geyer, C. J., Wagenius, S., and Shaw, R. G. (2007)
Aster models for life history analysis.
\emph{Biometrika}, \bold{94}, 415--426.
\doi{10.1093/biomet/asm030}.
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
aout <- aster(resp ~ varb + level : (nsloc + ewloc) + hdct : pop,
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
hdct <- grepl("hdct", as.character(renewdata$varb))
renewdata <- data.frame(renewdata, hdct = as.integer(hdct))
level <- gsub("[0-9]", "", as.character(renewdata$varb))
renewdata <- data.frame(renewdata, level = as.factor(level))
nind <- nrow(newdata)
nnode <- length(vars)
amat <- array(0, c(nind, nnode, nind))
for (i in 1:nind)
    amat[i , grep("hdct", vars), i] <- 1
foo <- predict(aout, varvar = varb, idvar = id, root = root,
    newdata = renewdata, se.fit = TRUE, amat = amat)
bar <- cbind(foo$fit, foo$se.fit)
dimnames(bar) <- list(as.character(newdata$pop), c("Estimate", "Std. Error"))
print(bar)
}
\keyword{models}
\keyword{regression}
