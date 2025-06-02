\name{summary.aster}
\alias{summary.aster}
\alias{print.summary.aster}
\title{Summarizing Aster Model Fits}
\usage{
\method{summary}{aster}(object, info = c("expected", "observed"),
    info.tol = sqrt(.Machine$double.eps),  show.graph = FALSE, \dots)

\method{print}{summary.aster}(x, digits = max(3, getOption("digits") - 3),
      signif.stars = getOption("show.signif.stars"), \dots)
}
\arguments{
  \item{object}{an object of class \code{"aster"}, usually, a result of a
    call to \code{\link{aster}}.}
  \item{info}{the type of Fisher information use to compute standard errors.}
  \item{info.tol}{tolerance for eigenvalues of Fisher information.
    If \code{eval} is the vector of eigenvalues of the information matrix,
    then \code{eval < cond.tol * max(eval)} are considered zero.  Hence the
    corresponding eigenvectors are directions of constancy or recession of
    the log likelihood.}
  \item{show.graph}{if \code{TRUE}, show the graphical model.}
  \item{x}{an object of class \code{"summary.aster"}, usually, a result of a
    call to \code{summary.aster}.}
  \item{digits}{the number of significant digits to use when printing.}
  \item{signif.stars}{logical. If \code{TRUE}, \dQuote{significance stars}
    are printed for each coefficient.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  These functions are all \code{\link{methods}} for class \code{aster} or
  \code{summary.aster} objects.
}
\value{
  \code{summary.aster} returns an object of class \code{"summary.aster"}
  list with the same components as \code{object}, which is of class
  \code{"aster"}.
}
\section{Directions of Recession}{
This function may give an error message
\code{"cannot compute standard errors, apparent directions of recession"}.
There are two reasons why this can happen.
\itemize{
\item There may actually be a direction of recession (DOR).  Then the
    maximum likelihood estimate does not exist; increasing the likelihood
    drives (some of) the coefficients to infinity or minus infinity.

    This function's guess at the DOR can be extracted
    from the error object obtained by wrapping this function
    in \code{\link{try}} and then extracting the \code{dor} component
    of the \code{condition} attribute of the error object.
    An example of this is on the help page for the \code{\link{foobar}}
    data set.

    This function's guessed DOR are apparent null eigenvector(s) of the Fisher
    information matrix.  Due to inaccuracy of computer arithmetic, this is
    only a guess.  What are deemed null eigenvectors is controlled by the
    \code{info.tol} argument of this function.  Reducing \code{info.tol}
    to perhaps \code{1e-9} or \code{1e-10} or even a little lower
    may make the putative DOR go away.  In this case they were probably
    bogus (see next item).  Reducing \code{info.tol} to near or below the
    machine epsilon \code{.Machine$double.eps} (\code{\link{.Machine}})
    instructs this function to feed you garbage with no error or warning.

    Putative DOR are probably true DOR if they are highly patterned with
    many zero or nearly zero components and other components that are
    nearly (small) integer multiples of each other.
    Putative DOR are probably bogus if they look like random noise.

    DOR, if true, cannot simply be ignored.  For more information,
    including how to do more rigorous investigation of whether putative
    DOR are true or bogus,
    see the example on the help page for the \code{\link{foobar}}
    data set and the reference cited on that help page.

\item All of the putative directions of recession may be bogus.
    Due to inaccuracy of computer arithmetic, ill-conditioning
    of predictor variables, or ill-conditioning of the aster model itself
    (large graphs cause more inaccurate computation), what appear to be
    null eigenvectors of the Fisher information matrix need not be true
    null eigenvectors.

    In this case, the problem will go away when \code{info.tol} is decreased
    slightly.  Only when one has proved that there is no DOR, should one
    use \code{info.tol = 1e-20} which says to ignore the problem altogether
    (whether putative DOR are true or bogus).
}
}
\seealso{
  \code{\link{aster}}, \code{\link{summary}}.
}
\keyword{models}
