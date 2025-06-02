\name{astertransform}
\alias{astertransform}
\title{Transform between Aster Model Parameterizations}
\usage{
astertransform(arg, obj, from = c("unconditional", "conditional"),
    to.cond = c("unconditional", "conditional"),
    to.mean = c("mean.value", "canonical"))
}
\description{
Transform between different parameterizations of the aster model.
In effect, this function is called inside \code{\link{predict.aster}}.
Users generally do not need to call it directly.
}
\arguments{
  \item{arg}{canonical parameter vector of length \code{nrow(obj$data)},
    either unconditional (\eqn{\varphi}{phi})
    or conditional (\eqn{\theta}{theta}) depending on the value of
    argument \code{from}.}
  \item{obj}{aster model object, the result of a call to \code{\link{aster}}.}
  \item{from}{the type of canonical parameter which argument \code{arg} is.}
  \item{to.cond}{the type of parameter we want.}
  \item{to.mean}{the type of parameter we want.}
}
\value{
  a vector of the same length as \code{arg}, the transformed parameter vector.
}
\keyword{misc}
