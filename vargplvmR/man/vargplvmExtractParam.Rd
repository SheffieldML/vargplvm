\name{vargplvmExtractParam}
\Rdversion{1.0}
\alias{vargplvmExtractParam}
\title{Extract a parameter vector from a variational GP-LVM model.}
\description{
  does the same as above, but also returns parameter names.
}
\usage{
values <- vargplvmExtractParam(model, model2)
}
\arguments{
  \item{model}{the model from which parameters are to be extracted.}
  \item{model2}{the model structure containing the information about the model.}
}
\value{
  \item{params}{the parameter vector extracted from the model.}
  \item{params}{a vector of parameters from the model.}
  \item{names}{cell array of parameter names.}
}
\seealso{
\code{\link{vargplvmCreate}, \link{vargplvmExpandParam}, \link{modelExtractParam}}.
}
\keyword{model}
