\name{vargplvmPosteriorVar}
\Rdversion{1.0}
\alias{vargplvmPosteriorVar}
\title{variances of the posterior at points given by X.}
\description{
  returns the posterior mean and variance for a given set of points.
}
\usage{
values <- vargplvmPosteriorVar(model, x)
}
\arguments{
  \item{model}{the model for which the posterior will be computed.}
  \item{x}{the input positions for which the posterior will be computed.}
}
\value{
  \item{mu}{the mean of the posterior distribution.}
  \item{sigma}{the variances of the posterior distributions.}
}
\seealso{
\code{\link{gpPosteriorMeanVar}, \link{vargplvmCreate}}.
}
\keyword{model}
