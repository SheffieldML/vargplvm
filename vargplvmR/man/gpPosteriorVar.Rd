\name{gpPosteriorVar}
\Rdversion{1.0}
\alias{gpPosteriorVar}
\title{Variances of the posterior at points given by X.}
\description{
  returns the posterior variance for a given set of points.
}
\usage{
sigma <- gpPosteriorVar(model, x)
}
\arguments{
  \item{model}{the model for which the posterior will be computed.}
  \item{x}{the input positions for which the posterior will be computed.}
}
\value{
  \item{sigma}{the variances of the posterior distributions.}
}
\seealso{
\code{\link{gpCreate}, \link{gpPosteriorMeanCovar}, \link{gpPosteriorGradMeanVar}, \link{gpPosteriorMeanVar}}.
}
\keyword{model}
