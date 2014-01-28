\name{modelLogLikeGradients}
\Rdversion{1.0}
\alias{modelLogLikeGradients}
\title{Compute a model's gradients wrt log likelihood.}
\description{
  is a wrapper function to compute the gradients of the log likelihood of a given model.
}
\usage{
g <- modelLogLikeGradients(model)
}
\arguments{
  \item{model}{the model for which likelihoods are computed.}
}
\value{
  \item{g}{teh gradients of the likelihood with respect to the parameters.}
}
\seealso{
\code{\link{modelCreate}}.
}
\keyword{model}
