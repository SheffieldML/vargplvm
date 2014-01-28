\name{vargplvmReduceModel}
\Rdversion{1.0}
\alias{vargplvmReduceModel}
\title{prunes out dimensions of the model.}
\description{
  order the latent dimensions acorrding to the inputScales and reduces the model to have smaller number of latent dimensions.
}
\usage{
model <- vargplvmReduceModel(model, P, dims)
}
\arguments{
  \item{model}{the model to be reduced.}
  \item{P}{the number of dimensions to move to (setting to model.q will just reorder the dimensions in the model).}
  \item{dims}{(optional) explicit set of dimensions to use}
}
\value{
  \item{model}{the model with the reduced number of dimensions.}
}
\seealso{
\code{\link{vargplvmCreate}}.
}
\keyword{model}
