\name{gaussianPriorParamInit}
\Rdversion{1.0}
\alias{gaussianPriorParamInit}
\title{Gaussian prior model's parameter initialisation.}
\description{
  initialises the parameters of the Gaussian prior with some default parameters.
}
\usage{
prior <- gaussianPriorParamInit(prior)
}
\arguments{
  \item{prior}{prior structure to be initialised.}
}
\value{
  \item{prior}{prior structure with initial values in place. }
}
\seealso{
\code{\link{priorCreate}, \link{gammaPriorParamInit}}.
}
\keyword{model}
