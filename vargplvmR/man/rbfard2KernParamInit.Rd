\name{rbfard2KernParamInit}
\Rdversion{1.0}
\alias{rbfard2KernParamInit}
\title{RBFARD2 kernel parameter initialisation.}
\description{
  initialises the automatic relevance determination radial basis function kernel structure with some default parameters.
}
\usage{
kern <- rbfard2KernParamInit(kern)
}
\arguments{
  \item{kern}{the kernel structure which requires initialisation.}
}
\value{
  \item{kern}{the kernel structure with the default parameters placed in.}
}
\seealso{
\code{\link{kernCreate}, \link{kernParamInit}}.
}
\keyword{model}
