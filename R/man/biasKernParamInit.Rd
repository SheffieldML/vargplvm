\name{biasKernParamInit}
\Rdversion{1.0}
\alias{biasKernParamInit}
\title{BIAS kernel parameter initialisation.}
\description{
  initialises the bias kernel structure with some default parameters.
}
\usage{
kern <- biasKernParamInit(kern)
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
