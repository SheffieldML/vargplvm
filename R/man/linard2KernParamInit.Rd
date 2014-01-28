\name{linard2KernParamInit}
\Rdversion{1.0}
\alias{linard2KernParamInit}
\title{LINARD2 kernel parameter initialisation.}
\description{
  initialises the automatic relevance determination linear kernel structure with some default parameters.
}
\usage{
kern <- linard2KernParamInit(kern)
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
