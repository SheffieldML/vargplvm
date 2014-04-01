\name{linard2KernExpandParam}
\Rdversion{1.0}
\alias{linard2KernExpandParam}
\title{Create kernel structure from LINARD2 kernel's parameters.}
\description{
  returns a automatic relevance determination linear kernel structure filled with the parameters in the given vector. This is used as a helper function to enable parameters to be optimised in, for example, the NETLAB optimisation functions.
}
\usage{
kern <- linard2KernExpandParam(kern, param)
}
\arguments{
  \item{kern}{the kernel structure in which the parameters are to be placed.}
  \item{param}{vector of parameters which are to be placed in the kernel structure.}
}
\value{
  \item{kern}{kernel structure with the given parameters in the relevant locations.}
}
\seealso{
\code{\link{linard2KernParamInit}, \link{linard2KernExtractParam}, \link{kernExpandParam}}.
}
\keyword{model}
