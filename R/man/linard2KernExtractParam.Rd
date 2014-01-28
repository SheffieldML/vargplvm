\name{linard2KernExtractParam}
\Rdversion{1.0}
\alias{linard2KernExtractParam}
\title{Extract parameters from the LINARD2 kernel structure.}
\description{
  Extract parameters from the automatic relevance determination linear kernel structure into a vector of parameters for optimisation.
}
\usage{
param <- linard2KernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\value{
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel matrix, the parameters will be transformed before optimisation (for example}
}
\keyword{model}
