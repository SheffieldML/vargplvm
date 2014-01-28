\name{rbfard2KernExtractParam}
\Rdversion{1.0}
\alias{rbfard2KernExtractParam}
\title{Extract parameters from the RBFARD2 kernel structure.}
\description{
  Extract parameters and parameter names from the automatic relevance determination radial basis function kernel structure.
}
\usage{
values <- rbfard2KernExtractParam(kern, kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\value{
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel matrix, the parameters will be transformed before optimisation (for example positive only parameters could be logged before being returned).}
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel matrix, the parameters will be transformed before optimisation (for example positive only parameters could be logged before being returned).}
  \item{names}{cell array of strings containg the parameter names.}
}
\seealso{
\code{\link{}}.
}
\keyword{model}
