\name{biasKernExtractParam}
\Rdversion{1.0}
\alias{biasKernExtractParam}
\title{Extract parameters from the BIAS kernel structure.}
\description{
  Extract parameters from the BIAS kernel structure.
}
\usage{
param <- biasKernExtractParam(kern)
values <- biasKernExtractParam(kern)
}
\arguments{
  \item{kern}{the kernel structure containing the parameters to be extracted.}
}
\details{
  \code{param <- biasKernExtractParam(kern)}
  extracts parameters from the bias kernel structure into a vector of parameters for optimisation.

  \code{param, names <- biasKernExtractParam(kern)}
  extracts parameters and parameter names from the bias kernel structure.

}
\value{
  \item{param}{vector of parameters extracted from the kernel. If the field 'transforms' is not empty in the kernel matrix, the parameters will be transformed before optimisation (for example positive only parameters could be logged before being returned).}
  \item{names}{cell array of strings containing parameter names.}
}
\seealso{
\code{\link{biasKernParamInit}, \link{biasKernExpandParam}, \link{kernExtractParam}, \link{scg}, \link{conjgrad}}.
}
\examples{
## missing
}
\keyword{model}
