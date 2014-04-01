\name{biasKernDiagCompute}
\Rdversion{1.0}
\alias{biasKernDiagCompute}
\title{Compute diagonal of BIAS kernel.}
\description{
  computes the diagonal of the kernel matrix for the bias kernel given a design matrix of inputs.
}
\usage{
k <- biasKernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{biasKernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{biasKernCompute}}.
}
\keyword{model}
