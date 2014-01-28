\name{rbfard2KernDiagCompute}
\Rdversion{1.0}
\alias{rbfard2KernDiagCompute}
\title{Compute diagonal of RBFARD2 kernel.}
\description{
  computes the diagonal of the kernel matrix for the automatic relevance determination radial basis function kernel given a design matrix of inputs.
}
\usage{
k <- rbfard2KernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{rbfard2KernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{rbfard2KernCompute}}.
}
\keyword{model}
