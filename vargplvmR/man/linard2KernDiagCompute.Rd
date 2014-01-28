\name{linard2KernDiagCompute}
\Rdversion{1.0}
\alias{linard2KernDiagCompute}
\title{Compute diagonal of LINARD2 kernel.}
\description{
  computes the diagonal of the kernel matrix for the automatic relevance determination linear kernel given a design matrix of inputs.
}
\usage{
k <- linard2KernDiagCompute(kern, x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\value{
  \item{k}{a vector containing the diagonal of the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{linard2KernParamInit}, \link{kernDiagCompute}, \link{kernCreate}, \link{linard2KernCompute}}.
}
\keyword{model}
