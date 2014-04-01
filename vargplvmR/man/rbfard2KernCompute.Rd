\name{rbfard2KernCompute}
\Rdversion{1.0}
\alias{rbfard2KernCompute}
\title{Compute the RBFARD kernel given the parameters and X.}
\description{
  Compute the RBFARD kernel given the parameters and X.
}
\usage{
k <- rbfard2KernCompute(kern, x1, x2)
k <- rbfard2KernCompute(x)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x1}{the input matrix associated with the rows of the kernel.}
  \item{x2}{the input matrix associated with the columns of the kernel.}
  \item{x}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k <- rbfard2KernCompute(kern, x1, x2)}
  computes the kernel parameters for the automatic relevance determination radial basis function kernel given inputs associated with rows and columns.

  \code{k <- rbfard2KernCompute(x)}
  computes the kernel matrix for the automatic relevance determination radial basis function kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{rbfard2KernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{rbfard2KernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
