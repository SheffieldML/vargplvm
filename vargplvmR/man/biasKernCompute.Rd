\name{biasKernCompute}
\Rdversion{1.0}
\alias{biasKernCompute}
\title{Compute the BIAS kernel given the parameters and X.}
\description{
  Compute the BIAS kernel given the parameters and X.
}
\usage{
k <- biasKernCompute(kern, x1, x2)
k <- biasKernCompute(kern, x3)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x1}{the input matrix associated with the rows of the kernel.}
  \item{x2}{the inpute matrix associated with the columns of the kernel.}
  \item{x3}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k <- biasKernCompute(kern, x1, x2)}
  computes the kernel parameters for the bias kernel given inputs associated with rows and columns.

  \code{k <- biasKernCompute(kern, x3)}
  computes the kernel matrix for the bias kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{biasKernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{biasKernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
