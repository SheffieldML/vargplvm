\name{linard2KernCompute}
\Rdversion{1.0}
\alias{linard2KernCompute}
\title{Compute the LINARD2 kernel given the parameters and X.}
\description{
  Compute the LINARD2 kernel given the parameters and X.
}
\usage{
k <- linard2KernCompute(kern, x, x2)
k <- linard2KernCompute(kern, x3)
}
\arguments{
  \item{kern}{the kernel structure for which the matrix is computed.}
  \item{x}{the input matrix associated with the rows of the kernel.}
  \item{x2}{the input matrix associated with the columns of the kernel.}
  \item{x3}{input data matrix in the form of a design matrix.}
}
\details{
  \code{k <- linard2KernCompute(kern, x, x2)}
  computes the kernel parameters for the automatic relevance determination linear kernel given inputs associated with rows and columns.

  \code{k <- linard2KernCompute(kern, x3)}
  computes the kernel matrix for the automatic relevance determination linear kernel given a design matrix of inputs.

}
\value{
  \item{k}{the kernel matrix computed at the given points.}
}
\seealso{
\code{\link{linard2KernParamInit}, \link{kernCompute}, \link{kernCreate}, \link{linard2KernDiagCompute}}.
}
\examples{
## missing
}
\keyword{model}
