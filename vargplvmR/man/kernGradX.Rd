\name{kernGradX}
\Rdversion{1.0}
\alias{kernGradX}
\title{Compute the gradient of the kernel wrt X.}
\description{
  Compute the gradient of the kernel wrt X.
}
\usage{
g <- kernGradX(kern, x)
g <- kernGradX(kern, x1, x2)
}
\arguments{
  \item{kern}{kernel structure for which gradients are being computed.}
  \item{x}{locations against which gradients are being computed.}
  \item{x1}{row locations against which gradients are being computed.}
  \item{x2}{column locations against which gradients are being computed.}
}
\details{
  \code{g <- kernGradX(kern, x)}
  computes the gradient of the kernel with respect to the input positions.

  \code{g <- kernGradX(kern, x1, x2)}
  computes the gradient of the kernel with respect to the input positions where both the row positions and column positions are provided separately.

}
\value{
  \item{g}{the returned gradients. The gradients are returned in a matrix which is numData x numInputs x numData. Where numData is the number of data points and numInputs is the number of input dimensions in X.}
}
\seealso{
\code{\link{kernDiagGradX}, \link{kernGradient}}.
}
\examples{
## missing
}
\keyword{model}
