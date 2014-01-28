\name{linard2KernDiagGradX}
\Rdversion{1.0}
\alias{linard2KernDiagGradX}
\title{Gradient of LINARD2 kernel's diagonal with respect to X.}
\description{
  computes the gradient of the diagonal of the automatic relevance determination linear kernel matrix with respect to the elements of the design matrix given in X.
}
\usage{
gX <- linard2KernDiagGradX(kern, X)
}
\arguments{
  \item{kern}{the kernel structure for which gradients are being computed.}
  \item{X}{the input data in the form of a design matrix.}
}
\value{
  \item{gX}{the gradients of the diagonal with respect to each element of X. The returned matrix has the same dimensions as X.}
}
\seealso{
\code{\link{linard2KernParamInit}, \link{kernDiagGradX}, \link{linard2kernGradX}}.
}
\keyword{model}
