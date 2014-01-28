\name{pdinv}
\Rdversion{1.0}
\alias{pdinv}
\title{Invert a positive definite matrix.}
\description{
  Invert a positive definite matrix.
}
\usage{
values <- pdinv(A)
values <- pdinv(A, U)
values <- pdinv(A, U)
}
\arguments{
  \item{A}{the input positive definite matrix to be inverted.}
  \item{U}{the Cholesky decomposition of A.}
}
\details{
  \code{Ainv, U <- pdinv(A)}
  inverts a positive definite matrix. If the matrix isn't quite positive definite the function adds 'jitter' to make it positive definite and gives out a warning message (this is done through JITCHOL).

  \code{Ainv, U <- pdinv(A, U)}
  inverts a positive definite matrix given the Cholesky decomposition of A.

  \code{Ainv, U, jitter <- pdinv(A, U)}
  inverts a positive definite matrix given the Cholesky decomposition of A. If jitter is used then the amount of jitter used is returned.

}
\value{
  \item{Ainv}{the inverse of A computed using Cholesky decomposition.}
  \item{U}{the Cholesky decomposition of A.}
  \item{jitter}{the amount of jitter added.}
}
\seealso{
\code{\link{jitChol}, \link{logdet}, \link{chol}}.
}
\examples{
## missing
}
\keyword{model}
