\name{jitChol}
\Rdversion{1.0}
\alias{jitChol}
\title{Do a Cholesky decomposition with jitter.}
\description{
  Do a Cholesky decomposition with jitter.
}
\usage{
U <- jitChol(A, maxTries)
values <- jitChol(A, maxTries)
}
\arguments{
  \item{A}{the matrix for which the Cholesky decomposition is required.}
  \item{maxTries}{the maximum number of times that jitter is added before giving up (default 10).}
}
\details{
  \code{U <- jitChol(A, maxTries)}
  attempts a Cholesky decomposition on the given matrix, if matrix isn't positive definite the function gives a warning, adds 'jitter' and tries again. At the first attempt the amount of jitter added is 1e-6 times the mean of the diagonal. Thereafter the amount of jitter is multiplied by 10 each time it is added again. This is continued for a maximum of 10 times.

  \code{U, jitter <- jitChol(A, maxTries)}
  attempts a Cholesky decomposition on the given matrix, if matrix isn't positive definite the function adds 'jitter' and tries again. Thereafter the amount of jitter is multiplied by 10 each time it is added again. This is continued for a maximum of 10 times.  The amount of jitter added is returned.

}
\value{
  \item{U}{the Cholesky decomposition for the matrix.}
  \item{jitter}{the amount of jitter that was added to the matrix.}
}
\seealso{
\code{\link{chol}, \link{pdinv}, \link{logdet}}.
}
\examples{
## missing
}
\keyword{model}
