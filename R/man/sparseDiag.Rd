\name{sparseDiag}
\Rdversion{1.0}
\alias{sparseDiag}
\title{Create a diagonal matrix that is sparse from a vector.}
\description{
  creates a diagonal matrix that is sparse from a vector.
}
\usage{
D <- sparseDiag(d)
}
\arguments{
  \item{d}{the diagonal vector from which the sparse diagonal matrix is formed.}
}
\value{
  \item{D}{the sparse diagonal matrix containing the vector as its diagonal.}
}
\seealso{
\code{\link{diag}, \link{spdiags}}.
}
\keyword{model}
