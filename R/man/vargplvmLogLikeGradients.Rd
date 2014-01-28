\name{vargplvmLogLikeGradients}
\Rdversion{1.0}
\alias{vargplvmLogLikeGradients}
\title{Compute the gradients for the variational GPLVM.}
\description{
  Compute the gradients for the variational GPLVM.
}
\usage{
g <- vargplvmLogLikeGradients(model)
values <- vargplvmLogLikeGradients(model)
}
\arguments{
  \item{model}{the FGPLVM structure containing the parameters and the latent positions.}
}
\details{
  \code{g <- vargplvmLogLikeGradients(model)}
  returns the gradients of the log likelihood with respect to the parameters of the GP-LVM model and with respect to the latent positions of the GP-LVM model.

  \code{gX, gParam <- vargplvmLogLikeGradients(model)}
  returns the gradients of the log likelihood with respect to the parameters of the GP-LVM model and with respect to the latent positions of the GP-LVM model in seperate matrices.

}
\value{
  \item{g}{the gradients of the latent positions (or the back constraint's parameters) and the parameters of the GP-LVM model.}
  \item{gX}{the gradients of the latent positions (or the back constraint's parameters).}
  \item{gParam}{gradients of the parameters of the GP-LVM model.}
}
\seealso{
\code{\link{vargplvmLogLikelihood}, \link{vargplvmCreate}, \link{modelLogLikeGradients}}.
}
\examples{
## missing
}
\keyword{model}
