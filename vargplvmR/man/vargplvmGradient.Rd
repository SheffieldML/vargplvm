\name{vargplvmGradient}
\Rdversion{1.0}
\alias{vargplvmGradient}
\title{Variational GP-LVM gradient wrapper.}
\description{
  is a wrapper function for the gradient of the negative log likelihood of an variatioanl GP-LVM model with respect to the latent postions and parameters.
}
\usage{
g <- vargplvmGradient(params, model)
}
\arguments{
  \item{params}{vector of parameters and latent postions where the gradient is to be evaluated.}
  \item{model}{the model structure into which the latent positions and the parameters will be placed.}
}
\value{
  \item{g}{the gradient of the negative log likelihood with respect to the latent positions and the parameters at the given point. }
}
\seealso{
\code{\link{vargplvmLogLikeGradients}, \link{vargplvmExpandParam}}.
}
\keyword{model}
