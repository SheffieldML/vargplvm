\name{vargplvmOptimise}
\Rdversion{1.0}
\alias{vargplvmOptimise}
\title{Optimise the VARGPLVM.}
\description{
  takes a given GP-LVM model structure and optimises with respect to parameters and latent positions.
}
\usage{
model <- vargplvmOptimise(model, display, iters)
}
\arguments{
  \item{model}{the model to be optimised.}
  \item{display}{flag dictating whether or not to display optimisation progress (set to greater than zero) (default value 1).}
  \item{iters}{number of iterations to run the optimiser for (default value 2000).}
}
\value{
  \item{model}{the optimised model.}
}
\seealso{
\code{\link{vargplvmCreate}, \link{vargplvmLogLikelihood}, \link{vargplvmLogLikeGradients}, \link{vargplvmObjective}, \link{vargplvmGradient}}.
}
\keyword{model}
