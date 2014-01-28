\name{vargplvmLogLikelihood}
\Rdversion{1.0}
\alias{vargplvmLogLikelihood}
\title{Log-likelihood for a variational GP-LVM.}
\description{
  returns the log likelihood for a given GP-LVM model.
}
\usage{
ll <- vargplvmLogLikelihood(model)
}
\arguments{
  \item{model}{the model for which the log likelihood is to be computed. The model contains the data for which the likelihood is being computed in the 'y' component of the structure.}
}
\value{
  \item{ll}{the log likelihood of the data given the model.}
}
\keyword{model}
