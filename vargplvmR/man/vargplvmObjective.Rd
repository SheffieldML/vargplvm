\name{vargplvmObjective}
\Rdversion{1.0}
\alias{vargplvmObjective}
\title{Wrapper function for variational GP-LVM objective.}
\description{
  provides a wrapper function for the variational GP-LVM, it takes the negative of the log likelihood, feeding the parameters correctly to the model.
}
\usage{
f <- vargplvmObjective(params, model)
}
\arguments{
  \item{params}{the parameters of the variational GP-LVM model.}
  \item{model}{the model structure in which the parameters are to be placed.}
}
\value{
  \item{f}{the negative of the log likelihood of the model. }
}
\seealso{
\code{\link{vargplvmCreate}, \link{vargplvmLogLikelihood}, \link{vargplvmExpandParam}}.
}
\keyword{model}
