\name{vargplvmExpandParam}
\Rdversion{1.0}
\alias{vargplvmExpandParam}
\title{Expand a parameter vector into a GP-LVM model.}
\description{
  takes an VARGPLVM structure and a vector of parameters, and fills the structure with the given parameters. Also performs any necessary precomputation for likelihood and gradient computations, so can be computationally intensive to call.
}
\usage{
vargplvmExpandParam(model, params)
}
\arguments{
  \item{model}{the VARGPLVM structure to put the parameters in.}
  \item{params}{parameter vector containing the parameters to put in the VARGPLVM structure. }
}
\value{
}
\seealso{
\code{\link{vargplvmCreate}, \link{vargplvmExtractParam}, \link{modelExpandParam}}.
}
\keyword{model}
