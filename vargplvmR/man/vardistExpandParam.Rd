\name{vardistExpandParam}
\Rdversion{1.0}
\alias{vardistExpandParam}
\title{Expand a parameter vector into a vardist structure.}
\description{
  takes an VARDIST structure and a vector of parameters, and fills the structure with the given parameters. Also performs any necessary precomputation for likelihood and gradient computations, so can be computationally intensive to call.
}
\usage{
vardistExpandParam(model, params)
}
\arguments{
  \item{model}{the VARDIST structure to put the parameters in.}
  \item{params}{parameter vector containing the parameters to put in the VARDIST structure. }
}
\value{
}
\seealso{
\code{\link{vardistCreate}, \link{vardistExtractParam}, \link{modelExpandParam}}.
}
\keyword{model}
