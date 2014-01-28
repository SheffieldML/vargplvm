\name{priorCreate}
\Rdversion{1.0}
\alias{priorCreate}
\title{Create a prior structure given a type.}
\description{
  creates a prior structure given a type.
}
\usage{
prior <- priorCreate(type)
}
\arguments{
  \item{type}{Type of prior to be created,  some standard types are 'gamma', 'gaussian', 'laplace' and 'invgamma'.}
}
\value{
  \item{prior}{The prior structure.}
}
\seealso{
\code{\link{priorParamInit}}.
}
\keyword{model}
