\name{vargplvmOptions}
\Rdversion{1.0}
\alias{vargplvmOptions}
\title{Return default options for VARGPLVM model.}
\description{
  returns the default options in a structure for a VARGPLVM model.
}
\usage{
options <- vargplvmOptions(approx)
}
\arguments{
  \item{approx}{approximation type, either 'ftc' (no approximation), 'dtc' (deterministic training conditional), 'dtcvar', variational sparse approximation, 'fitc' (fully independent training conditional) or 'pitc' (partially independent training conditional.}
}
\value{
  \item{options}{structure containing the default options for the given approximation type.}
}
\seealso{
\code{\link{vargplvmCreate}}.
}
\keyword{model}
