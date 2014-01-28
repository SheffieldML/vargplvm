\name{optimiDefaultOptimiser}
\Rdversion{1.0}
\alias{optimiDefaultOptimiser}
\title{Returns the default optimiser to be used.}
\description{
  returns the default optimiser, placing the command here makes it easier to change the default globally for all toolboxes.
}
\usage{
str <- optimiDefaultOptimiser()
}
\arguments{
}
\value{
  \item{str}{string which represents the default optimiser (currently 'scg', the NETLAB scaled conjugate gradient}
}
\seealso{
\code{\link{scg}, \link{conjgrad}, \link{minimize}, \link{optimiDefaultConstraint}}.
}
\keyword{model}
