\name{vargplvmCreate}
\Rdversion{1.0}
\alias{vargplvmCreate}
\title{Create a GPLVM model with inducing variables.}
\description{
  creates a GP-LVM model with the possibility of using inducing variables to speed up computation.
}
\usage{
model <- vargplvmCreate(q, d, Y, options, enableDgtN)
}
\arguments{
  \item{q}{dimensionality of latent space.}
  \item{d}{dimensionality of data space.}
  \item{Y}{the data to be modelled in design matrix format (as many rows as there are data points).}
  \item{options}{options structure as returned from VARGPLVMOPTIONS. This structure determines the type of approximations to be used (if any).}
  \item{enableDgtN}{if set to true, the model will be created in D >> N mode (if this is indeed the case).}
}
\value{
  \item{model}{the GP-LVM model.}
}
\seealso{
\code{\link{vargplvmOptions}}.
}
\keyword{model}
