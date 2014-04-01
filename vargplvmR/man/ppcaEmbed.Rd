\name{ppcaEmbed}
\Rdversion{1.0}
\alias{ppcaEmbed}
\title{Embed data set with probabilistic PCA.}
\description{
  returns latent positions for a given data set via probabilistic PCA.
}
\usage{
values <- ppcaEmbed(Y, dims)
}
\arguments{
  \item{Y}{the data set which you want the latent positions for.}
  \item{dims}{the dimensionality of the latent space.}
}
\value{
  \item{X}{the latent positions.}
  \item{sigma2}{the variance not explained by the latent positions.}
  \item{W}{the matrix required to invert the transformation, Y=X*W'.}
}
\seealso{
\code{\link{lleEmbed}, \link{isomapEmbed}}.
}
\keyword{model}
