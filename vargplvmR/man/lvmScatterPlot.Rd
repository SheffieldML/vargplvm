\name{lvmScatterPlot}
\Rdversion{1.0}
\alias{lvmScatterPlot}
\title{2-D scatter plot of the latent points.}
\description{
  produces a visualisation of the latent space for the given model, using the provided labels to distinguish the latent points.
}
\usage{
values <- lvmScatterPlot(model, model, lbls, model, lbls, ax)
}
\arguments{
  \item{model}{the model for which the scatter plot is being produced.}
  \item{model}{the model for which the scatter plot is being produced.}
  \item{lbls}{labels for each data point so that they may be given different symbols. Useful when each data point is associated with a different class.}
  \item{model}{the model for which the scatter plot is being produced.}
  \item{lbls}{labels for each data point so that they may be given different symbols. Useful when each data point is associated with a different class.}
  \item{ax}{the axes where the plot is to be placed.}
}
\value{
  \item{ax}{the axes handle where the scatter plot was placed.}
  \item{ax}{the axes handle where the scatter plot was placed. }
  \item{ax}{the axes handle where the scatter plot was placed.}
}
\seealso{
\code{\link{lvmVisualise}, \link{lvmTwoDPlot}, \link{lvmScatterPlotColor}}.
}
\keyword{model}
