% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Isobole_gridRefinementAlgorithm.R
\name{fastIsoboles_iteration0}
\alias{fastIsoboles_iteration0}
\title{Initialize the grid refinement algorithm}
\usage{
fastIsoboles_iteration0(
  gridmin,
  gridmax,
  gridlen_i = (gridmax - gridmin)/2,
  i = 0
)
}
\arguments{
\item{gridmin, gridmax}{as in fastIsoboles_iteration}

\item{gridlen_i, i}{as in fastIsoboles_iteration. gridlen_i and i are only there for backward compatibility}
}
\value{
list with grid and path as explained in the Arguments-section of fastIsoboles_iteration
}
\description{
In iteration 0, a dummy-grid and a dummy-path object are returned which will be used in iteration 1 to calculate objvalue values.
The algorithm evaluates only objvalue points which are close to the path of the isobole and haven't yet been evaluated.
The dummy path contains all grid points ensuring that all points of the grid in iteration 1 lie close to the path.
The dummy grid contains contains all grid points but evaluatedCum is set to zero, therefore the points will be evaluated next time.
}
\examples{
populationIsoboles:::fastIsoboles_iteration0(c(x = 0,y = 0), c(x = 16,y = 16), c(x = 8,y = 8))
}
\seealso{
Other Grid refinement functions for isoboles: 
\code{\link{arrange_isobole_internal}()},
\code{\link{fastIsoboles_iteration}()},
\code{\link{fastIsoboles}()},
\code{\link{get_gridlens}()},
\code{\link{get_isobole_path}()}
}
\author{
Daniel Lill, IntiQuan \email{daniel.lill@intiquan.com}
}
\concept{Grid refinement functions for isoboles}
