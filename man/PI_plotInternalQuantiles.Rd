% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Isobole_plotting.R
\name{PI_plotInternalQuantiles}
\alias{PI_plotInternalQuantiles}
\title{Plot spaghetti with quantiles}
\usage{
PI_plotInternalQuantiles(
  grid_fine,
  finalPaths,
  levels = c(0.025, 0.5, 0.975),
  Compound1 = getCompoundInfo("Compound1"),
  Compound2 = getCompoundInfo("Compound2"),
  filename = NULL,
  ...
)
}
\arguments{
\item{grid_fine}{Result from \link{aggregateIsoboles}}

\item{finalPaths}{Result from \link{PI_getFinalPaths}}

\item{levels}{Vector of confidence levels, e.g. c(0.5,0.95). For plotRibbon, named vector with min and max, e.g. c(min = 0.025, max = 0.975)}

\item{Compound1}{Character. Compound name, for cases when it's not available in the Malaria project}

\item{Compound2}{Character. Compound name, for cases when it's not available in the Malaria project}

\item{filename}{Path to save plot}

\item{...}{Arguments goint to \code{\link[ggplot2]{ggsave}}}
}
\value{
ggplot
}
\description{
Plot spaghetti with quantiles
}
\seealso{
Other Internal isobole plotting functions: 
\code{\link{PI_plotInternalAlgoExample}()},
\code{\link{PI_plotInternalMonoDistributions}()},
\code{\link{PI_plotInternalSingleExample}()},
\code{\link{PI_plotInternalSpaghetti}()}
}
\concept{Internal isobole plotting functions}
