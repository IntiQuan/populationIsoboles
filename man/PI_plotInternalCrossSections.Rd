% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Isobole_plotting.R
\name{PI_plotInternalCrossSections}
\alias{PI_plotInternalCrossSections}
\title{See PI_plotInternalMonoDistributions}
\usage{
PI_plotInternalCrossSections(
  grid_fine,
  cross_x = NULL,
  cross_y = NULL,
  filename = NULL,
  ...
)
}
\arguments{
\item{grid_fine}{Output from aggregateIsoboles}

\item{cross_x, cross_y}{Doses at which the cross section should be evaluated}

\item{filename, ...}{going to \code{\link[ggplot2:ggsave]{ggplot2::ggsave()}}}
}
\value{
ggplot
}
\description{
See PI_plotInternalMonoDistributions
}
\author{
Daniel Lill (daniel.lill@intiquan.com)
}
