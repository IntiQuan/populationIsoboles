% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Isobole_plotting.R
\name{aggregateIsoboles}
\alias{aggregateIsoboles}
\title{Calculate Confidence Level Response Surface}
\usage{
aggregateIsoboles(
  .inputFolder,
  gridInfo = NULL,
  finalPaths = NULL,
  FLAGoverwriteQuantiles = determineRecalcNeccessary(.inputFolder = .inputFolder,
    fileToRecalc = "GLOBAL-101-percentages-grid_fine.rds")
)
}
\arguments{
\item{.inputFolder}{path/to/Simulations}

\item{gridInfo}{Output of \link{PI_readGridInfo}}

\item{finalPaths}{Result from \link{PI_getFinalPaths}}

\item{FLAGoverwriteQuantiles}{Redo the calculation of the quantiles?}
}
\value{
data.table(x,y, percentages)
}
\description{
For each point on the full grid and each population,
check if it lies left or right of the population isobole and
assign 0 or 1 to it.
Then, calculate the mean for each point.
This is equivalent to the confidence level response surface
}
\details{
This function will look in .inputFolder for a file called "GLOBAL-101-percentages-grid_fine.rds".
If it is found and FLAGoverwriteQuantiles = FALSE, it will load this result
}
\author{
Daniel Lill (daniel.lill@intiquan.com)
}
\concept{Isobole simulation post processing}
