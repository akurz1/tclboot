% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCplot_simdata.R
\name{MCplot_simdata}
\alias{MCplot_simdata}
\title{Plot of power function.}
\usage{
MCplot_simdata(
  N,
  k,
  nstats = c("W", "LR", "RS", "G", "abs", "sqs", "St"),
  alpha = 0.05,
  xlim = 5,
  legend_title = ggplot2::element_blank(),
  tag_title = ggplot2::element_blank(),
  name_title = ggplot2::element_blank(),
  legend.position = c(0.5, 0.75),
  legend.direction = "vertical",
  tcolor = FALSE
)
}
\arguments{
\item{N}{Sample size \eqn{n} of selected case.  N =  10, 15, 20, 25, 30, 50}

\item{k}{Number of items \eqn{k} of selected case. k = 6, 8, 10, 12, 14, 16, 18, 20, 25, 30}

\item{nstats}{Name(s) of test statistics to be plotted. Note a maximum of 4 tests allowed. By default c("W", "LR", "RS", "G") will be plotted.
Note "U1" = "sqs" or sum of squared elements of the score function and "U2" = "abs" or sum of absolute values.}

\item{alpha}{Probability of error of first kind (in plot shown as dotted horizontal line).}

\item{xlim}{A numeric value, specifying the left/lower limit and the right/upper limit of the x scale. Maximum value is xlim = 10.}

\item{legend_title}{The text for Legend. By default it is set to "Test".}

\item{tag_title}{The text for tag label, which will be displayed at the top-left of the plot, by default of Item 1 only. By default it is element_blank().}

\item{name_title}{The text for title of plot. By default it is element_blank(),}

\item{legend.position}{the position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)}

\item{legend.direction}{layout of items in legends ("horizontal" or "vertical")}

\item{tcolor}{Logic value. If TRUE then ggplot plot will be in color otherwise black and white.}
}
\value{
An object of class "\code{MCplot.case}" containing:
 \item{plotlist}{A list of length k containing ggplot objects of power plots of each item and selected test statistics.}
 \item{call}{The matched call.}
}
\description{
This function plots the power function of selected test statistics
}
\examples{
\dontrun{
# load MC sim data of case G-V-30x18
# with N=30, k=18
# and generate power plots of selected tests

res <- MCplot_simdata(N=30, k=18, nstats = c("LR", "G"))

}
}
\references{
{
Draxler, C., & Dahm, S. (2020). Conditional or Pseudo Exact Tests with an Application in the Context of Modeling
Response Times. Psych, 2(4), 198-208. \url{https://doi.org/10.3390/psych2040017}
}
}
\keyword{mcmc}
