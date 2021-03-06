% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCsimRasch.R
\name{MCsimRasch}
\alias{MCsimRasch}
\title{Generating initial data sample for MCMC simulation.}
\usage{
MCsimRasch(N, splitcr, items)
}
\arguments{
\item{N}{Total sample size for which data matrix should be generated.}

\item{splitcr}{Split criterion which is a numeric vector x with length equal to number of persons and contains zeros and
ones. It indicates group membership for every person. If missing splitcr, by default group 1 with covariate value 0 is set to
\code{ng1 = floor(N/2)}, and group 2 with covariate value 1 set to difference  \code{N - ng1}.}

\item{items}{A vector of item parameters.}
}
\value{
A list of results.
\item{X}{Generated binary data matrix according to the Rasch model.}
\item{splitcr}{see Arguments.}
\item{check}{Logical value. If FALSE data matrix complies to specified conditions. See Details.}
\item{i_total}{Number of iterations needed to generate data matrix.}
\item{call}{The matched call.}
}
\description{
This utility function generates an initial data sample according to the Rasch model.
}
\details{
Restricted and unrestricted models are checked for the conditions of the existence
and uniqueness of the CML estimates according to Fischer (1981) and for appropriate
response pattern, i.e. no only 0/full responses in any item are allowed.
 Vector of person parameters (drawn from a specified distribution) are drawn by default
 at random from the standard normal distribution.
}
\examples{
\dontrun{
# Rasch model
items <- c(0, -2, -1, 0, 1, 2)
x <- c( rep(0, 15), rep(1,15))

res <- MCsimRasch(N = 30, splitcr = x, items = items)



}
}
\references{
{
none
}
}
\seealso{
\code{\link{MCboot}}, \code{\link{MCplot}} and \code{\link{MCplot_crit}}
}
\keyword{mcmc}
