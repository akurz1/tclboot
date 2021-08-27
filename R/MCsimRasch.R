######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# MCsimRasch
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# This file contains a routine that generates initial data sample
#
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 10/05/2021
######################################################################
#' Generating initial data sample for MCMC simulation.
#'
#' This utility function generates an initial data sample according to the Rasch model.
#'
#' Restricted and unrestricted models are checked for the conditions of the existence
#' and uniqueness of the CML estimates according to Fischer (1981) and for appropriate
#' response pattern, i.e. no only 0/full responses in any item are allowed.
#'  Vector of person parameters (drawn from a specified distribution) are drawn by default
#'  at random from the standard normal distribution.
#'
#' @param N Total sample size for which data matrix should be generated.
#' @param splitcr Split criterion which is a numeric vector x with length equal to number of persons and contains zeros and
#' ones. It indicates group membership for every person. If missing splitcr, by default group 1 with covariate value 0 is set to
#' \code{ng1 = floor(N/2)}, and group 2 with covariate value 1 set to difference  \code{N - ng1}.
#' @param items A vector of item parameters.
#'
#' @return A list of results.
#' \item{X}{Generated binary data matrix according to the Rasch model.}
#' \item{splitcr}{see Arguments.}
#' \item{check}{Logical value. If FALSE data matrix complies to specified conditions. See Details.}
#' \item{i_total}{Number of iterations needed to generate data matrix.}
#' \item{call}{The matched call.}
#'
#' @references{
#' none
#' }
#' @keywords mcmc
#' @export
#' @seealso \code{\link{MCboot}}, \code{\link{MCplot}} and \code{\link{MCplot_crit}}
#' @examples
#' \dontrun{
#' # Rasch model
#' items <- c(0, -2, -1, 0, 1, 2)
#' x <- c( rep(0, 15), rep(1,15))
#'
#' res <- MCsimRasch(N = 30, splitcr = x, items = items)
#'
#'
#'
#'}


MCsimRasch  <- function(N,
                        splitcr,
                        items){
  call<-match.call()

  if (missing(splitcr)) {
    ng1 <- floor(N/2)
    splitcr <- c(rep(0, ng1), rep(1, N - ng1)) # binary covariate
  }

  mci <- mc_inputmat(splitcr = splitcr, beta = items)
  mci$call <- call

  return(mci)
}
