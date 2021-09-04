######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# MCplot
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# This file contains a routine that plots the power function
#
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 12/05/2021
######################################################################
#' Plot of power function.
#'
#'This function plots the power function of selected test statistics
#'
#' @param object An object of class "\code{MCplot}"
#' @param nstats Name(s) of test statistics to be plotted. Note a maximum of 4 tests allowed. By default c("W", "LR", "RS", "G") will be plotted.
#' Note "U1" = "sqs" or sum of squared elements of the score function and "U2" = "abs" or sum of absolute values.
#' @param alpha Probability of error of first kind (in plot shown as dotted horizontal line).
#' @param xlim  A numeric value, specifying the left/lower limit and the right/upper limit of the x scale. Maximum value is xlim = 10.
#' @param legend_title The text for Legend. By default it is set to "Test".
#' @param tag_title The text for tag label, which will be displayed at the top-left of the plot, by default of Item 1 only. By default it is element_blank().
#' @param name_title The text for title of plot. By default it is element_blank().
#' @param tcolor Logic value. If TRUE then ggplot plot will be in color otherwise black and white.
#'
#' @return An object of class "\code{MCplot}" containing:
#'  \item{plotlist}{A list of length k containing ggplot objects of power plots of each item and selected test statistics.}
#'  \item{call}{The matched call.}
#'
#' @references{
#' Draxler, C., & Dahm, S. (2020). Conditional or Pseudo Exact Tests with an Application in the Context of Modeling
#' Response Times. Psych, 2(4), 198-208. \url{https://doi.org/10.3390/psych2040017}
#' }
#' @keywords mcmc
#' @export
#' @examples
#' \dontrun{
#' # Rasch model with sample size N = 60, k = 8
#' library(tclboot)
#'
#' # model parameters
#' items <- c(0, -2, -1, -0.5, 0, 0.5, 1, 2)
#' x <- c(rep(0,25), rep(1,35))
#'
#' # generate initial data matrix
#' y <- MCsimRasch(N = 60, splitcr = x, items = items)$X
#'
#' obj <- MCboot(X=y, splitcr = x, n_iter = 1000)
#'
#' res <- MCplot(object = obj, nstats = c("W", "G"))
#'
#' }

MCplot <- function(object,
                   nstats = c("W", "LR", "RS", "G", "abs", "sqs", "St"),
                   alpha = 0.05,
                   xlim = 5,
                   legend_title = ggplot2::element_blank(), #"Test",
                   legend.direction = "vertical",
                   tag_title = ggplot2::element_blank(),
                   name_title = ggplot2::element_blank(),
                   tcolor=FALSE){
  call<-match.call()
  if (length(nstats) > 4)  nstats <- c("W", "LR", "RS", "G")

  # plotlist <- plotlist(object, nstats, alpha, xlim, legend_title, tag_title, tcolor)
  plotlist <- plotlist(object = object,
                       nstats = nstats,
                       alpha = alpha,
                       xlim = xlim,
                       legend_title = legend_title,
                       legend.direction = legend.direction,
                       tag_title = tag_title,
                       name_title = name_title,
                       tcolor = tcolor)

  results <- list(plotlist = plotlist,
                  "call" = call)

  class(results) <- "MCplot"

  return(results)

}
