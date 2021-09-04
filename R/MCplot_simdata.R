######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# MCplot_simdata
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# This file contains a routine that plots the power function of
# simulated data.
#
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 12/05/2021
######################################################################
#' Plot of power function.
#'
#'This function plots the power function of selected test statistics
#'
#' @param N Sample size \eqn{n} of selected case.  N =  10, 15, 20, 25, 30, 50
#' @param k Number of items \eqn{k} of selected case. k = 6, 8, 10, 12, 14, 16, 18, 20, 25, 30
#' @param nstats Name(s) of test statistics to be plotted. Note a maximum of 4 tests allowed. By default c("W", "LR", "RS", "G") will be plotted.
#' Note "U1" = "sqs" or sum of squared elements of the score function and "U2" = "abs" or sum of absolute values.
#' @param alpha Probability of error of first kind (in plot shown as dotted horizontal line).
#' @param xlim  A numeric value, specifying the left/lower limit and the right/upper limit of the x scale. Maximum value is xlim = 10.
#' @param legend_title The text for Legend. By default it is set to "Test".
#' @param tag_title The text for tag label, which will be displayed at the top-left of the plot, by default of Item 1 only. By default it is element_blank().
#' @param name_title The text for title of plot. By default it is element_blank(),
#' @param legend.position the position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param legend.direction layout of items in legends ("horizontal" or "vertical")
#' @param tcolor Logic value. If TRUE then ggplot plot will be in color otherwise black and white.
#'
#' @return An object of class "\code{MCplot.case}" containing:
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
#' # load MC sim data of case G-V-30x18
#' # with N=30, k=18
#' # and generate power plots of selected tests
#'
#'res <- MCplot_simdata(N=30, k=18, nstats = c("LR", "G"))
#'
#' }

MCplot_simdata <- function(N, k,
                           nstats = c("W", "LR", "RS", "G", "abs", "sqs", "St"),
                           alpha = 0.05,
                           xlim = 5,
                           legend_title = ggplot2::element_blank(),
                           tag_title = ggplot2::element_blank(),
                           name_title = ggplot2::element_blank(),
                           legend.position = c(.5, .75), # "top",
                           legend.direction = "vertical",
                           tcolor=FALSE){
  call<-match.call()
  if (length(nstats) > 4)  nstats <- c("W", "LR", "RS", "G")
  # legend_title = "Test" # expression(italic("n"))

  obj <- load_data(N = N, k = k) # load simulation data und results
  class(obj) <- "MCplot_case"

  plotlist <- plotlist(object=obj,
                       nstats=nstats,
                       alpha=alpha,
                       xlim=xlim,
                       legend_title=legend_title,
                       tag_title=tag_title,
                       name_title=name_title,
                       legend.position = legend.position,
                       legend.direction = legend.direction,
                       tcolor=tcolor)

  results <- list(plotlist = plotlist,
                  object = obj,
                  "call" = call)

  class(results) <- "MCplot_case"

  return(results)

}
