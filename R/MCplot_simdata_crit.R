######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# MCplot_simdata_crit
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
#' @param nstats Name of test statistics to be plotted. Note only one test allowed. By default "sqs" will be plotted.
#' @param alpha Probability of error of first kind (in plot shown as dotted horizontal line).
#' @param xlim  A numeric value, specifying the left/lower limit and the right/upper limit of the x scale. Maximum value is xlim = 10.
#' @param legend_title The text for Legend. By default it is set to "Test".
#' @param tag_title The text for tag label, which will be displayed at the top-left of the plot, by default of Item 1 only. By default it is element_blank().
#' @param name_title The text for title of plot. By default it is element_blank().
#' @param legend.position the position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param legend.direction layout of items in legends ("horizontal" or "vertical")
#' @param tcolor Logic value. If TRUE then ggplot plot will be in color otherwise black and white.
#' @param alpha_mod A numeric value. The crit. region \eqn{C} of selected test statistic will be \eqn{100(alpha*alpha_mod)} percent largest values of this test statistics. Default value is 1.
#' @param ctr_mod Character. Which type of modification of critical region should be used. Adding low most extreme values of the respective of \eqn{T} statistics ("low"),
#' adding high most extreme values of the respective of \eqn{T} statistics ("high"), or adding both ("both"). The latter is also Default value.
#'
#' @return An object of class "\code{MCplot.case}" containing:
#' \item{power_item}{A list of length k containing ggplot objects of power plots of each item and selected test statistics and
#'  modification of critical region (cond).}
#' \item{hist_item}{Histograms of critical values of sufficient statistics \eqn{t}, without ("test name") and with (mod) modification of critical region.}
#'  \item{plotlist}{A list of length k containing combined ggplot objects of (1) power plots of each item and selected test statistics and
#'  modification of critical region, and (2) Histograms of critical values of sufficient statistics \eqn{t}, without ("test name") and with (mod) modification
#'  of critical region.}
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
#'
#' # load MC sim data of case G-V-30x18 (with N=30, k=18)
#' #  and generate power plots of selected test and
#' #  with modification of critical region
#'
#'res <- MCplot_simdata_crit(N=30, k=18, nstats = "G")
#'
#'
#' }

MCplot_simdata_crit <- function( N, k,
                                 nstats = c("W", "LR", "RS", "G", "abs", "sqs", "St"),
                                 alpha = 0.05,
                                 xlim = 5,
                                 legend_title = ggplot2::element_blank(), # "Test",
                                 tag_title = ggplot2::element_blank(),
                                 name_title = ggplot2::element_blank(),
                                 legend.position = c(.5, .75), # "top",
                                 legend.direction = "vertical",
                                 tcolor=FALSE,
                                 alpha_mod = 1,
                                 ctr_mod = c("low", "high", "both")){
  call<-match.call()
  if (length(nstats) > 1)  nstats <- c("sqs")
  if (length(ctr_mod) > 1)  ctr_mod <- c("both")
  # legend_title = "Test" # expression(italic("n"))

  obj <- load_data(N = N, k = k) # load simulation data und results

  results <- plotlist_crit( object=obj,
                            nstats=nstats,
                            alpha=alpha,
                            xlim=xlim,
                            legend_title=legend_title,
                            tag_title=tag_title,
                            name_title=name_title,
                            legend.position = legend.position,
                            legend.direction = legend.direction,
                            tcolor=tcolor,
                            alpha_mod = alpha_mod,
                            ctr_mod =ctr_mod)

  results$object <- obj
  results$call <- call

  class(results) <- "MCplot_case_crit"

  return(results)

}
