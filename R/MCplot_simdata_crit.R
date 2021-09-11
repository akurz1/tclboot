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
# copyright (c) 2021, Last Modified 30/08/2021
######################################################################
#' Plot of power function.
#'
#'This function plots the power function of selected test statistics
#'
#' @param N Sample size \eqn{n} of selected case.  N =  10, 15, 20, 25, 30, 50
#' @param k Number of items \eqn{k} of selected case. k = 6, 8, 10, 12, 14, 16, 18, 20, 25, 30
#' @param d_item A vector of length k containing delta values of each item. If missing, by default all deltas are set to 0.
#' @param nstats Name of test statistics chosen out of possible tests c("W", "LR", "RS", "G", "U1", "U2", "St") to be modified. Note only one test allowed. By default "U1" will be plotted.
#' Note "U1" = "sqs" or sum of squared elements of the score function and "U2" = "abs" or sum of absolute values.
#' @param alpha Probability of error of first kind (in plot shown as dotted horizontal line).
#' @param xlim  A numeric value, specifying the left/lower limit and the right/upper limit of the x scale. Maximum value is xlim = 10.
#' @param legend_title The text for Legend. By default it is set to "Test".
#' @param tag_title The text for tag label, which will be displayed at the top-left of the plot, by default of Item 1 only. By default it is element_blank().
#' @param name_title The text for title of plot. By default it is element_blank().
#' @param legend.position the position of legends ("none", "left", "right", "bottom", "top", or two-element numeric vector)
#' @param legend.direction layout of items in legends ("horizontal" or "vertical")
#' @param tcolor Logic value. If TRUE then ggplot plot will be in color otherwise black and white.
#' @param ctr_alpha A numeric value. The crit. region \eqn{C} of selected test statistic will be \eqn{100(alpha*ctr_alpha)} percent largest values of this test statistics. Default value is 1.
#' @param ctr_mod Character. Which type of modification of critical region should be used. Adding low most extreme values of the respective of \eqn{T} statistics ("low"),
#' adding high most extreme values of the respective of \eqn{T} statistics ("high"), or adding both ("both"). The latter is also Default value.
#'
#' @return An object of class "\code{MCplot_case_crit}" containing:
#' \item{power_item}{A list of length k containing ggplot objects of power plots of each item and selected test statistics and
#'  modification of critical region (cond).}
#' \item{hist_item}{Histograms of critical values of sufficient statistics \eqn{t}, without ("test name") and with (mod) modification of critical region.}
#'  \item{plotlist}{A list of length k containing combined ggplot objects of (1) power plots of each item and selected test statistics and
#'  modification of critical region, and (2) Histograms of critical values of sufficient statistics \eqn{t}, without ("test name") and with (mod) modification
#'  of critical region.}
#'  \item{resultlist}{A list of length k containing several sublists: values of sufficient statistics \eqn{t} "t_item",
#'  indices of \eqn{100(alpha*ctr_alpha)} percent largest values ("icrit"), and \eqn{100(alpha)} percent largest values of test statistic ("icrit0"),
#'  incices of modification ("imod", "imod_low", "imod_high"), difference of indices ("idiff"), union of indices in crit. reg. after modification ("icrit_mod"),
#'  local power values with all \eqn{delta}=0 ("lpwr_d0"), local power with vector d_item and d=0 for item i ("lpwr_d_item"), and a table of rel. frequency of values of sufficient statistics \eqn{t} in critical region ("t_table") }
#'  \item{hist_crit_reg}{A list of length k containing histograms of critical values of sufficient statistics \eqn{t}, without modification ("test name") and with (mod) modification of critical region.}
#'   \item{call}{The matched call.}
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
                                 d_item,
                                 # nstats = c("W", "LR", "RS", "G", "abs", "sqs", "St"),
                                 nstats = c("W", "LR", "RS", "G", "U1", "U2", "St"),
                                 alpha = 0.05,
                                 xlim = 5,
                                 legend_title = ggplot2::element_blank(), # "Test",
                                 tag_title = ggplot2::element_blank(),
                                 name_title = ggplot2::element_blank(),
                                 legend.position = c(.5, .75), # "top",
                                 legend.direction = "vertical",
                                 tcolor=FALSE,
                                 ctr_alpha = 1,
                                 ctr_mod = c("low", "high", "both")){
  call<-match.call()

  if (missing(d_item)) d_item <- rep(0, k) # vector d_item  equal zero
  if (d_item[length(d_item)] != 0) d_item[length(d_item)] <- 0 # kth vector element set equal 0
  names(d_item) <- paste0("I", seq(1:length(d_item)))

  if (length(nstats) > 1) {
    warning("Test for modification set by default to U1!")
    nstats <- c("U1")
  }

  if (nstats == "U1") nstats <- "sqs" # for internal use
  if (nstats == "U2") nstats <- "abs" # for internal use

  if (length(ctr_mod) > 1)  {
    warning("Default value for modification set (both)!")
    ctr_mod <- c("both")
  }
  # legend_title = "Test" # expression(italic("n"))

  obj <- load_data(N = N, k = k) # load simulation data und results
  class(obj) <- "MCplot_case_crit"

  results <- plotlist_crit( object=obj,
                            d_item = d_item, # vector d_item
                            nstats=nstats,
                            alpha=alpha,
                            xlim=xlim,
                            legend_title=legend_title,
                            tag_title=tag_title,
                            name_title=name_title,
                            legend.position = legend.position,
                            legend.direction = legend.direction,
                            tcolor=tcolor,
                            ctr_alpha = ctr_alpha,
                            ctr_mod = ctr_mod)

  results$object <- obj
  results$d_item <- d_item
  results$call <- call

  class(results) <- "MCplot_case_crit"

  return(results)

}
