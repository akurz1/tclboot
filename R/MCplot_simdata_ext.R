######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# MCplot_simdata_ext
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# This file contains a routine that plots the power function of
# simulated data.
#
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 08/06/2021
######################################################################
#' Plot of power function.
#'
#'This function plots the power function of selected test statistics
#'
#' @param npers A vector of sample sizes \eqn{n} for cases of interest.  N =  10, 15, 20, 25, 30, 50. Default value are N = 10, 20, 30, 50.
#' @param k Number of items \eqn{k} of selected case. k = 6, 8, 10, 12, 14, 16, 18, 20, 25, 30
#' @param nstats Name of test statistics to be plotted. Note a maximum of 1 test allowed. By default "sqs" will be plotted.
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
#'  \item{plotlist}{A list of length k containing ggplot objects of power plots of each item and selected test statistics for sample sizes.}
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
#' # load MC sim data of cases
#' # with N=10,20,30,50, k=14
#' # and generate power plots of sqs test
#'
#'res <-  MCplot_simdata_ext(npers = c(10,20, 30, 50), k=14, nstats = "sqs", xlim = 4, tcolor = TRUE)
#'
#' }

MCplot_simdata_ext <- function(npers = c(10,15,20,25,30,50),
                               k,
                               nstats = c("W", "LR", "RS", "G", "abs", "sqs", "St"),
                               alpha = 0.05,
                               xlim = 5,
                               legend_title = "n",
                               tag_title = ggplot2::element_blank(),
                               name_title = ggplot2::element_blank(),
                               legend.position = c(.5, .75), # "top",
                               legend.direction = "vertical",
                               tcolor=FALSE){
  call<-match.call()
  if (length(npers) > 4) npers <- c(10,20,30,50)
  if (length(nstats) > 1)  nstats <- c("sqs")

  i_pers <- seq(npers)

  # empty_list <- vector(mode = "list", length = length(npers) + 1)
  # names(empty_list) <- c("x", as.character(npers[i_pers])  )
  empty_list <- vector(mode = "list", length = length(npers))
  names(empty_list) <- c(as.character(npers[i_pers])  )


obj_list <- empty_list

  for (ip in i_pers) { ######## N

    obj_list[[ip]] <- load_data(N = npers[ip], k = k) # load simulation data und results

  } # end for
  class(obj_list) <- "MCplot_case_ext"

  # obj <- load_data(N = N, k = k) # load simulation data und results
  plotlist <- plotlist_ext(object_list=obj_list,
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
                  object = obj_list,
                  "call" = call)

  class(results) <- "MCplot_case_ext"

  return(results)

}

