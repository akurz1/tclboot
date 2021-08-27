# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 06/05/2021
######################################################################
#' Dataset \code{data.results} in \pkg{tclboot} Package
#'
#' Results of MCMC study of 60 cases of simulated Rasch data of size n x k included in the \pkg{tclboot} package.
#'
#'
#' @format
#'
#'\describe{
#' \item{lpwr_d0}{Type I error rates for different tests for each case.}
#' \item{larg_min}{Optimal input arguments of items from minimization of power function for different tests for each case.}
#' \item{pvalue}{\eqn{p} values for different tests for each case.}
#' \item{plot}{Power rates for different items and tests as a function of \eqn{d} for each case. Nominal level \eqn{\alpha} set to 0.05}
#' \item{pwr_xtable}{Tables representing power rates for different items and tests as a function of \eqn{d} for each case. Nominal level \eqn{\alpha} set to 0.05}
#' \item{descr_table}{Descriptive statistics for different tests under the null hypothesis based on
#' 10.000 bootstrap replications for each case. Nominal level \eqn{\alpha} set to 0.05}
#' }
#'
#' \code{data.results} has the following format:
#'
#'\enumerate{
#'\code{List of 6}:
#'
#'\item \code{I : List of 10}
#'
#'\itemize{
#'\item \code{$ A-I-10x6 :List of 6} \cr
#'
#'\code{$ lpwr_d0    : Named num [1:7] 0.0498 0.0448 0.0458 0.0476 0.0456 0.0398 0.0448} \cr
#'\code{$ larg_min   : num [1:7, 1:6] -5.45 -7.84 -7.78 -20.69 -7.78 ...} \cr
#'\code{$ pvalue     : Named num [1:7] 0.392 0.353 0.383 0.41 0.353 ...} \cr
#'\code{$ plot       :List of 6} \cr
#'\code{$ pwr_xtable :List of 6} \cr
#'\code{$ descr_table:List of 5} \cr
#'
#'\item \code{$ B-I-10x8 :List of 6} \cr
#'\code{ ... } \cr
#'\item \code{$ C-I-10x10:List of 6} \cr
#'\code{ ... } \cr
#'\item \code{$ D-I-10x12:List of 6} \cr
#'\code{ ... } \cr
#'\item \code{$ E-I-10x14:List of 6} \cr
#'\code{ ... } \cr
#'\item \code{$ F-I-10x16:List of 6} \cr
#'\code{ ... } \cr
#'\item \code{$ G-I-10x18:List of 6} \cr
#'\code{ ... } \cr
#'\item \code{$ H-I-10x20:List of 6} \cr
#'\code{ ... } \cr
#'\item \code{$ I-I-10x25:List of 6} \cr
#'\code{ ... } \cr
#'\item \code{$ J-I-10x30:List of 6} \cr
#'\code{ ... } \cr
#'}
#'\item \code{II : List of 10}
#'
#'
#'\item \code{III : List of 10}
#'
#'
#'\item \code{VI : List of 10}
#'
#'
#'\item \code{V : List of 10}
#'
#'
#'\item \code{VI : List of 10}
#'
#'} %% end itemize
#'Note: Roman numbers \code{I} to \code{VI} refer to cases with sample sizes \eqn{n = 10, 15, 20, 25, 30, 50}. Letters \code{A} to \code{J} refer to
#'cases with number of items \eqn{k = 6, 8, 10, 12, 14, 16, 18, 20, 25, 30}. Total 60 cases.
#'
"data.results"
