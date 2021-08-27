######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# MCboot
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# This file contains a routine that computes sample distribution
# of test statistics given an input matrix and covariate vector
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 11/05/2021
######################################################################
#' Markov Chain Monte Carlo resampling method for Rasch model
#'
#' This function computes sample distribution of test statistics for hypothesis of equality of item parameters
#' between two groups of persons against a two-sided alternative that at least one item parameter differs
#' between the two groups given an input matrix and a numeric covariate vector.
#'
#' MCMC binary matrices with given margin sums are computed using function \code{rsampler} of \pkg{eRm} package.
#' Matrices fulfilling conditions of Fischer (1981) and showing an appropriate response pattern within subgroups
#' are used to compute sample distributions of test statistics. Conditions on the existence and uniqueness of
#' maximum-likelihood estimates (Fischer, 1981) are checked using graph theory implementing FORTRAN subroutine
#' \code{digraph_adj_components.f90} (implemented from Thulasiraman & Swamy, 1992, by Burkardt, 2020).
#' Likelihood based test statistics are Wald (W), likelihood ratio (LR), Rao score (RS) and
#' gradient (GR). Nonparametric test statistics are sum of squared elements of the score function and sum of
#' the absolute values of the elements of the score function and an alternative version of the score test.
#'
#' Based on sample distribution of test statistics Type I error rates and power values are computed for the
#' hypothesis to be tested and a deviation from it.
#'
#' @param inputmat A binary data matrix with \eqn{n} rows and \eqn{k} columns.
#' @param splitcr Split criterion which is a numeric vector x with length equal to number of persons and contains zeros and
#' ones. It indicates group membership for every person.
#' @param n_iter Number of generated matrices (effective matrices) complying to the Rasch model.
#' @param ctr An object of \pkg{eRm} class \code{\link[eRm]{RSctr}} with components \code{burn_in}, \code{n_eff}, \code{step}, \code{seed}
#' used as \code{controls} argument in \code{\link[eRm]{rsampler}}.
# #' @param n_eff Number of generated matrices (effective matrices) in function \code{rsampler} of \pkg{eRm} package. By default
# #' \code{n_eff} is set to 8000. Note \code{n_eff} must be positive and not larger than 8191. See Details \code{rsctrl} \{\pkg{eRm}\}.
#' @param alpha Probability of error of first kind.
#'
#' @return \code{MCboot} returns an object of class \code{"MCboot"} containing two lists:
#' \item{MCobject}{A list of results containing:}
#' \itemize{
#'  \item \code{$eta_rest}: A numeric vector containing CML estimates of item parameters of data matrix \code{X} (full model).
#'  \item \code{$inputmat}: The input matrix.
#'  \item \code{$score}: A numeric matrix of size k x \code{n_iter}, containing in each column the value of score function in each sample drawn.
#'  \item \code{$lstats}: A numeric matrix of size k x \code{n_iter}, containing in each column values of test statistics in each sample drawn. See details.
#'  \item \code{$t}: A numeric matrix of size k x \code{n_iter}, containing in each column in each column the observed values of sufficient statistic
#'  for \eqn{d} computed for each sample drawn.
#'  \item \code{$splitcr}: Split criterion which is a numeric vector x with length equal to number of persons and contains zeros and
#' ones. It indicates group membership for every person.
#'  }
#' \item{result_list}{A list of results containing:}
#' \itemize{
#'  \item \code{$lpwr_d0}: Type I error rates for different tests.
#'  \item \code{$larg_min}: Optimal input arguments of items from minimization of power function for different tests.
#'  \item \code{$pvalue}: \eqn{p} values for different tests.
#'  \item \code{$plot}: A list for different items containing power rates of different tests as a function of \eqn{d}.
#'  \item \code{$pwr_xtable}: Tables representing power rates for different items and tests as a function of \eqn{d}.
#'  \item \code{$descr_table}: Descriptive statistics for different tests under the null hypothesis based on \code{n_iter} bootstrap replications.
#' }
#'   \item{call}{The matched call.}
#'
#' @references{
#' Burkardt, John. (2020, November). GRAFPACK - Graph Computations [FORTRAN90 library]. Retrieved from
#' \url{https://people.sc.fsu.edu/~jburkardt/f_src/grafpack/grafpack.html}
#'
#' Draxler, C., Kurz, A., & Lemonte, A. J. (2020). The gradient test and its finite sample size properties in
#' a conditional maximum likelihood and psychometric modeling context. Communications in Statistics-Simulation and
#' Computation, 1-19. \url{https://doi.org/10.1080/03610918.2019.1710193}
#'
#' Draxler, C., & Dahm, S. (2020). Conditional or Pseudo Exact Tests with an Application in the Context of Modeling
#' Response Times. Psych, 2(4), 198-208. \url{https://doi.org/10.3390/psych2040017}
#'
#' Fischer, G. H. (1981). On the existence and uniqueness of maximum-likelihood estimates in the Rasch model.
#' Psychometrika, 46(1), 59-77.
#'
#' Thulasiraman, K., & Swamy, M. N. (2011). Graphs: Theory and Algorithms. New York, John Wiley & Sons.
#'
#'  }
#' @keywords mcmc
#' @export
#' @seealso \code{\link{MCplot}}, \code{\link{MCsimRasch}}
#' @examples
#' \dontrun{
#' # Choose input matrix of size 20 x 10 from dataset data.sim.rasch
#' y <- data.sim.rasch$inputmat_list$III$`C-III-20x10`
#'
#' res <- MCboot(inputmat=y, splitcr = c(rep(0,10), rep(1,10)), n_iter = 1000)
#'
#' # Generation of initial data sample
#' items <- c(0, -2, -1, 0, 1, 2)
#' x <- c( rep(0, 10), rep(1,10))
#'
#' y <- MCsimRasch(N = 20, splitcr = x, items = items)$X
#'
#' res2 <- MCboot(inputmat=y, splitcr = x, n_iter = 1000)
#'
#' }


MCboot <- function(inputmat,
                   splitcr,
                   n_iter,
                   # n_eff = 8000,
                   ctr = eRm::rsctrl(burn_in=100, n_eff=8000, step=32, seed=0),
                   alpha = 0.05) {

  t1t <- Sys.time()  # start time
  call<-match.call()

  y <- inputmat # for internal computation

  N <- dim(y)[1]
  k <- dim(y)[2]

  if (missing(splitcr)) {
    ng1 <- floor(N/2)  # ng2 <- N - ng1
    splitcr <- c(rep(0, ng1), rep(1, N - ng1)) # binary covariate
  }
  x <- splitcr # for internal computation

  #-----------------------------------------------------------------
  # MCMC samples
  #-----------------------------------------------------------------
  mc_stats <- MCMC_boot(y=y, x=x, n_iter=n_iter, ctr=ctr) # n_eff=n_eff)


  #-----------------------------------------------------------------
  # power calculation
  #-----------------------------------------------------------------
  lpower <- pwr_calc_t0b(object = mc_stats, alpha = alpha)


  t2t <- Sys.time() # end time
  tdiff <- round( t2t - t1t, 2) # time difference = End - start
  print (tdiff)


  results <- list(MCobject = mc_stats,
                  result_list = lpower,
                  "call" = call)

  class(results) <- "MCboot"

  return(results)
}


