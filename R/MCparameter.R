######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# MCparameter
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# This file contains a routine that loads MC model parameters
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 14/05/2021
######################################################################
#' MC model parameter
#'
#' Utilty function to load MC model parameter
#'
#' @param N Sample size \eqn{n} of selected case.  N =  10, 15, 20, 25, 30, 50
#' @param k Number of items \eqn{k} of selected case. k = 6, 8, 10, 12, 14, 16, 18, 20, 25, 30
#' @param t_inputmat Logical value. If TRUE an initial data sample according to
#' the Rasch model will be generated using function \code{MCsimRasch}.

#' @return List of MC model parameters and initial data sample..
#' \item{N}{Sample size.}
#' \item{k}{Number of items.}
#' \item{splitcr}{Split criterion which is a numeric vector x with length equal
#' to number of persons and contains zeros and ones. It indicates group membership for every person.}
#' \item{items}{A vector of item parameters of length k. First item set to 0.}
#' \item{filename}{Name of MC case. Note: Roman numbers \code{I} to \code{VI}
#' refer to cases with sample sizes \eqn{n = 10, 15, 20, 25, 30, 50}.
#' Letters \code{A} to \code{J} refer to cases with number of items
#' \eqn{k = 6, 8, 10, 12, 14, 16, 18, 20, 25, 30} }
#' \item{inputmat}{If \code{t_inputmat = TRUE}: Initial data sample according to the Rasch model, generated with \code{\link{MCsimRasch}}.}
#' \item{call}{The matched call.}
#'
#' @references{
#' see also dataset \code{data.sim.rasch}
#'  }
#' @keywords mcmc
#' @export
#' @examples
#' \dontrun{
#' # Load MC model parameters and generate a new intial data matrix
#'
#' res <- MCparameter(N=10, k=10, t_inputmat = TRUE)
#'
#' # $N
#'# [1] 10
#'#
#'# $k
#'# [1] 10
#'#
#'# $splitcr
#'# [1] 0 0 0 0 0 1 1 1 1 1
#'#
#'# $items
#'# I1   I2   I3   I4   I5   I6   I7   I8   I9  I10
#'# 0.0 -1.0 -0.5  0.0  0.0  0.0  0.0  0.0  0.5  1.0
#'#
#'# $filename
#'# [1] "C-I-10x10"
#'#
#'# $inputmat
#'#     I1 I2 I3 I4 I5 I6 I7 I8 I9 I10
#'# P1   0  1  1  0  1  1  1  0  0   0
#'# P2   1  1  1  1  1  1  1  1  1   1
#'# P3   1  1  1  0  0  0  0  0  0   1
#'# P4   1  0  0  1  0  1  1  1  1   1
#'# P5   1  0  0  1  1  0  1  0  0   0
#'# P6   0  0  0  0  0  0  0  1  1   0
#'# P7   0  0  1  1  1  0  0  1  1   1
#'# P8   0  0  0  0  0  0  0  0  0   1
#'# P9   1  1  0  1  1  1  1  1  1   1
#'# P10  0  1  0  0  0  0  1  0  0   1
#'#
#'# $call
#'# MCparameter(N = 10, k = 10, t_inputmat = TRUE)
#'#
#'# attr(,"class")
#'# [1] "MCparameter"
#'
#'
#' }

MCparameter <- function(N = c(10, 15, 20, 25, 30, 50),
                        k = c(6, 8, 10, 12, 14, 16, 18, 20, 25, 30),
                        t_inputmat = FALSE) {

  call<-match.call()

  if (length(N) > 1) N[1]
  if (length(k) > 1) k[1]

  ng1 <- floor(N/2)  # ng2 <- N - ng1
  splitcr <- c(rep(0, ng1), rep(1, N - ng1)) # binary covariate as used in MC simulation study

  par <- mc_parameter() # model parameters
  npers <- unlist(par$N_list)
  itemlist <- par$item_list4
  kitem <- sapply(par$item_list4, length)

  i_pers <- which(npers == N)
  j_item <- which(kitem == k)
  items <- unlist(itemlist[j_item])
  names(items) <- paste0("I", seq(k))

  if (t_inputmat) {
    y <- MCsimRasch(N = N, splitcr=splitcr, items = items)$X
  } else {
    y <- NA
  }

  maincase <- paste (utils::as.roman(1:length (npers)) )
  subcase <- paste ( LETTERS[1 : length(itemlist) ] )

  filename <-   paste0(as.character (subcase[j_item]), "-", maincase[i_pers] ,"-",
                     as.character(N), "x",
                     as.character(k))

  results <- list(N=N, k=k, splitcr=splitcr, items=items, filename=filename, inputmat=y, "call" = call)

  class(results) <- "MCparameter"

  return(results)
}
