# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 05/05/2021
######################################################################
#' Dataset \code{data.sim.rasch} in \pkg{tclboot} Package
#'
#' Datasets of 60 cases of simulated Rasch data of size \eqn{n * k} included in the \pkg{tclboot} package.
#'
#' Bootstrap data was generated by repeated resampling with replacement from an
#' initial data set by a non-parametric bootstrap simulation design. From an initially
#' generated data sample of size \eqn{n * k}, denoted as \code{inputmat}, a number of 10.000
#' bootstrap samples are generated and test statistics are computed.
#'
#'The initial data sample of size \eqn{n * k} was generated according to the Rasch model,
#'given a vector of initial item parameters. Person parameters were drawn from N(0,1).
#'A sample had to meet the requirements of Fischer (1981) and were not allowed to have any
#'items with only zero or full responses, both in the restricted and unrestricted models.
#'
#'The initial data sample is then used as input matrix \code{inputmat} for the Rasch sampler (Mair & Hatzinger, 2007; Verhelst et al., 2007).
#'A number of \eqn{B = 10,000} bootstrap samples of size \eqn{n * k} are generated. Each bootstrap sample had to meet the
#'requirements of Fischer (1981) and were not allowed to have any items with only zero or full responses,
#'both in the restricted and unrestricted models.
#'
#'Finally, the bootstrap replications of the test statistics are obtained for
#'Wald (W), likelihood ratio (LR), Rao score (RS), gradient (GR), test statistic based on absolute values of the score function (abs) and squared
#'elements of the score function (sqs), and an alternative version of the score test (St).
#'
#' @format \code{data.sim.rasch} has the following format:
#'
#'\itemize{
#'\code{List of 5}
#'
#'\item \code{eta_list}: A list of numeric matrices of size 2 x (k-1), with item parameters as rows and items as columns. The first row
#'refers to CML estimates of item parameters of the restricted (full) model, the second row refers to initial item parameters used to generate
#'the initial input matrix for each case. Note the item parameter of first item is set to 0.
#'
#'\code{ eta_list : List of 6} \cr
#'\code{ $ I  :List of 10} \cr
#'\code{ $ II :List of 10} \cr
#'\code{ $ III:List of 10} \cr
#'\code{ $ IV :List of 10} \cr
#'\code{ $ V  :List of 10} \cr
#'\code{ $ VI :List of 10} \cr
#'
#'\item \code{inputmat_list}: A list of numeric input matrices of size n x k for each case with persons as rows, items as columns,
#'
#'\code{inputmat_list : List of 6} \cr
#'
#'\code{ $ I  :List of 10} \cr
#'\code{ A-I-10x6 : num [1:10, 1:6] 1 0 1 0 0 1 1 0 0 1 ...} \cr
#'\code{ B-I-10x8 : num [1:10, 1:8] 0 0 1 1 1 1 0 0 0 0 ...} \cr
#'\code{ C-I-10x10: num [1:10, 1:10] 0 1 1 0 1 1 1 0 0 1 ...} \cr
#'\code{ D-I-10x12: num [1:10, 1:12] 1 0 0 1 1 0 1 1 0 1 ...} \cr
#'\code{ E-I-10x14: num [1:10, 1:14] 0 0 0 1 1 0 1 0 0 1 ...} \cr
#'\code{ F-I-10x16: num [1:10, 1:16] 1 1 1 0 0 1 0 1 0 0 ...} \cr
#'\code{ G-I-10x18: num [1:10, 1:18] 1 1 1 0 0 1 0 1 0 0 ...} \cr
#'\code{ H-I-10x20: num [1:10, 1:20] 1 1 0 0 1 1 1 0 1 0 ...} \cr
#'\code{ I-I-10x25: num [1:10, 1:25] 1 0 1 0 1 1 0 1 1 0 ...} \cr
#'\code{ J-I-10x30: num [1:10, 1:30] 0 1 1 0 1 0 1 1 0 0 ...} \cr
#'
#'\code{ $ II :List of 10} \cr
#'\code{ A-II-15x6 : num [1:15, 1:6] 0 1 0 1 1 1 1 1 1 0 ...} \cr
#'\code{ ... } \cr
#'\code{ J-II-15x30: num [1:15, 1:30] 1 0 1 0 1 0 1 0 1 1 ...} \cr
#'
#'\code{ $ III:List of 10} \cr
#'\code{ A-III-20x6 : num [1:20, 1:6] 1 0 0 1 1 0 1 1 0 0 ...} \cr
#'\code{ ... } \cr
#'\code{ J-III-20x30: num [1:20, 1:30] 0 0 1 1 0 0 1 1 1 0 ...} \cr
#'
#'\code{ $ IV :List of 10} \cr
#'\code{ A-IV-25x6 : num [1:25, 1:6] 0 1 1 1 0 1 0 0 0 0 ...} \cr
#'\code{ ... } \cr
#'\code{ J-IV-25x30: num [1:25, 1:30] 0 1 1 0 1 0 1 0 0 0 ...} \cr
#'
#'\code{ $ V  :List of 10} \cr
#'\code{ A-V-30x6 : num [1:30, 1:6] 0 0 1 0 0 0 1 0 1 1 ...} \cr
#'\code{ ... } \cr
#'\code{ J-V-30x30: num [1:30, 1:30] 1 0 0 1 1 0 1 0 0 0 ...} \cr
#'
#'\code{ $ VI :List of 10} \cr
#'\code{ A-VI-50x6 : num [1:50, 1:6] 0 0 0 0 1 0 0 1 1 1 ...} \cr
#'\code{...} \cr
#'\code{ J-VI-50x30: num [1:50, 1:30] 0 1 1 0 1 0 0 0 1 1 ...} \cr
#'
#'\item \code{score_list}: A list of numeric matrices of size k x 10000 for each case, containing in each column the value of score
#'function in each sample drawn.
#'
#'\code{ score_list : List of 6} \cr
#'\code{ $ I  :List of 10} \cr
#'\code{ $ II :List of 10} \cr
#'\code{ $ III:List of 10} \cr
#'\code{ $ IV :List of 10} \cr
#'\code{ $ V  :List of 10} \cr
#'\code{ $ VI :List of 10} \cr
#'
#'\item \code{stat_list}:  A list of numeric matrices of size 7 x 10000 for each case, containing in rows values of test statistics for each sample drawn. See Details.
#'
#'\code{ stat_list : List of 6} \cr
#'\code{ $ I  :List of 10} \cr
#'\code{ $ II :List of 10} \cr
#'\code{ $ III:List of 10} \cr
#'\code{ $ IV :List of 10} \cr
#'\code{ $ V  :List of 10} \cr
#'\code{ $ VI :List of 10} \cr
#'
#'\item \code{t_list}: A list of numeric matrices of size k x 10000 for each case, containing in each column the observed values of sufficient statistic for \eqn{d} computed for each sample drawn.
#'
#'\code{ t_list : List of 6} \cr
#'\code{ $ I  :List of 10} \cr
#'\code{ $ II :List of 10} \cr
#'\code{ $ III:List of 10} \cr
#'\code{ $ IV :List of 10} \cr
#'\code{ $ V  :List of 10} \cr
#'\code{ $ VI :List of 10} \cr
#'
#'} %% end itemize
#'Note: Roman numbers \code{I} to \code{VI} refer to cases with sample sizes \eqn{n = 10, 15, 20, 25, 30, 50}. Letters \code{A} to \code{J} refer to
#'cases with number of items \eqn{k = 6, 8, 10, 12, 14, 16, 18, 20, 25, 30}. Total 60 cases.
#'
"data.sim.rasch"
