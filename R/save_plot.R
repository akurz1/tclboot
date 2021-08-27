######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
#  save_plot
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# This file contains a utility function that saves plots to disc
#
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 12/06/2021
######################################################################
#' Utility function for saving a plotlist of ggplot objects
#'
#' This utility function saves individual plots of a plotlist of ggplot objects
#'  using \pkg{ggplot2} function \code{ggsave}.
#'
#'
#' @param plotlist A list of length k containing ggplot objects of power plots.
#' @param sub_dir Sub directory to save plot to. \code{sub_dir} is added by default to current working directory.
#' @param scale Multiplicative scaling factor used in \code{ggplot}. Default = 1.
#' @param px Plot size in "mm". Default = 80 mm.
#' @param py Plot size in "mm". Default = 80 mm.
#'
#' @return A list of results.
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
#' #
#'
#' res <- MCplot_simdata_ext(npers = c(10,20, 30, 50), k=14, nstats = "sqs", xlim = 4, tcolor = TRUE)
#'
#' res <- save_plot(plotlist = res$plotlist, sub_dir = "plot_k=14")

#'
#'}


save_plot <- function (plotlist, sub_dir = "plot",scale = 1,  px = 80, py = 80) {

  call<-match.call()

 main_dir <- getwd() # "/path/to/main/dir"
 output_dir <- file.path(main_dir, sub_dir) # "/path/to/main/dir/sub_directory"
 if (!dir.exists(output_dir)) {dir.create(output_dir)}

 x = length(plotlist)
 # px <- py <- 40
 rows <- cols <- 1
 pheight <- rows * py
 pwidth <- cols * px

 for (i in seq(x)) {
   filename = paste0("item_",i,".pdf")

  ggplot2::ggsave(filename,
          plot = plotlist[[i]],
          path =  output_dir,
          scale = scale,
          width = pwidth,
          height = pheight,
          units = "mm")
 }

 return("call" = call)
}
