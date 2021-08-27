######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# MCMC_boot
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# This file contains a routine that computes sample distribution
# of test statistics given an input matrix and covariate vector
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 07/05/2021
######################################################################


MCMC_boot <- function(y, x, n_iter, ctr){ # n_eff) {

  # t1t <- Sys.time()  # start time
  call<-match.call()

  # ---------------------------------------------------------------------
  # check of data matrix X for
  #    - inappropriate response patterns
  #    - ill-conditioned data matrix
  # ---------------------------------------------------------------------
  XWcheck <- tcl_components(X=y)
  if (XWcheck){
    warning(paste0(
      "\n",
      "\n",
      prettyPaste("Estimation stopped due to ill-conditioned data matrix X! Suspicious items in full model):"),
      "\n",
      "\n"
      # paste("No Estimation not possible in full model!", collapse=" ")
    ),
    call. = FALSE, immediate.=TRUE)
    return ("none")
  } # end if

  k <- dim(y)[2]

  #-----------------------------------------------------------------
  # MCMC samples
  #-----------------------------------------------------------------

  r_rest <- psychotools::raschmodel(y, hessian = FALSE)
  # eta_rest <- -r_rest$coefficients #  estimated item (easiness) parameters!

  empty_list <- vector(mode = "list", length = n_iter)  # empty list with length
  names( empty_list) <- paste0("MC",1:n_iter)
  lmat2 <-  empty_list

  lstats <- matrix(0,nrow = n_iter, ncol = 7,
                   dimnames = list(paste0("MC",1:n_iter), c("W", "LR", "RS", "G", "abs", "sqs", "St") ) ) # matrix for test statistics

  tmat <- matrix(0,nrow = k, ncol = n_iter,
                 dimnames = list(paste0("I",1:k),paste0("MC",1:n_iter)) ) # matrix for test statistics

  ieff <-0 # total sample count
  imat <-1 # succesful matrix
  ivec <-1
  m <- 0

  while(imat <= n_iter) {

    mcs <- MC_samples(y=y, x=x, ctr=ctr) # n_eff = n_eff, seed = 0)
    o = mcs$o
    t <- mcs$t

    m <- m + 1 # number of MC samples

    lmat0 <- eRm::rstats(o, cbind)[-1] # list of MC matrix samples without input matrix

    for (i in 1:length(lmat0)) {

      #---------------------------------------------------------------------
      # check of data matrix X for
      #    - inappropriate response patterns within subgroups
      #    - ill-conditioned data matrix (restricted and unrestricted models)
      #---------------------------------------------------------------------
      # e <- tcl_splitcr(X = X, splitcr = x, model = "RM")
      # Xlist_check <- tcl_datcheck(X=X, Xlist = e$X.list, model = model)

        check_g <- tcl_checkdata2(X=lmat0[[i]], splitcr = x)

      if(!check_g) {
        lmat2[[imat]] <- lmat0[[i]]

        lstats[imat,1:4] <- mc_stats_n(X=lmat0[[i]], splitcr = x, r = r_rest)

        tmat[,imat] <- t[,i]
        imat <- imat + 1
        # ivec <- c(ivec,i)

        if (imat > n_iter) break # stop for-loop
      } # end if

      ieff <- ieff + 1

    } # end for i

    if (m > 100)  {
      print(paste("Coming out from While loop, number of samples drawn =  ", m))
      break
    }


  } # end while

  # value of score function in each sample drawn contained in columns
  score <- tmat - rowMeans(tmat)

  tstats <- mc_stats_t(score)
  lstats[,5] <- tstats$abs
  lstats[,6] <- tstats$sqs

  lstats[,7] <- tstats$St

  lstats <-  t(lstats)


  results <- list( eta_rest = round (-r_rest$coefficients, digits = 2), #  estimated item (easiness) parameters!
                  inputmat = y,
                  score = score ,
                  lstats = lstats,
                  t = tmat,
                  splitcr = x)

  class(results) <- "MCstats"

  return(results)
}


