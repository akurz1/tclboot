######################################################################
# UMIT - Private University for Health Sciences,
#        Medical Informatics and Technology
#        Institute of Psychology
#        Statistics and Psychometrics Working Group
#
# R utility functions
#
# Part of R/tlcboot - Conditional inference in small sample scenarios
#
# Licensed under the GNU General Public License Version 3 (June 2007)
# copyright (c) 2021, Last Modified 06/07/2021
######################################################################

### R Code

#require(eRm); require(psychotools); require(Matrix)
#-----------------------------------------------------------------
# Model parameter for generation of intitial data matrix
#-----------------------------------------------------------------
mc_parameter <- function () {
  N_list <- list( 10,15,20,25, 30,50)
  item_list4 <- list(
    k6 = c(0, -1, -0.5, 0, 0.5 , 1), k8 = c(0, -1,-0.5, rep(0,3), 0.5,1),
    k10 = c(0, -1,-0.5,rep(0,5), 0.5,1),  k12 = c(0, -1,-0.5,rep(0,7), 0.5,1),
    k14 = c(0, -1,-0.5,rep(0,9), 0.5,1),  k16 = c(0, -1,-0.5,rep(0,11), 0.5,1),
    k18 = c(0, -1,-0.5,rep(0,13), 0.5,1), k20 = c(0, -1,-0.5,rep(0,15), 0.5,1),
    k25 = c(0, -1,-0.5,rep(0,20), 0.5,1), k30 = c(0, -1,-0.5,rep(0,25), 0.5,1) )
  k <- sapply(item_list4, length)
  return (list("N_list" = N_list, "item_list4" = item_list4, "k"=k, "df"=k-1 ))
}

#-----------------------------------------------------------------
#  Generation of initial input matrix Y_0
#-----------------------------------------------------------------
mc_inputmat <- function(splitcr, beta) {
  # ill_conditioned <- TRUE
  # del.pos.g <- TRUE
  check_g <- TRUE
  i.total <- 0

  groups <- as.list(sort(unique(splitcr)))
  ng1 <- length(splitcr[splitcr == groups[[1]]])
  ng2 <- length(splitcr[splitcr == groups[[2]]])

  while(check_g) {
    r.sim <- mc_sim_data2(ng1 = ng1, ng2 = ng2, b = beta)
    y <- r.sim$y
    # x <- r.sim$x  # binary covariate length N

    #-----------------------------------------------------------------
    # check of X for ill-conditioned and 0/full responses
    #-----------------------------------------------------------------
    check_g <- tcl_checkdata2(X=y, splitcr = splitcr)

    i.total <- i.total + 1
    if (i.total > 10000)  {
      print(paste("Coming out from While loop number iterations =  ", i.total))
      y <- NA
      # break
      results <- list(X=NA,splitcr=splitcr, check=check_g, i.total=i.total)
      return(results)
    }

  } # end while

  cat("\n","Iteration", i.total,"\n")
  results <- list(X=y, splitcr=splitcr, check =check_g, i_total=i.total)
  return(results)
}



mc_sim_data2 <- function(ng1, ng2, b) {
  d <- rep(0,length(b))  # data generated under null hypothesis, all delta = 0
  n <- ng1 + ng2
  t <- stats::rnorm(n)  # t <- c(rnorm(ng1), rnorm(ng2))
  k <- length(b)
  y <- matrix(0, nrow = n, ncol = k, dimnames = list(paste0("P", seq_len(n)),paste("I",1:k,sep="") ) )
  x <- c(rep(0, ng1), rep(1, ng2)) # binary covariate
# Simulation of Rasch homogeneous data
  for (i in seq_len(k) ) {
    y[,i] <- stats::rbinom(n, size = 1, prob = exp(t + b[i] + d[i] * x) / (1 + exp(t + b[i] + d[i] * x)) )
  }
  return(list (y=y, x=x))
}

#-----------------------------------------------------------------
# generate binary random matrices based on an input matrix = data
#-----------------------------------------------------------------
MC_samples <- function(y, x, ctr){ # n_eff=8000, seed=0) {
  n_eff <- ctr$n_eff
  # burn_in <- 100
  # step <- 32
  # ctr <- eRm::rsctrl(burn_in, n_eff, step, seed)  # Arguments for rsampler
  o <- eRm::rsampler(inpmat = y, controls = ctr) # draws samples

# function to compute observed value of sufficient statistic for d
  f <- function(y) colSums(y * x)

# list of observed values of sufficient statistic for d computed for each sample drawn
  list <- eRm::rstats(o, f)

# t numerical matrix containing in each column the observed values of sufficient statistic for d for each sample drawn
  t <- sapply(list, FUN = cbind)[,2:(n_eff+1)] # t without input matrix
  rownames(t) <- paste("I",1:nrow(t),sep="")

  e <- rowMeans(t)
  # value of score function in each sample drawn contained in columns
  score <- t - e
  return(list(o=o, t=t, e=e, score=score))
} # end MC_samples

#-----------------------------------------------------------------
# Likelihood Test statistics W, LR, RS, G
#-----------------------------------------------------------------
mc_stats_n <- function(X, splitcr, r)  {
  # adopted from mc-modWald.n.R
  y <- X
  x <- splitcr
  rvind <- splitcr
  Xlist <- by(X, rvind, function(x) x)
  names(Xlist) <- as.list(sort(unique(splitcr)))
  y1 <- Xlist[[2]]  # x==1
  y2 <- Xlist[[1]]  # x==0
# estimation of MLE itemparameter
  if(missing(r)) r <- psychotools::raschmodel(y, hessian = FALSE)
  r1 <- psychotools::raschmodel (y1, hessian = FALSE)
  r2 <- psychotools::raschmodel (y2, hessian = FALSE)
  b  <- r$coefficients; b <- c(0,-b)
  b1 <- r1$coefficients;  b1 <- c(0,-b1)
  b2 <- r2$coefficients;  b2 <- c(0,-b2)
  t1 <- table(factor(rowSums(y1),levels=1:(ncol(y1)-1)))
  t2 <- table(factor(rowSums(y2),levels=1:(ncol(y2)-1)))
  itemparameter <- rbind(b, b)
  itemparameter.w2 <- rbind(b1 - b1[ncol(y * x)], b2 - b2[ncol(y * abs(x - 1))])
  gruppen <- rbind(t1, t2)
# Fisher information matrix
  i <- informationsmatrix(itemparameter, gruppen)  # Score Test and mod. Wald Test
  i.w2 <- informationsmatrix(itemparameter.w2, gruppen) # classical Wald Test
# classical Wald Test:
  # SW = t( (Theta_Dach - Theta_0) ) K(Theta_Dach) (Theta_Dach - Theta_0)
  d <- (b1 - b1[ncol(y1)]) - (b2 - b2[ncol(y2)])
  w2 <- sum ( colSums ( ( (d[1:(ncol(y) - 1)] ) )
                        * solve ( solve(i.w2)[1:(ncol(y) - 1), 1:(ncol(y) - 1) ]
                                  + solve(i.w2)[ncol(y):(2 * (ncol(y)-1) ), ncol(y):(2 * (ncol(y)-1))] ) )
              * (d[1:(ncol(y)-1)] ) )
# likelihood ratio test: SLR = 2 [ l(Theta_Dach) - l(Theta_0) ]
  lrt <- -2 * ( r$loglik - ( r1$loglik + r2$loglik ) )
# Score Vektor: U(Theta_0) = dl(Theta) / dTheta
  g <- psychotools::elementary_symmetric_functions ( par = -b, order = 1) $'0'
  dg <- psychotools::elementary_symmetric_functions ( par = -b, order = 1 )$'1'
  scorevector <- colSums(y * x) - colSums((dg/g) * c(table(factor(rowSums(y * x), levels = 0:ncol(y) ) ) ) )
# Score test:  SR = t(U(Theta_0)) K(Theta_0) (U(Theta_0)
  scoretest <- sum(colSums ( scorevector [1:(ncol(y) - 1) ] *
                               (solve (i)[1:(ncol(y) - 1), 1:(ncol(y) - 1) ]
                                + solve(i)[ncol(y):(2 * (ncol(y) - 1) ), ncol(y):(2 * (ncol(y) - 1) ) ] ) )
                   * scorevector[1:(ncol(y) - 1)])
# Gradient test: ST = t ( U(Theta_0) ) (Theta_Dach - Theta_0)
  gtest <- sum ( scorevector * (b1 - b2))
  results <- c( W=w2, LR=lrt, RS=scoretest, G=gtest )
  return(results)
}

#-----------------------------------------------------------------
#  Computation of Informationsmatrix
#-----------------------------------------------------------------
informationsmatrix <- function(itemparameter, gruppen){
  d <- dim(itemparameter)
  g <- d[1] # Gruppenanzahl
  k <- d[2] - 1 # Itemanzahl
  l <- d[2]
  ma <- matrix(rep(0, l * l), ncol = l)
  info <- matrix(rep(0, k * k), ncol = k)
  t <- fisherGruppen <- list(); items <- 1:k
  bedingteWi <- function(itemparameter){
    for (i in 1:l) ma[ ,i] <- psychotools::elementary_symmetric_functions(par = itemparameter, order = 1)[[1]][-1]
    return(psychotools::elementary_symmetric_functions(par = itemparameter, order = 1)[[2]][-1, ] / ma)
  }
  diagonalen <- function(itemparameter, gruppen, gruppe){
    #Berechnen von pi
    bedingte <- bedingteWi(as.vector(itemparameter[gruppe, ]))[-l, -l]
    #Berechnen aller pi und pij
    for (i in 1:k){
      t[[i]] <- (psychotools::elementary_symmetric_functions(par = as.vector(itemparameter[gruppe, ]), order = 2)[[3]][,,items[i]][-1, ] /
                   psychotools::elementary_symmetric_functions(par = as.vector(itemparameter[gruppe, ]), order = 1)[[1]][-1])[-l, -l]
    }
    p <- t
    #Berechnen der Diagonale und Nebendiagonalen
    for (i in 1:k){
      p[[i]][ ,items[i]]  <- t[[i]][ ,items[i]] * (1 - t[[i]][ ,items[i]]) * gruppen[gruppe, ]
      p[[i]][ ,items[-i]] <- (t[[i]][ ,items[-i]] - bedingte[ ,items[i]] * bedingte[ ,items[-i]]) * gruppen[gruppe, ]
      info[i, ] <- colSums(p[[i]])
    }
    return(info)
  }
  #Anwenden auf alle Gruppen
  for (j in 1:g){
    gruppe <- j
    fisherGruppen[[j]] <- diagonalen(itemparameter = itemparameter, gruppen = gruppen, gruppe = gruppe)
  }
  #Ausgabe der Blockmatrix
  return(Matrix::bdiag(fisherGruppen))
}

#-----------------------------------------------------------------
# Test statistics abs, sqs, St
#-----------------------------------------------------------------
mc_stats_t <- function(score) {
  # score <- t - rowMeans(t)
  k <- dim(score)[1]
  n_eff <- dim(score)[2]
# test statistic based on absolute values of the score computed for each sample drawn
  abs <- colSums(abs(score))
# squared test statistic computed for each sample drawn
  sqs <- colSums(score^2)
  #-----------------------------------------------------------------
  # Variance Covariance matrix from score vectors
  #-----------------------------------------------------------------
  lcov <- vcov_t(score)
  varcov <- lcov$varcov
  t_vcov <- lcov$is.singular
  St <- rep(NA, n_eff)
  # Check if varcov matrix is sigular
  if (!lcov$is.singular) {
    COV <- lcov$varcov
    for (i in seq(n_eff)) {
      St[i] <- t(score[1:(k-1), i]) %*% solve(COV) %*% score[1:(k-1), i]
    }
  } # end if
  return(list(abs = abs, sqs = sqs, St = St, lcov = lcov))
}

#-----------------------------------------------------------------
# Variance Covariance matrix from score vectors
#-----------------------------------------------------------------
vcov_t <- function(score) {
  # score <- t - rowMeans(t)
  k <- dim(score)[1]
  # Svcov <- score[1:(k-1),] # one item not free
  varcov <- matrix(0, (k-1), (k-1))
  for (j in seq(1 : (k-1))) {
    for (l in seq(1: (k-1))) {

      if (j==l) {
        # varcov[j,j] <- mean(Svcov[j,] * Svcov[j,])   # Var(Sj)
        varcov[j,j] <- mean(score[j,] * score[l,])   # Var(Sjj)
      } else if (j!= l) {
        varcov[j, l] <- mean (score[j,] * score[l,]) # COV(Sij)
      }
    } # end for l
  } # end for j
  # Check if varcov matrix is singular
  # require(matrixcalc)
  # tvcov <- is.singular.matrix(varcov)
  f <- function(m) "matrix" %in% class(try(solve(m),silent=TRUE))
  tvcov <- !f(varcov)
  return (list("varcov"=varcov, "is.singular"=tvcov) )
}

#-----------------------------------------------------------------
# Check for condtions according to Fischer (1981)
# based on GRAFPACK (Burkardt , 2020)
#-----------------------------------------------------------------
tcl_checkdata2 <- function(X, splitcr){
  # source('digraph_adj_components.R')
  # tcl_components <- function(X){
  #   k <- ncol(X)
  #   # build adjacency matrix
  #   adj <- matrix(0, ncol = k, nrow = k)
  #   for (i in 1:k) for(j in 1:k) {
  #     adj[i,j] <- 1*any(X[,i]>X[,j],na.rm=TRUE)
  #   }
  #   # number of components according to graph theory
  #   ncomp <- digraph_adj_components(adj)$ncomp
  #   ill_conditioned <- FALSE
  #   if (ncomp > 1) ill_conditioned <- TRUE
  #   csum <- colSums(X)
  #   del.pos <- FALSE
  #   if ( (any(csum == nrow(X)) || (any(csum == 0))) ) del.pos <- TRUE
  #   check_g <- FALSE
  #   if (del.pos || ill_conditioned) check_g <- TRUE
  #   return(check_g)
  # }
  rvind <- splitcr
  Xlist <- by(X, rvind, function(x) x)
  names(Xlist) <- as.list(sort(unique(splitcr)))
  check_g <- FALSE
  if( tcl_components(X=X) ||  tcl_components(X=Xlist[[1]]) ||  tcl_components(X=Xlist[[2]]) ) check_g <-TRUE
  return(check_g)
}

tcl_components <- function(X){
  k <- ncol(X)
  # build adjacency matrix
  adj <- matrix(0, ncol = k, nrow = k)
  for (i in 1:k) for(j in 1:k) {
    adj[i,j] <- 1*any(X[,i]>X[,j],na.rm=TRUE)
  }
  # number of components according to graph theory
  ncomp <- digraph_adj_components(adj)$ncomp
  ill_conditioned <- FALSE
  if (ncomp > 1) ill_conditioned <- TRUE

  csum <- colSums(X)
  del.pos <- FALSE
  if ( (any(csum == nrow(X)) || (any(csum == 0))) ) del.pos <- TRUE

  check_g <- FALSE
  if (del.pos || ill_conditioned) check_g <- TRUE
  return(check_g)
}

#-----------------------------------------------------------------
# R wrapper for FORTRAN90 subroutine DIGRAPH_ADJ_COMPONENTS
# from GRAFPACK (Burkardt, 2020)
#------- ----------------------------------------------------------
# Version digraph_adj_components_V2.f90 with ier added
digraph_adj_components <- function(adj) {
  ### R CMD SHLIB -o digraph_adj_components.so digraph_adj_components.f90
  lda <- dim(adj)[1];  nnode <- dim(adj)[2]
  ncomp <-0; comp <- rep(0,nnode)
  dad <- rep(0,nnode); order <- rep(0,nnode)
  ier<-0
  # dyn.load("digraph_adj_components.so")
  r <- .Fortran("digraph_adj_components", adj=as.integer(adj),
                lda = as.integer(lda),
                nnode = as.integer(nnode),
                ncomp = as.integer(ncomp),
                comp = as.integer(comp),
                dad = as.integer(dad),
                order = as.integer(order),
                ier=as.integer(ier))
  if (r$ier == 1) r <- NA
  return(r)
}

#-----------------------------------------------------------------
# p-value, tables
#-----------------------------------------------------------------
mc_pvalue_test <- function(y, x, lstats, t, score) {

  N <- dim(y)[1]
  # ng1 <- floor(N/2) # subgroup 0
  # ng2 <- N - ng1 # subgroup 1
  k <- dim(y)[2]
  # x <- c(rep(0, ng1), rep(1, ng2)) # binary covariate
  x = x

  ## calculate test statistics of input matrix y
  empty_vec <- vector(mode = "double", length = 7)  # vector of test statistics of input matrix y
  names(empty_vec) <- c("W", "LR", "RS", "G", "abs", "sqs", "St")

  stats_y <-pvalue <- empty_vec

  stats_y[1:4] <-  mc_stats_n(X=y, splitcr = x) # "W", "LR", "RS", "G",


  # observed values of sufficient statistic for d for input matrix y
  t_y <- colSums(y * x)
  score_y <- as.matrix(t_y - rowMeans(t))

  # test statistic based on absolute values of the score computed for input matrix
  stats_y[5] <- colSums(abs(score_y))

  # squared test statistic computed for input matrix
  stats_y[6] <- sum(score_y^2)

  # stats_y[7] <- NA
  lcov <- vcov_t(score)

  if (!lcov$is.singular) {
    COV <- lcov$varcov
    stats_y[7] <- t(score_y[1:(k-1), 1]) %*% solve(COV) %*% score_y[1:(k-1), 1]
  } else {
    stats_y[7] <- NA
  } # end if

  stats_tot <- cbind(MC0=stats_y,lstats)

  pvalue_test <- function(stats) {
    # i_stats_y <- stats_tot[,1]
    # i_stats <- stats_tot[,2:ncol(stats_tot)]

    n_iter <- length(stats) - 1
    pv <- table(stats[1] < stats)[2] / n_iter # p value of test 1 using t1
    return (pv)
  }

  for (i in 1:7) {

    if (i < 7) {
      pvalue[i] <- pvalue_test(stats =  stats_tot[i,])
    } else if (i == 7) {
      if (any(is.na(stats_tot[i,])) ) {
        pvalue[i] <- NA
      } else {
        pvalue[i] <- pvalue_test(stats = stats_tot[i,])
      }
    }

  }

  return(pvalue)
}


#-----------------------------------------------------------------
# power calculation
#-----------------------------------------------------------------
# power <- function(d) sum(exp(colSums(t[,stats > quantile(stats, pquant)] * d))) / sum(exp(colSums(t * d)))

# utility functions
pwr_d0_stats <- function(stats,k,t,pquant) {
  stats <-stats
  k <- k
  t <- t
  pquant <- pquant

  # if( any(is.na(stats)) ) return( list(pwr_d0 = NA, arg_min = NA) )
  if( any(is.na(stats)) ) return( NA)

  power <- function(d) sum(exp(colSums(t[,stats > stats::quantile(stats, pquant)] * d))) / sum(exp(colSums(t * d)))

  pwr_d0 <- power(d = rep(0, k))
  #opt_d0 <- nlm(power, rep(0, k))

  #arg_min <- opt_d0$estimate - opt_d0$estimate[k]
  #arg_min <- round(arg_min, digits = 3)

  # results <- list(pwr_d0 = pwr_d0, arg_min = arg_min)
  # return(results)
  return(pwr_d0)
}

pwr_arg_min_stats <- function(stats,k,t,pquant) {
  stats <-stats
  k <- k
  t <- t
  pquant <- pquant

  if( any(is.na(stats)) ) return( rep(NA, k) )

  power <- function(d) sum(exp(colSums(t[,stats > stats::quantile(stats, pquant)] * d))) / sum(exp(colSums(t * d)))

  #pwr_d0 <- power(d = rep(0, k))
  opt_d0 <- stats::nlm(power, rep(0, k))

  arg_min <- opt_d0$estimate - opt_d0$estimate[k]
  arg_min <- round(arg_min, digits = 3)

  # results <- list(pwr_d0 = pwr_d0, arg_min = arg_min)
  # return(results)
  return(arg_min)
}


plot_pwr_item <- function(item, stats, icrit, k, t, pquant, xlim = 10) {
  item <- item
  stats<-stats
  k <- k
  t <- t
  pquant <- pquant  # 1 - alpha

  d_item <- rep(0, k) # vector d_item  equal zero

  # if (missing(icrit)) icrit <- which(stats > quantile(stats, pquant) )   # index for stats
  if (missing(icrit)) icrit <- which(stats > stats::quantile(stats, pquant, na.rm = TRUE) )   # index for stats

  power_i <- function(d, item, icrit) {
    # d_item[item] <- d  # set ith d_item to d
    # return(sum(exp(colSums( t[, icrit] * d_item))) / sum(exp( colSums( t * d_item ))))
    pwr <- sum(exp( t[item, icrit] * d )) / sum(exp(t[item,] * d ))
    if (pwr == 0) pwr <- NA # case na.rm = TRUE
    return(pwr)
  }

  # calcultate power curve
  # dvec<- seq(-5, 5, by=.1)
 # dvec<- seq(-10, 10, by=.1)
  dvec <- seq(-xlim, xlim, by=.1)
  pwr <- sapply(dvec, FUN = power_i, item=item, icrit=icrit)
  pwr_mat <- as.data.frame(cbind(dvec, pwr))

  dvec<- seq(-1.75, 1.75, by=.25)
  pwr <- sapply(dvec, FUN = power_i, item=item, icrit=icrit)
  pwr_table <-  as.data.frame(cbind(dvec, pwr))

  results <- list(pwr_mat=pwr_mat,
                  pwr_table=pwr_table)
  return(results)
}

# power_i <- function(d, item, icrit) {
#   #d_item[item] <- d  # set ith d_item to d
#   pwr <- sum(exp( t[item, icrit] * d )) / sum(exp(t[item,] * d ))
#   return(pwr)
# }

############# Descriptive statistics ##########
### with Coefficient of Variation (CV) = SD / MEAN
descr_calc2 <- function(lstats, k) {

  df <- k - 1
  name.test <- rownames(lstats)  # c("W", "LR", "RS", "G", "abs", "sqs", "St")

  # values of chi^2 limiting distribution with df
  chi_mean <- df
  chi_var <- 2*df
  chi_cv <- sqrt(chi_var) / chi_mean
  chi_skew <- 2 * sqrt(2) / sqrt(df)
  chi_kurt <- 12 / df

  chi_limit <- c(mean=chi_mean, var=chi_var, cv = chi_cv,
                 skew=chi_skew, kurt=chi_kurt,
                 q50=stats::qchisq(.5,df), q90=stats::qchisq(.9,df),
                 q95=stats::qchisq(.95,df),q99=stats::qchisq(.99,df))
  names(chi_limit) <- "chi"

  # descriptive statistics of test statistics
  r_mean <- apply(lstats,1,mean)
  r_var <- apply(lstats,1,stats::var)
  r_cv <- sqrt(r_var) / r_mean
  r_skew <- apply(lstats,1, psych::skew)
  r_kurtosis <- apply(lstats,1, psych::kurtosi)

  r_q50 <- apply(lstats,1,stats::quantile, probs = 0.5, na.rm = TRUE)  # q50 = median
  r_q90 <- apply(lstats,1,stats::quantile, probs = 0.9, na.rm = TRUE)  # q90
  r_q95 <- apply(lstats,1,stats::quantile, probs = 0.95, na.rm = TRUE)  # q95
  r_q99 <- apply(lstats,1,stats::quantile, probs = 0.99, na.rm = TRUE)  # q99

  d_stats <- as.data.frame(rbind(chi=chi_limit, cbind(r_mean, r_var, r_cv,r_skew, r_kurtosis, r_q50, r_q90, r_q95, r_q99) ) )
  d_stats <- round(d_stats, digits = 3)
  colnames(d_stats) <- c("mean", "var","cv",  "skew", "kurt", "q50", "q90","q95","q99")

  # (d_stats-chi_limit)/chi_limit*100

  deviation <- function(x) {
    dev <- (x[2:length(x)]-x[1]) / x[1]*100
    return(c(NA,dev))
  }

  # deviation from chi^2 limiting distribution with df = k-1
  # dev <- apply(d_stats,2, FUN = deviation)
  dev <- as.data.frame( round(apply(d_stats,2, FUN = deviation), digits = 2 ) )
  colnames(dev) <- paste0(colnames(d_stats),"_%")

  d_stats_dev <- cbind(d_stats, dev)
  d_stats_dev <- d_stats_dev[,c(1,10,2,11,3,12,4,13,5,14,6,15,7,16,8,17,9, 18)]

  results <- list(mean = r_mean,
                  var = r_var,
                  # d_stats = d_stats,
                  # dev = dev,
                  d_stats_dev = d_stats_dev)

  return(results)
} # end descr_calc2


#############
pwr_calc_t0b <- function(object, alpha, tcond_item = FALSE) {
# pwr_calc_t0 <- function(y, lstats,k, t, score, alpha, alpha_cond, istats = c(1), tcond_item = FALSE) {
  # pwr_calc_t0 <- function(lstats,k, t, score, pquant, istats = c(1,2,3,4)) {

  y <- object$inputmat
 # N <- dim(y)[1]
 # k <- dim(y)[2]
  x <- object$splitcr

  lstats <- object$lstats   # lstats
  k <- dim(y)[2] # k <- k
  t <- object$t
  score <- object$score
  # pquant <- 1 - alpha

  name.test <- rownames(lstats)  # c("W", "LR", "RS", "G", "abs", "sqs", "St")
  # istats <- seq(length(name.test)) # NEW


  # # empty_list <- vector(mode = "list", length = length (istats))  # empty list with length
  # empty_list <- vector(mode = "list", length = length (name.test))  # empty list with length
  # names(empty_list) <- name.test #[istats]
  # lpwr_d0_stats <- larg_min <- empty_list
  # # lpwr_d0 <- lpwr_d0_stats <- larg_min <- lpwr <- empty_list


  ### check St for any NA, del index for St in istats vector ###
  # if (any(istats==7)) {
  #   if (any(is.na(lstats[7,])) ) istats <- istats[-which(istats==7)]
  # }
  ##############################################################

  #-----------------------------------------------------------------
  # calculation of pvalues
  #-----------------------------------------------------------------
  pvalue <- mc_pvalue_test(y=y, x=x, lstats = lstats, t=t, score = score)
  # pvalue <- pvalue[name.test[istats]]

  #-----------------------------------------------------------------
  # power calculation
  #-----------------------------------------------------------------

  lpwr_d0 <- apply(lstats,1, pwr_d0_stats, k=k, t=t, pquant = 1-alpha)

  larg_min <- t(apply(lstats,1, pwr_arg_min_stats, k=k, t=t, pquant = 1-alpha))
  colnames(larg_min) <- paste0("I", seq(k))

  # for (i in 1 : dim(lstats)[1] ) {
  # #  for (i in 1:length(istats)) {
  #   # lpwr_d0[[i]] <- pwr_d0_stats(stats = lstats[name.test[istats][i], ], k=k, t=t, pquant = 1-alpha)
  #
  #   # split into sublists
  #   lpwr_d0_stats[[i]] <- lpwr_d0[[i]]$pwr_d0
  #   larg_min[[i]] <- lpwr_d0[[i]]$arg_min
  # }

  empty_sublist <- vector(mode = "list", length = length(name.test) + 1 )  # empty sublist with length number of test + 1
  names(empty_sublist) <- c("x", name.test) #[istats])

  empty_list <- rep(list(empty_sublist), length = k)
  names(empty_list) <- paste("I",1:k,sep="")
  plot <- pwr_table <- empty_list

  empty_list <- vector(mode = "list", length = k)  # empty list with length k
  names(empty_list) <- paste("I",1:k,sep="")
  pwr_xtable <- empty_list

  #### cond_item critical region ##########################
  # if (tcond_item) {
  #   sqs <-  lstats["sqs", ]
  #
  #   #  isqs <- which(sqs > quantile(sqs, pquant) )   # index for sqs
  #   # if (missing(alpha_cond)) alpha_cond <- alpha/2
  #
  #   isqs <- which(sqs > quantile(sqs, 1 - alpha_cond) )   # index for sqs
  #
  #   lcond <- list_cond_item(score = score, alpha = alpha)   # list of indices item 1..k with condition
  #   icsqs <- lapply(lcond, function(x) union(x,isqs) )  # list of union sqs and condtion
  # } # end if
  #########################################################

  for (l in seq(k)) { # item l

    for (i in seq(length(name.test))){ # 1:length(istats)) { # test stats i

      lpwr_i <-  plot_pwr_item(item = l,
                               stats = lstats[name.test[i], ],
                               k=k, t=t, pquant = 1-alpha)

      if (i==1) {
        plot[[l]][[1]] <- lpwr_i$pwr_mat$dvec
        pwr_table[[l]][[1]] <- lpwr_i$pwr_table$dvec
      }

      plot[[l]][[i+1]] <- lpwr_i$pwr_mat$pwr
      pwr_table[[l]][[i+1]] <- lpwr_i$pwr_table$pwr

    } # end i

    #### cond_item critical region ##########################
    # if (tcond_item) {
    #
    #   lpwr_i_cond <-  plot_pwr_item(item = l, stats = sqs, icrit = icsqs[[l]], k=k, t=t, pquant=1-alpha)
    #
    #   plot[[l]] <- append(plot[[l]], list(csqs=lpwr_i_cond$pwr_mat$pwr))
    #   pwr_table[[l]] <- append(pwr_table[[l]], list(csqs=lpwr_i_cond$pwr_table$pwr))
    #
    #   # } # end if
    #
    # } # end if tcond_item
    #########################################################


    df <- do.call(cbind.data.frame, pwr_table[[l]])
    rownames(df) <- formatC(df[,1], digits = 2, format = "f")
    df <- round(df, digits= 3)
  #  pwr_xtable[[l]] <- list(xtable = xtable::xtable( df[-1], digits = 3), table_df=df)
    pwr_xtable[[l]] <- df[-1]

  } # end item l

  # descr_table <- descr_calc(lstats = lstats, k=k)
  descr_table <- descr_calc2(lstats = lstats, k=k) # with cv


  results <- (list(lpwr_d0 = lpwr_d0,
                   larg_min = larg_min,
                   pvalue = pvalue,
                   plot = plot,
                   pwr_table = pwr_xtable,
                   # pwr_xtable = pwr_xtable,
                   descr_table = descr_table))

  class(results) <- "MCresults"
  return(results)

} # end powr_calc_t0b

#########################################
#        Power plot functions           #
#########################################
### function plot power item iplot

pwr_plot_item <- function(res, iplot, alpha,
                          legend_title,   # legend_title = "Sample size (n)"
                          tag_title ,
                          name_title ,
                          legend.position,
                          legend.direction,
                          tcolor) {


  # beta_rest <- c(0, res$eta_rest)
  # y <- object$inputmat
  # N <- dim(y)[1]
  # k <- dim(y)[2]

  if (missing(legend.position)) legend.position <- c(.5, .65)


  df <- as.data.frame(res$plot)
  xymelt <- reshape2::melt( data = df,
                            variable.name = "index",
                            id.vars = "x")

  ltype <- c("solid","longdash", "dashed",  "twodash", "dotted")
  lcols <- 0

  # URL: http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
  # The palette with grey:
  cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


  nlabels <- colnames(df)[-1]

  if (class(res) == "crit") {

  nstats <- colnames(df)[2]
  if (nstats == "sqs") {
    l1 <- expression(italic("U")[1]) # sqs
    l2 <- expression( paste(italic("U")[1], " modified"))
    nlabels <- c(l1, l2)
  } else if (nstats == "abs") {
    l1 <- expression(italic("U")[2]) # abs
    l2 <- expression( paste(italic("U")[2], " modified"))
    nlabels <- c(l1, l2)
  } else {
    nlabels <- c(nstats, "mod")
  }

  } # end if class



  if (tcolor) {
    p <-  ggplot2::ggplot(data = xymelt, ggplot2::aes(x = x, y = value, group = index)) +
      ggpubr::theme_pubr() +
      ggplot2::geom_line( ggplot2::aes(linetype = index, col = index)) + #
      ggplot2::xlab(rlang::expr( paste(delta[!!iplot])) )+  # xlab(paste0("delta_",iplot)) +  # xlab("delta") +
      ggplot2::ylab("prob. of rejection") + #ggplot2::ylab("power") +
      # ggplot2::ylim(0,1) +
      ggplot2::labs(title = name_title, tag = tag_title) + # labs(title = name_title, subtitle = name_subtitle) +
      ggplot2::geom_hline(yintercept=alpha, linetype="dotted", color = "black")  +
      # ggplot2::scale_linetype_discrete(name = legend_title) + # scale_linetype_manual(name = legend_title, values = ltype )+
      ggplot2::scale_linetype_discrete(name = legend_title, labels = nlabels) +
      ggplot2::scale_color_discrete(name = legend_title, type = cbPalette)  +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5,size = ggplot2::rel(1.4)), # 15),
            legend.title = ggplot2::element_text(hjust = 1),
            legend.text=ggplot2::element_text(size = ggplot2::rel(1.2)),
            axis.title = ggplot2::element_text( size = ggplot2::rel(1.2)),
            axis.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
            # legend.position = c(.5, .65),  #legend.position = c(.5, .8),
            legend.position = legend.position,
            legend.direction = legend.direction)
            # legend.direction = "vertical") # legend.direction = "horizontal")
  } else {
    p <-  ggplot2::ggplot(data = xymelt, ggplot2::aes(x = x, y = value, group = index)) +
      ggpubr::theme_pubr() +
      ggplot2::geom_line( ggplot2::aes(linetype = index)) + #
      ggplot2::xlab(rlang::expr( paste(delta[!!iplot])) )+  # xlab(paste0("delta_",iplot)) +  # xlab("delta") +
      ggplot2::ylab("prob. of rejection") + # ggplot2::ylab("power") +
      # ylim(0,1) +
      ggplot2::labs(title = name_title, tag = tag_title) + # labs(title = name_title, subtitle = name_subtitle) +
      ggplot2::geom_hline(yintercept=alpha, linetype="dotted", color = "black")  +
      # ggplot2::scale_linetype_discrete(name = legend_title) + # scale_linetype_manual(name = legend_title, values = ltype )+
      ggplot2::scale_linetype_discrete(name = legend_title, labels = nlabels) +
      ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(1.4)), # 15),
            legend.title = ggplot2::element_text(hjust = 1),
            legend.text=ggplot2::element_text(size = ggplot2::rel(1.2)),
            axis.title = ggplot2::element_text( size = ggplot2::rel(1.2)),
            axis.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
            # legend.position = c(.5, .65),  #legend.position = c(.5, .85),
            legend.position = legend.position,
            legend.direction = legend.direction)
            # legend.direction = "vertical") # legend.direction = "horizontal")
  }

  return(p)
}


####
# utility function to select specific case from data.sim.rasch and data.results
#
# i_pers index 1 to 6 : N =  10 15 20 25 30 50
# j_item index 1 to 10: k = 6   8  10  12  14  16  18  20  25  30

load_data <- function(N, k) {

  call<-match.call()

  par <- mc_parameter() # model parameters
  npers <- unlist(par$N_list)
  kitem <- sapply(par$item_list4, length)

  i_pers <- which(npers == N)
  j_item <- which(kitem == k)

  y <- tclboot::data.sim.rasch$inputmat_list[[i_pers]][[j_item]]
  N <- dim(y)[1]
  ng1 <- floor(N/2)  # ng2 <- N - ng1
  splitcr <- c(rep(0, ng1), rep(1, N - ng1)) # binary covariate as used in MC simulation study

  mc_stats <- list( eta_rest =  tclboot::data.sim.rasch$eta_list[[i_pers]][[j_item]][1,], #  estimated item (easiness) parameters!
                    inputmat = y,  #data.sim.rasch$inputmat_list[[i_pers]][[j_item]],
                    score = tclboot::data.sim.rasch$score_list[[i_pers]][[j_item]],
                    lstats = tclboot::data.sim.rasch$stat_list[[i_pers]][[j_item]],
                    t = tclboot::data.sim.rasch$t_list[[i_pers]][[j_item]],
                    splitcr=splitcr)

  lpower <- (list( lpwr_d0 = tclboot::data.results[[i_pers]][[j_item]]$lpwr_d0,
                   larg_min = tclboot::data.results[[i_pers]][[j_item]]$larg_min,
                   pvalue = tclboot::data.results[[i_pers]][[j_item]]$pvalue,
                   plot = tclboot::data.results[[i_pers]][[j_item]]$plot,
                   pwr_table = tclboot::data.results[[i_pers]][[j_item]]$pwr_xtable,
                   descr_table = tclboot::data.results[[i_pers]][[j_item]]$descr_table))

    results <- list(MCobject = mc_stats,
                    result_list = lpower,
                    "call" = call) #,  # plot_list = Filter(length, plot_list))

    class(results) <- "MCboot_case"

    return(results)
}


plotlist <- function(object,
                     nstats,
                     alpha,
                     xlim,
                     legend_title,
                     tag_title,
                     name_title,
                     legend.position,
                     legend.direction,
                     tcolor) {

  # legend_title = "Test" # expression(italic("n"))

  y <- object$MCobject$inputmat
  N <- dim(y)[1]
  k <- dim(y)[2]
  eta_rest <- object$MCobject$eta_rest
  lstats <- as.matrix(object$MCobject$lstats)

  empty_list_k <- vector(mode = "list", length = k)
  names(empty_list_k) <- paste("I",1:k,sep="")
  plotlist <- empty_list_k

  for ( item in seq(k)) {  ###### item

    tab <- object$result_list$plot[[item]] [ c("x",nstats)]
    plot <- as.data.frame (do.call(cbind, tab))
    if (!missing(xlim)) {
      ilow <- which(plot$x == - xlim)
      ihigh <- which(plot$x ==  xlim)
      plot <- plot[ilow:ihigh ,]
    }

    res <- list( eta_rest = eta_rest, y=y,k=k, plot=plot )

    if(item == 1) tag_title <- tag_title # paste0("I",item)
    if(item > 1) tag_title <- ggplot2::element_blank()

    p <- pwr_plot_item(res=res, iplot = item, alpha = alpha,
                       legend_title =legend_title,
                       tag_title = tag_title,
                       name_title = name_title,
                       legend.position=legend.position,
                       legend.direction = legend.direction,
                       tcolor=tcolor)
    plotlist[[item]] <- p

  } # end item ###########

  return(plotlist)

} # end plotlist


##### extension power plots of different sample sizes for given k

plotlist_ext <- function(object_list,
                         nstats,
                         alpha,
                         xlim,
                         legend_title,
                         tag_title,
                         name_title,
                         legend.position,
                         legend.direction,
                         tcolor) {

  # legend_title = "Test" # expression(italic("n"))
  npers <- names(object_list)
  i_pers <- seq(npers)

  empty_list <- vector(mode = "list", length = length(i_pers) + 1)
  names(empty_list) <- c("x", as.character(npers[i_pers])  )

  y <- object_list[[1]]$MCobject$inputmat
  k <- dim(y)[2]

  empty_list_k <- vector(mode = "list", length = k)
  names(empty_list_k) <- paste("I",1:k,sep="")
  plotlist <- empty_list_k

  for ( item in seq(k)) {  ###### item

    tab <- empty_list

    for (ip in i_pers) { ######## N

      if (ip==1) tab[[1]] <- unlist(object_list[[ip]]$result_list$plot[[item]][c("x")])
      tab[[ip+1]] <- unlist(object_list[[ip]]$result_list$plot[[item]][c(nstats)])
    } # end ip


    # tab <- object$result_list$plot[[item]] [ c("x",nstats)]
    plot <- as.data.frame (do.call(cbind, tab))
    if (!missing(xlim)) {
      ilow <- which(plot$x == - xlim)
      ihigh <- which(plot$x ==  xlim)
      plot <- plot[ilow:ihigh ,]
    }

    res <- list( eta_rest = NA, y=NA, k=NA, plot=plot )

    if(item == 1) tag_title <- tag_title # paste0("I",item)
    if(item > 1) tag_title <- ggplot2::element_blank()

    p <- pwr_plot_item(res=res, iplot = item, alpha = alpha,
                       legend_title =legend_title,
                       tag_title = tag_title,
                       name_title = name_title,
                       legend.position=legend.position,
                       legend.direction = legend.direction,
                       tcolor=tcolor)
    plotlist[[item]] <- p

  } # end item ###########

  return(plotlist)

} # end plotlist_ext



#####

#########################################
#  modification  of critical region C   #
#########################################
plotlist_crit <- function(object,
                          d_item, # vector d_item
                          nstats,
                          alpha,
                          xlim,
                          legend_title,
                          tag_title,
                          name_title,
                          legend.position,
                          legend.direction,
                          tcolor,
                          ctr_alpha,
                          ctr_mod) {

  # legend_title = "Test" # expression(italic("n"))
  if (length(nstats) > 1) nstats <- nstats[1] # choose first test if > 1 test
  nstats_c <- "mod"

  y <- object$MCobject$inputmat
  N <- dim(y)[1]
  k <- dim(y)[2]
  eta_rest <- object$MCobject$eta_rest
  lstats <- as.matrix(object$MCobject$lstats)
  t <- as.matrix(object$MCobject$t)

  # sample size in group with covariate = 1
  # if (class(object) == "MCboot_case") {
  #   ng2 <- N -  floor(N/2)
  # } else {
    ng2 <- sum(object$MCobject$splitcr)
  # }
  n_iter <- dim(t)[2]

  empty_list_k <- vector(mode = "list", length = k)
  names(empty_list_k) <- paste("I",1:k,sep="")
  plotlist <- plotlist_t <- plotlist_tot <- resultlist <-  plotlist_t2 <- empty_list_k

  #-----------------------------------------------------------------
  # power calculation
  #-----------------------------------------------------------------
  power_i <- function(d, item, icrit, d_item) {
    d_item[item] <- d  # set ith element of d_item to d
    pwr <- sum(exp(colSums( t[, icrit] * d_item))) / sum(exp( colSums( t * d_item )))
    # pwr <- sum(exp( t[item, icrit] * d )) / sum(exp(t[item,] * d ))
    return(pwr)
  }

  stats <- lstats[nstats,]   #  sqs <-  lstats["sqs", ]
  icrit <- which(stats > stats::quantile(stats, (1 - alpha * ctr_alpha)) )   # index for stats, ctr_alpha = factor
  icrit0 <- which(stats > stats::quantile(stats, (1 - alpha)) ) # crit. region of test statistics without modification

  for ( item in seq(k)) {  ###### item

    tab <- object$result_list$plot[[item]] [ c("x",nstats)]
    dvec <- tab$x

    # if (alpha != .05) tab[[2]] <- sapply(dvec, FUN = power_i, item=item, icrit = icrit0)  # saved result data with alpha = 0.05!
    if((alpha != .05) || (sum( abs(d_item)) > 0)) tab[[2]] <- sapply(dvec, FUN = power_i, item=item, icrit = icrit0, d_item=d_item)

    ### conditionioning critical region ##################################
    imod_low <- which(t[item,] < stats::quantile(t[item,], alpha))
    imod_high <- which(t[item,] > stats::quantile(t[item,], (1 - alpha) ))

    if (ctr_mod == "both") {
        imod <- c(imod_low, imod_high)
    } else if (ctr_mod == "high") {
      imod <- imod_high
    } else if (ctr_mod == "low") {
      imod <- imod_low
    } # end if

    icrit_mod <- union(icrit, imod)
    # tab[[nstats_c]] <- sapply(dvec, FUN = power_i, item=item, icrit = icrit_mod)
    tab[[nstats_c]] <- sapply(dvec, FUN = power_i, item=item, icrit = icrit_mod, d_item = d_item)

    #idiff <- base::setdiff(imod, icrit)
    idiff <- base::setdiff(imod, icrit0)  # difference with and without modification
    resultlist[[item]] <- list(t_item = t[item,],
                              icrit = icrit,
                              icrit0 = icrit0,
                              imod = imod,
                              imod_low = imod_low,
                              imod_high = imod_high,
                              idiff = idiff,
                              icrit_mod = icrit_mod)
    ######################################################################

    plot <- as.data.frame (do.call(cbind, tab))
    if (!missing(xlim)) {
      ilow <- which(plot$x == - xlim)
      ihigh <- which(plot$x ==  xlim)
      plot <- plot[ilow:ihigh ,]
    }

    res <- list( eta_rest= eta_rest, y=y,k=k, plot=plot )
    class(res) <- "crit"

    if(item == 1) tag_title <- tag_title # paste0("I",item)
    if(item > 1) tag_title <- ggplot2::element_blank()

    p <- pwr_plot_item(res=res, iplot = item, alpha = alpha,
                       legend_title =legend_title,
                       tag_title = tag_title,
                       name_title = name_title,
                       legend.position=legend.position,
                       legend.direction = legend.direction,
                       tcolor=tcolor)
    plotlist[[item]] <- p

    #### histogram t, crit. region
    tab2 <-list()

    table_t <- function(x,ng2) {
      lvec <- list()
      for (i in 1 : (ng2-1)) {
        lvec[[i]]  <- table(x)[names(table(x)) == i]
        if (rlang::is_empty(lvec[[i]])) lvec[[i]] <- 0
      }
      vec <- do.call(cbind, lvec)
      colnames(vec) <- seq(ng2-1)
      return(vec)
    }

    tab2$score <- c(1:(ng2-1))
    tab2$t <- table_t(t[item, ], ng2) / n_iter * 100
    tab2$crit <- table_t(t[item, icrit0], ng2) /n_iter  * 100 # rel. freq. without modification
    tab2$mod <- table_t(t[item, icrit_mod ], ng2) /n_iter * 100

    plot2 <- as.data.frame( t(matrix(unlist(tab2), nrow=4, byrow = TRUE)) )
    colnames(plot2) <- c("x","t", nstats, nstats_c)

    p2 <- pwr_hist_item(plot = plot2, iplot = item)

    plotlist_t[[item]] <- p2
    #############################

    # name.labs <- paste0("n = ", N) #expre( paste(italic("n")," = ",N))
    plotlist_tot[[item]] <- ggpubr::ggarrange(p,p2, ncol = 1, nrow = 2, heights = c(20,10), widths = 40)

    #### histogram sufficient statistics
    t_i <- t[item, ]

    tlist <- list()

    tlist[[1]] <-  table(t_i[icrit0])  # crit region t without modification (alpha)
    tlist[[2]] <-  table(t_i[icrit])   # crit region t with modification (alpha * ctr_alpha)
    tlist[[3]] <-  table(t_i[imod])    # modification crit. reg. (alpha)
    # names(tlist) <- c(nstats, paste0(nstats,"_mod"), nstats_c)
    names(tlist) <- c(nstats, paste0(nstats,"_mod"), "alpha_mod")

    df <- mtable(t_i=t_i, tlist = tlist) # data frame of table of T statistcs
    df$mod <- df[,3] + df[,4]

    # lpwr_d0 <- c( sum(df[,2]), sum(df[,3]) , sum(df[,4]), sum(df[,3] + df[,4]) , sum(df[,3] + df[,4] - df[,2])  ) / 100
    lpwr_d0 <- c( sum(df[,2]), sum(df[,3]) , sum(df[,4]), sum(df[,5] ) , sum(df[,5] - df[,2]) ) / 100
    names(lpwr_d0) <- c(nstats, paste0(nstats,"_mod"), paste0("alpha_",nstats_c),  nstats_c, "diff")
    resultlist[[item]]$lpwr_d0 <- lpwr_d0

    # local power with vector d_item at d_i = 0
    lpwr_d_item <- c(tab[[2]][which(tab[[1]] == 0)], tab[[3]][which(tab[[1]] == 0)])
    names(lpwr_d_item) <- c(nstats, nstats_c)
    resultlist[[item]]$lpwr_d_item <- lpwr_d_item

    plotlist_t2[[item]] <- t_hist_item(df = df)

    resultlist[[item]]$t_table <- df


  } # end item ###########

  results <- list( power_item = plotlist,
                   hist_item = plotlist_t,
                   plotlist = plotlist_tot,
                   resultlist = resultlist,
                   hist_crit_reg = plotlist_t2)

  return(results)

} # end plotlist_crit


# histogram sufficient statistic and critical region
# @importFrom grDevices italic

pwr_hist_item <- function(plot, iplot) {
  # requireNamespace(ggplot2)

  name.ylab <- "RF"
  # title <- expression(italic("n"))

  df <-  as.data.frame(plot[,-2])
  ymax <-  5 #

  if(round(max(df[,c(2,3)])) > 5 ) ymax <- round(max(df[,c(2,3)]))
  # if(round(max(dfÇ[,c(2,3)])) < 2.5 ) ymax <- 2.5

  xymelt <- reshape2::melt( data = df,
                            variable.name = "index",
                            id.vars = "x")

  # name.xlab <- rlang::expr( paste(italic(t)[!!iplot]) ) # rlang::expr( paste(italic(t)[!!iplot]) )
  name.xlab <- rlang::expr( paste(t[!!iplot]) )
  x <- value <- index <- delta <- NULL # NULLing out strategy: R CMD check: no visible binding for global variable

  p <- ggplot2::ggplot(data = xymelt, ggplot2::aes(x = x, y = value, fill = index)) +
    ggpubr::theme_pubr() +
    ggplot2::geom_bar( stat="identity", position = ggplot2::position_dodge()) +
    ggplot2::xlab(name.xlab) +
    ggplot2::ylab(name.ylab) +
    ggplot2::scale_fill_grey(start=0.4,end =0.8) +
    ggplot2::scale_x_continuous(breaks= scales::pretty_breaks()) +
    ggplot2::scale_y_continuous(breaks=seq(0.0, ymax, 2.5), limits=c(0, ymax))+
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, size = ggplot2::rel(1.4)),
          axis.title = ggplot2::element_text( size = ggplot2::rel(1.3)),
          axis.text = ggplot2::element_text(size = ggplot2::rel(0.8)),
         # axis.title.x = ggplot2::element_text(face = "italic", size = ggplot2::rel(1.3)),
          legend.title = ggplot2::element_blank(),
          legend.position="none",
          panel.spacing = ggplot2::unit(1, "lines")) +
    ggplot2::facet_grid(. ~ index)

  return(p)

} # end pwr_hist_item


####
mtable <- function(t_i, tlist) {

  ctable <- function(x) {
    m <- matrix(c(as.integer(names(x)), x), ncol = 2, dimnames = list(NULL, c("t", "Freq")))
    return(m)
  }
  n_iter <- length(t_i)

  x <- c(min(t_i):max((t_i)))

  mat <- matrix(0, nrow = length(x), ncol = length(tlist)+1)

  colnames(mat) <- c("t", names(tlist))
  rownames(mat) <- x
  mat[,1] <- x

  for (i in 1: length(tlist)) {

    if ( !(rlang::is_empty(tlist[[i]]))) {

    vec <- ctable(tlist[[i]])

    for (j in 1:length(x)) {

      for (l in 1:dim(vec)[1]) {

        if (vec[l,1] == x[j]) {
          mat[j,i+1] <- vec[l,2] / n_iter *100
        }

      } # end l

    } # end j

    } # end if


  } # end i

  return(as.data.frame(mat))
}


####
t_hist_item <- function(df) {

  # n_iter <- length(t_i)


  # mtable <- function(t_i, tlist) {
  #
  #   ctable <- function(x) {
  #     m <- matrix(c(as.integer(names(x)), x), ncol = 2, dimnames = list(NULL, c("t", "Freq")))
  #     return(m)
  #   }
  #
  #   x <- c(min(t_i):max((t_i)))
  #
  #   mat <- matrix(0, nrow = length(x), ncol = length(tlist)+1)
  #
  #   colnames(mat) <- c("t", names(tlist))
  #   rownames(mat) <- x
  #   mat[,1] <- x
  #
  #   for (i in 1: length(tlist)) {
  #     vec <- ctable(tlist[[i]])
  #
  #     for (j in 1:length(x)) {
  #
  #       for (l in 1:dim(vec)[1]) {
  #
  #         if (vec[l,1] == x[j]) {
  #           mat[j,i+1] <- vec[l,2] / n_iter *100
  #         }
  #
  #       } # end ä
  #
  #     } # end j
  #
  #   } # end i
  #
  #   return(as.data.frame(mat))
  # }

  # df <- mtable(t_i=t_i, tlist = tlist)

  xymelt <- reshape2::melt( data = df,
                            variable.name = "index",
                            id.vars = "t")

  iplot <- 2
  name.xlab <- rlang::expr( paste(t[!!iplot]) )
  name.ylab <- "RF" # "frequency"

  p <- ggplot2::ggplot(data = xymelt, ggplot2::aes(x = t, y = value, fill = index)) +
    ggpubr::theme_pubr() +
    ggplot2::geom_bar( stat="identity", position = ggplot2::position_dodge()) + #geom_bar( position='stack', stat="identity") +
    ggplot2::scale_fill_grey(start=0.8,end =0.5) +
    ggplot2::scale_x_continuous(breaks= df$t) + # scale_x_continuous(breaks= scales::pretty_breaks()) +
    ggplot2::xlab(name.xlab) +
    ggplot2::ylab(name.ylab) +
    ggplot2::theme(legend.title = ggplot2::element_blank(),
         axis.title = ggplot2::element_text( size = ggplot2::rel(1.3)),
         axis.text = ggplot2::element_text(size = ggplot2::rel(0.9)))

  return(p)

}


#### functions from eRm package ####
prettyPaste <- function(...){
  paste(strwrap(paste0(..., collapse = ""), width = getOption("width")), sep="\n", collapse="\n")
}
