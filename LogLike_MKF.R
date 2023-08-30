# first created: 30 Apr 2023
# last updated: 19 Aug 2023
# author: Andreas Uthemann

LogLike_MKF <- function(y, paras, ord, p_class, M){
  
  # computes the log-likelihood of y for model of order "ord" and parameters "paras" using M type draws in the Mixture Kalman Filter 
  
  library(matrixStats)  # only used for logExpSum function (can be hand-coded if package does not work, e.g. https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/)
  
  source("SampleTypes.R")
  source("MixtureKalmanFilter.R")
  
  smp <- SampleTypes(y, paras, ord, p_class, M)
  kf_draws <- MKF(y, smp$types, smp$models, smp$ind_mod)
  
  n_paths <- length(kf_draws)
  
  ln_wT <- sapply(1:n_paths, function(i) kf_draws[[i]]$ln_wt)
  LL_T <- sapply(1:n_paths, function(i) kf_draws[[i]]$LL)
  
  LL <- logSumExp(ln_wT + LL_T)  # LL is the log-likelihood where LL = log( sum_m exp(ln.wT[m]) * exp(LL.T[m]) ) = log( sum_m exp(ln.wT[m] + LL.T[m])
  return(LL)
  
}