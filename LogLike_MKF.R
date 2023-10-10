# first created: 30 Apr 2023
# last updated: 28 Sep 2023
# author: Andreas Uthemann

LogLike_MKF <- function(y, paras, init, types, ord, crit_eff){
  
  # computes the log-likelihood of y for model of order "ord" and parameters "paras" and "init" using matrix types in the Mixture Kalman Filter 
  # crit_eff sets critical value for resampling (as percentage of total sample - see MixtureKalmanFilter.R)

  library(matrixStats)  # only used for logExpSum function (can be hand-coded if package does not work, e.g. https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/)
  
  source("MakeModels.R")
  source("MixtureKalmanFilter.R")
  
  mod <- MakeModels(y, types, paras, init, ord)
  kf_draws <- MKF(y, types, mod$models, mod$ind_mod, crit_eff = crit_eff)
  
  n_paths <- length(kf_draws)
  
  ln_wT <- sapply(1:n_paths, function(i) kf_draws[[i]]$ln_wt)
  LL_T <- sapply(1:n_paths, function(i) kf_draws[[i]]$LL)
  
  LL <- logSumExp(ln_wT + LL_T)  # LL is the log-likelihood where LL = log( sum_m exp(ln.wT[m]) * exp(LL.T[m]) ) = log( sum_m exp(ln.wT[m] + LL.T[m])
  
  return(LL)
  
}