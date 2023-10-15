# first created: 20 Aug 2023
# last updated: 12 Oct 2023
# author: Andreas Uthemann

ll_mcmc <- function(par, y, init, types, ord, crit_eff, l, h, seed = 1){
  # transform parameters to unconstrained problem
  # see mcmc R package pdf for logic of ifelse (avoids numerical problems)
  source("LogLike_MKF.R")
  
  paras <- rep(0, 6)
  #par = c(par[1],0,par[2:5])
  paras[1] <- ifelse(par[1] > 0, 1 / (1 + exp(-par[1])), 
                     exp(par[1]) / (1 + exp(par[1])))
  paras[2] <- ifelse(par[2] > 0, 1 / (1 + exp(-par[2])), 
                     exp(par[2]) / (1 + exp(par[2])))
  #paras[2] <- W/(dim(yt)[1]-1)
  paras[3] <- ifelse(par[3] > 0, (h + l * exp(-par[3])) / (1 + exp(-par[3])),
                     (l + h * exp(par[3])) / (1 + exp(par[3])))
  paras[4] <- ifelse(par[4] > 0, (h + l * exp(-par[4])) / (1 + exp(-par[4])),
                     (l + h * exp(par[4])) / (1 + exp(par[4])))
  paras[5] <- ifelse(par[5] > 0, (h + l * exp(-par[5])) / (1 + exp(-par[5])),
                     (l + h * exp(par[5])) / (1 + exp(par[5])))
  paras[6] <- ifelse(par[6] > 0, (h + l * exp(-par[6])) / (1 + exp(-par[6])),
                     (l + h * exp(par[6])) / (1 + exp(par[6])))
  
  res <- LogLike_MKF(y = y, paras = paras, init = init, types = types, 
                        ord = ord, crit_eff = crit_eff, seed = seed)
  ll_out <- res$LL
  
  # calculate Jacobian for transformed parameters
  # (see https://stats.stackexchange.com/questions/330402/what-is-an-example-of-a-transformation-on-a-posterior-distribution-such-that-the )
  log_jac <- rep(0,6)
  log_jac[1] <- ifelse(par[1] > 0, -log1p(exp(-par[1])) - log1p(exp(par[1])),
                       par[1] - 2 * log1p(exp(par[1])))
  log_jac[2] <- ifelse(par[2] > 0, -log1p(exp(-par[2])) - log1p(exp(par[2])),
                       par[2] - 2 * log1p(exp(par[2])))
  #log.jac[2] <- 0
  log_jac[3] <- ifelse(par[3] > 0,
                       log(h - l) - log1p(exp(-par[3])) - log1p(exp(par[3])),
                       log(h - l) + par[3] - 2 * log1p(exp(par[3])))
  log_jac[4] <- ifelse(par[4] > 0,
                       log(h - l) - log1p(exp(-par[4])) - log1p(exp(par[4])),
                       log(h - l) + par[4] - 2 * log1p(exp(par[4])))
  log_jac[5] <- ifelse(par[5] > 0,
                       log(h - l) - log1p(exp(-par[5])) - log1p(exp(par[5])),
                       log(h - l) + par[5] - 2 * log1p(exp(par[5])))
  log_jac[6] <- ifelse(par[6] > 0,
                       log(h - l) - log1p(exp(-par[6])) - log1p(exp(par[6])),
                       log(h - l) + par[6] - 2 * log1p(exp(par[6])))
  
  ll <- ll_out + sum(log_jac) - 4 * log(h - l)
  
  return(ll)
}