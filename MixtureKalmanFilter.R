# first created: 30 Apr 2023
# last updated: 19 Aug 2023
# author: Andreas Uthemann

MKF <- function(y, types, models, ind_mod, crit_eff = 0.2){
  
  # Given intial draw of submitters types (types: M*S matrix with types[m,s]=1 implying submitter s in draw m is weak type) 
  # and corresponding state space models (models: ind.mod[m] gives the correct model for type draw m
  # MFK returns log-likelihood (LL) and posterior for states (at,Pt) and types (ln.wt) for types having observed data y
  # types are resampled at step t if effective sample size Mt.eff < crit.eff * Mt (i.e. final sample size Mt can be much smaller than initial type draw)
  
  # output: list kf.val with length(kf.val) = number of types that "survived" resampling
  # kf.val[[m]]$LL : log-likelihood of final(!) type draw m given data y
  # kf.val[[m]]$ln.wt : posterior log-weight of draw m given y
  # kf.val[[m]]at : posterior mean of state given y and m
  # kf.val[[m]]Pt : posterior variance of state given y and m
  
  library(matrixStats)  # only used for logExpSum function (can be hand-coded if package does not work, e.g. https://gregorygundersen.com/blog/2020/02/09/log-sum-exp/)
  
  source("KalmanStep.R")
  
  TT <- dim(y)[2]      # number of pricing periods
  S <- (dim(y)[1] - 1) # number of submitters
  
  crit <- crit_eff  # resampling step if effective sample size falls below crit * (sample size)
  
  # only keep unique type draws and weigh resulting draws equally (alternative: keep track of duplicates and reweigh corresponding draws -> computationally intensive for large M)
  ind_mod <- ind_mod[!duplicated(types)]
  types <- unique(types)
  Mt <- dim(types)[1]  # potentially reduced in size later via resampling (hence t in Mt)
  
  kf_val <- vector(mode = "list", length = Mt)  # list to keep track of ln.wt,LL,at,Pt
  
  for(t in 1:TT){
  
    for(m in 1:Mt){  # CAN BE PARALLELIZED
      
      # initial values for path m 
      # (ln.wt: see Sarkka "Bayesian Filtering and Smoothing" p. 122; alternative: W.0*log(omega.0) + (S-W.0)*log(1-omega) and then normalize)
      
      if (t==1){
        kf_val[[m]]$ln_wt <- -log(Mt)
        kf_val[[m]]$LL <- -0.5 * TT * (S + 1) * log(2 * pi)
        kf_val[[m]]$at <- models[[ind_mod[m]]]$a0 
        kf_val[[m]]$Pt <- models[[ind_mod[m]]]$P0
        kf_val[[m]]$types <- types[m,]
      }
      
      # rearrange submitters to fit state-space model
      aux <- y[2:(S + 1), t]
      yt <- c(y[1, t], aux[as.logical(kf_val[[m]]$types)], aux[!as.logical(kf_val[[m]]$types)])
      
      # update state distribution and log-weight for path m given ytt
      ans <- KalmanStep(yt = yt, ln_wt = kf_val[[m]]$ln_wt, at = kf_val[[m]]$at, Pt = kf_val[[m]]$Pt, LL = kf_val[[m]]$LL, 
                        dt = models[[ind_mod[m]]]$dt, Tt = models[[ind_mod[m]]]$Tt, HHt = models[[ind_mod[m]]]$HHt, 
                        ct = models[[ind_mod[m]]]$ct, Zt = models[[ind_mod[m]]]$Zt, GGt = models[[ind_mod[m]]]$GGt)
      
      kf_val[[m]]$ln_wt <- ans$ln_wt
      kf_val[[m]]$LL <- ans$LLt
      kf_val[[m]]$at <- ans$anext
      kf_val[[m]]$Pt <- ans$Pnext
      
    } # end for loop for m in 1:Mt
    
    # normalize weights to sum to 1
    ln_wt_vec <- sapply(1:Mt, function(i) kf_val[[i]]$ln_wt )
    ln_wt_vec <- ln_wt_vec - logSumExp(ln_wt_vec) 
    
    # calculate effective sample size  (1 / sum_m wt^2_m )
    ln_sum_wt2 <- logSumExp(2 * ln_wt_vec) 
    Mt_eff <- 1 / exp(ln_sum_wt2)
    
    # resample if effective sample is too small, else only update log-weights to normalized log-weights
    if (Mt_eff / Mt < crit){
      wt_vec <- exp(ln_wt_vec)
      aux <- sample(1:Mt, size = Mt, replace = TRUE, prob = wt_vec) # CHECK: find R sample fct that allows log-weights
      
      new_paths <- sort(unique(aux))
      Mt <- length(new_paths)
      
      # only keep data for m in new.paths and update log-weights and model indices  (CAREFUL: paths get relabeled) 
      kf_val <- kf_val[new_paths]
      types <- types[new_paths,]
      ind_mod <- ind_mod[new_paths]
      
      new_weights <- sapply(new_paths, function(p) sum(aux == p) / Mt)
      for (m in 1:Mt){
        kf_val[[m]]$ln_wt <- log(new_weights[m])
      }
    } else{
      for (m in 1:Mt){
        kf_val[[m]]$ln_wt <- ln_wt_vec[m]
      }
    }
    
  }  # end for loop for t in 1:TT
  
  return(kf_val)
  
}