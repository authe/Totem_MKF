# first created: 17 Apr 2024
# last updated: 24 Apr 2024
# author: Andreas Uthemann

InitialParas <- function(y, paras0 = rep(NA, 6), seed=1){

  # IntialParas(y) returns an initial guess for model parameters based on moments of submission data y
  # (can be used to pick initial value for optimization routine and to set "reasonable" priors (a0,P0) for the initial state)
  
  library(lmtest) # needed for DW test for autocorrelation
  set.seed(seed)

  S <- dim(y)[1] - 1  # number of submitters
  TT <- dim(y)[2]     # number of submission periods

  # at the moment these parameters have to be provided
  rho_hat <- paras0[1]
  sig_u_hat <- paras0[3]
  sig_n_hat <- paras0[5]
  
  # sig_e is the stddev of the demeaned consensus price
  eps <- t(y[1, ])
  if (is.na(paras0[4])) {
    sig_e_hat <- sd(eps)
  } else {
    sig_e_hat <- paras0[4]
  }

  eps <- eps / sd(eps) # e_t in model is normalized 
  
  # starting guess for weak v. strong submitters
  # demeaned submission load with opposite sign on (theta_t^1 - theta_t^0)
  # after removing this component (using PC1), strong submission should be z_{i,t}
  # and have 0 autocorrelation where weak submissions are x_{i,t} and correlated over time
  
  # extract PC1 from submissions
  subs <- t(y[2:(S + 1), ])
  pc <- prcomp(subs, center = TRUE, scale. = TRUE)
  pc1 <- pc$x[, 1]
  
  # identify two groups based on sign of loading on pc1
  # residual of strong submitters has lower autocorrelation than for weak
  
  # fct to get loading on pc1 and DW test for residuals for each submission (x) 
  pcreg <- function(x){
    reg_pc1 <- lm(x ~ pc1)
    b1 <- reg_pc1$coefficients[2]
    sd_res <- sd(reg_pc1$residuals)
    aux <- dwtest(x ~ pc1)
    dw_test <- aux$p.value
    return(c(b1, dw_test, sd_res))
  }
  
  group_test <- apply(subs, 2, pcreg)
  group1 <- group_test[1, ] > 0
  dw1 <- mean(group_test[2, group1])
  dw2 <- mean(group_test[2, !group1])
  
  # weak have smaller p-value for DW test (residuals are autocorrelated)
  weak_hat <- group1
  if (dw2 <= dw1) {
    weak_hat <- !weak_hat
  }
  
  # omega can be obtained from group mean submission
  mean_weak <- rowMeans(subs[, weak_hat])
  mean_strong <- rowMeans(subs[, !weak_hat])
  
  if (is.na(paras0[2])) {
    omega_hat <- mean(mean_strong / (mean_strong - mean_weak))
  } else {
    omega_hat <- paras0[2]
  }
  
  # sig_z is the sd of strongs' residuals
  if (is.na(paras0[6])) {
    sig_z_hat <- mean(group_test[3, !weak_hat])
  } else {
    sig_z_hat <- paras0[6]
  }
  
  paras_0 <- c(rho_hat, omega_hat, sig_u_hat, sig_e_hat, sig_n_hat, sig_z_hat)
  names(paras_0) <- c("rho0", "omega0", "sig_u0", "sig_e0", "sig_n0", "sig_z0")
  ret <- list(paras_0, weak_hat)
  names(ret) <- c("par", "weak")
  
  return(ret)
}