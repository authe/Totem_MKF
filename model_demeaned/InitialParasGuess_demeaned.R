# first created: 17 Apr 2024
# last updated: 18 Apr 2024
# author: Andreas Uthemann

InitialParas <- function(y, seed=1){

  # IntialParas(y) returns an initial guess for model parameters based on moments of submission data y
  # (can be used to pick initial value for optimization routine and to set "reasonable" priors (a0,P0) for the initial state)
  
  library(lmtest) # needed for DW test for autocorrelation
  set.seed(seed)

  S <- dim(y)[1] - 1  # number of submitters
  TT <- dim(y)[2]     # number of submission periods
  
  # sig_e is the stddev of the demeaned consensus price
  eps <- t(y[1, ])
  sig_e_hat <- sd(eps)
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
  
  omega_hat <- mean(mean_strong / (mean_strong - mean_weak))
  
  # sig_z is the sd of strongs' residuals
  sig_z_hat <- mean(group_test[3, !weak_hat])
  
  # extract theta_t^(1) - theta_t from subs
  aux <- - mean_strong / omega_hat
  aux2 <- mean_weak / (1-omega_hat)
  dtheta <- rowMeans( cbind(aux, aux2))
  
  # for 1st order model we have
  # dtheta_t = (1-k) rho dtheta_{t-1} + k12 eps_{t-1} + k11 u_t
  reg2 <- lm(dtheta[2:TT] ~ dtheta[1:(TT-1)] + eps[2:TT])
  
  reg2_coeff <- reg2$coefficients
  reg2_res <- reg2$residuals
  
  k_rho <- reg2_coeff[2]
  k12 <- reg2_coeff[3]
  k11 <- sd(reg2_res)
  
  # we can solve for rho knowing that k12 = rho (k - k11)
  # this implies: rho = [ (1-k) rho + k12 ] / (1 - k11)
  
  rho_hat <- (k_rho + k12) / (1 - k11)
  if (rho_hat >= 1){
    rho_hat <- 0.99
  } else if (rho_hat <= 0){
    rho_hat <- 0.01
  }
  
  # find xi from equation for k
  k_hat <- (rho_hat - k_rho) / rho_hat
  eq_xi <- function(xi){
    val <- 0.5 + (1/(2*rho_hat^2))*(sqrt((1-rho_hat)^2 + xi)*sqrt((1+rho_hat)^2 + xi) - (1+xi)) - k_hat
    return(val)
  }
  
  aux <- uniroot(eq_xi, c(0, 10))
  xi_hat <- aux$root
  
  # find psi from equation for k11
  eq_psi <- function(psi){
    val <- (xi_hat*psi + rho_hat^2 * k_hat)/(xi_hat*psi + rho_hat^2 + psi/(1 + psi)) - k11
    return(val)
  }
  
  aux <- uniroot(eq_psi, c(0, 1))
  psi_hat <- aux$root
  
  sig_n_hat <- sqrt(psi_hat/(1-psi_hat)) * sig_e_hat
  
  zeta_hat <- xi_hat * psi_hat
  sig_u_hat <- sqrt(zeta_hat) * sig_e_hat
  
  
  paras_0 <- c(rho_hat, omega_hat, sig_u_hat, sig_e_hat, sig_n_hat, sig_z_hat)
  names(paras_0) <- c("rho0", "omega0", "sig_u0", "sig_e0", "sig_n0", "sig_z0")
  ret <- list(paras_0, weak_hat)
  names(ret) <- c("par", "weak")
  
  return(ret)
}