ModelMoments <- function(ord, rho, omega, sig_u, sig_e, sig_n, sig_z, tol=1e-15){
  
  # Function to output model implied moments given parameters
  # Returns: 
  # - cross-sectional standard deviation of first order beliefs: stdDev_cross
  # - time-series standard deviation of consensus price: stdDev_price
  #
  # first created: 28 Feb 2024
  # last updated: 29 Feb 2024
  # author: Andreas Uthemann
  
  source("createStateSpaceMat2_Lagged_pubpriv_heterogenous.R")

  library(expm)
  
  # get matrices for state space system
  KalMat <- StateSpaceMatLag(ord=ord, rho=rho, omega=omega, sig_u=sig_u, sig_e=sig_e, sig_n=sig_n, tol=tol)

  # transition equation for weak submitters: theta_t = M_k theta_{t-1} + N_k w_t with w_t ~ N(0, I_2)
  # observation equation: y_t = D_{k,1} theta_t + D_{k,2} theta_{t-1} + R_w w_t + R_n eta_{j,t} with eta_{j,t} = N(0,1)
  # matrices for individual deviations from theta_t
  # x_(j,t) = Q x_(j,t-1) + V eta_(j,t)

  M <- KalMat$M
  N <- KalMat$N
  Q <- KalMat$Q   
  V <- KalMat$V

  # calculate cross-sectional dispersion of beliefs for weak submitters given by (see Nimark 2017, p.15, expression 5.2)
  # Lyapunov equation:
  # Sigma = Q Sigma Q^T + V V^T
  
  i <- 0 
  test_conv <- TRUE
  
  VV <- V %*% t(V)
  Sum_Old <- matrix(0, nrow=ord, ncol=ord)
  count_max <- 1000
  
  while (test_conv & i <= count_max){
    
    Sum_New <- Sum_Old + (Q %^% i) %*% VV %*% (t(Q) %^% i)
    # matrix distance between PP and P
    dist <- max(abs(Sum_New - Sum_Old))
    test_conv <- (dist > tol)
    
    Sum_Old <- Sum_New
    i <- i + 1
  }
  
  SigmaCross_weak <- Sum_New

  # unconditional cross-section variance of 1st order beliefs:
  # weighted sum of weak and strong variance (setting unconditional means for weak and strong to 0)

  sig_cross <- omega * SigmaCross_weak[1,1] + (1-omega) * sig_z^2
  stdDev_cross <- sqrt(sig_cross)

  ###############################  Get unconditional variance of states #################################

  # theta_t = M * theta_{t-1} + N * w_t , where w_t = { 1 , 1 } , , N = {u_t , e_t}
  # Var(theta_t) = Sigma_theta = M * Sigma_theta * M^(t) + N * t(N) 
  # given the  P = A^T * P * A + Q 
  # for us this A = M , Q = N * t(N)

  i <- 0 
  test_conv <- TRUE  

  NN <- N %*% t(N)
  Sum_Old <- matrix(0,nrow=(ord+1), ncol=(ord+1))
  count_max <- 1000
  
  while (test_conv & i <= count_max){ 
    
    Sum_New <- Sum_Old + (M %^% i) %*% NN %*% (t(M) %^% i)
    # matrix distance between PP and P
    dist <- max(abs(Sum_New - Sum_Old))
    test_conv <- (dist > tol)
    
    Sum_Old <- Sum_New
    i <- i + 1
  }
  
  Sigma_theta <- Sum_New

  # consensus price is p_t = omega theta_t^(1,weak) + (1-omega) theta_t^(0,weak) + sig_e epsilon_t
  # = (1-omega, omega, 0, ..., 0) theta_t^(weak) + sig_e epsilon_t = c^T theta_t^(weak) + sig_e epsilon_t
  # -> Var(p_t) = c^T Sigma_theta c + sig_e^2

  c <- matrix(rep(0, ord+1), ncol=1)
  c[1] <- 1 - omega
  c[2] <- omega

  var_price <- t(c) %*% Sigma_theta %*% c + sig_e^2
  stdDev_price <- sqrt(var_price) 


  ################################   Return results   #########################

  res <- list()
  res$"stdDev_cross" <- stdDev_cross
  res$"stdDev_price" <- stdDev_price

  return(res)

}