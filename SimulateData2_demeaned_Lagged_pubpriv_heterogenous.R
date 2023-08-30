SimulateDataLagged <- function(S, W, TT, rho, omega, sig_u, sig_e, sig_n, sig_z=0, ord, seed=1, tol=1e-15){
  
  # code to simulate submissions for S submitters and belief order ord
  
  # first created: 20 Jul 2022
  # last modified: 20 Aug 2023
  # author: Andreas Uthemann
  
  # import function to create state space system
  library(MASS)
  source("createStateSpaceMat2_Lagged_pubpriv_heterogenous.R")
  
  # simulation parameters
  set.seed(seed)    # fix seed
  # S <- 40         # number of submitters
  # W <- 15         # number of weak submitters
  # TT <- 200       # submission dates
  # rho <- 0.94      # persistence of fundamental theta_t^(0) = rho theta_(t-1)^(0) + u_t
  # omega <- 0.33    # population share of weak submitters
  # sig.u <- 0.09  # sd of fundamental shock u_t
  # sig.e <- 0.01   # sd of noise shock e_t to public signal p_t = (1-omega) theta_{t-1}^(0) + omega theta_{t-1}^(1) + e_t
  # sig.n <- 0.09   # sd of noise shock n_(j,t) to private signal of weak submitters j: s_(j,t) = theta_t^(0) + n_(j,t)
  # sig.z <- 0.01      # sd of measurement error for individual submissions
  # ord <- 2       # order of weak submitters' average belief hierachy theta_t = (theta_t^(0),theta_t^(1),...,theta_t^(ord))
  
  # get matrices for state space system
  SSMat <- StateSpaceMatLag(ord = ord, rho = rho, omega = omega, sig_u = sig_u, sig_e = sig_e, sig_n = sig_n, tol = tol)
  
  #############################  SIMULATION  ##########################################################
  
  ######## simulate time series for mean beliefs
  
  # draw time series of aggregate shocks w
  w <- mvrnorm(TT+1, rep(0,2), diag(2))
  # draw initial condition using prior N(0,PP)
  theta0 <- mvrnorm(1, rep(0, ord + 1), SSMat$PP)
  theta <- theta0
  
  aux <- theta0
  for (t in 2:(TT + 1)){
    aux <- SSMat$M %*% aux + SSMat$N %*% w[t,]
    theta <- rbind(theta, t(aux))
  }
  
  
  ###### simulate time series for individual submissions
  
  # generate consensus price p_t = (1-omega) theta^(0)_{t-1} + omega theta^(1)_{t-1} + sig.e e_t
  price = (1 - omega) * theta[(1:TT), 1] + omega * theta[(1:TT), 2] + sig_e * w[2:(TT + 1), 2]
  price = c(0, price)
  
  # generate beliefs for W "weak" submitters
  
  beliefs <- list()
  
  for (s in 1:W){
    
    y0 <- mvrnorm(1, rep(0, ord + 1), SSMat$PP)
    y0 <- y0[2:(ord + 1)]  # initial individual beliefs
    y <- y0
    
    aux <- y0
    for (t in 2:(TT + 1)){
      
      priv_signal <- theta[t, 1] + sig_n * rnorm(1)
      signals = matrix( c(priv_signal, price[t]), nrow = 2)
      
      aux <- SSMat$M_ind %*% aux + SSMat$KK %*% ( signals - SSMat$D1 %*% SSMat$M_ind %*% aux - SSMat$D2 %*% aux) 
      y <- rbind(y, t(aux))
    }
    
    beliefs[[s]] <- y[2:(TT+1),]
  }
  
  # create submission data
  submissions <- c()
  
  # each "weak" submitter j=1,...,W submits his best estimate theta_(j,t)^(1) plus measurement error for t=1,...,TT 
  for (s in 1:W){
    submissions <- cbind(submissions, beliefs[[s]][, 1] + sig_z * rnorm(TT))
  }
  
  # each "strong" submitter j=1,...,W submits the fundamental theta_(j,t)^(0) plus measurement error
  for (s in (W + 1):S){
    submissions <- cbind(submissions, theta[2:(TT + 1), 1] + sig_z * rnorm(TT))
  }
  
  # first row of sim.data is period's consensus price, row i > 1 corresponds to submitter i's submissions for t=1,...,TT
  sim_data <- t(cbind(price[2:(TT + 1)], submissions))
  
  return(sim_data)
}