SP_model_demeaned <- function(rho, omega, sig_u, sig_e, sig_n, ord, S, W, sig_z=0, tol=1e-15){

  #########################  CREATE STATE SPACE MODEL FOR FKF PACKAGE ########################################
  
  # first created: 18 Jul 2022
  # last updated: 19 AUG 2023
  # author: Andreas Uthemann
  
  # (i) bring model into state space form for FKF package Kalman filter
  # alpha_t = T alpha_{t-1} + H eps_t, eps_t ~ N(0,I) where eps_t = (u_t, e_t, eta_{1,t},...eta{W,t})'
  # y_t = Z alpha_t + G gamma_t , gamma_t ~ N(0,I)
  # Let i=1,...,W be the "weak" submitters with W <= S (number of submitters)
  # Then alpha_t = (theta_t^(0), ..., theta_t^(ord), x_{1,t}^(1),...,x_{1,t}^(ord), ..., x_{W,t}^(ord), e_t, (1-omega) theta_{t-1}^(0) + omega theta_{t-1}^(1))
  # where x_{i,t}^(l) = theta_{i,t}^(l) - theta_t^(l) and 1 <= l <= ord 
  # y_t = (p_t, theta_t^(1) + x_{1,t}^(1) + sig.z z_{1,t},...,theta_t^(1) + x_{W,t}^(1) + sig.z z_{W,t}, theta_t^(0) + sig.z z_{S-W+1,t},...,theta_t^(0) + sig.z z_{S,t} }
  # where p_t = (1-omega) theta_{t-1}^(0) + omega theta_{t-1}^(1) + sig.e e_t
  
  source("createStateSpaceMat2_Lagged_pubpriv_heterogenous.R")
  
  aux <- StateSpaceMatLag(ord = ord, rho = rho, omega = omega, sig_u = sig_u, sig_e = sig_e, sig_n = sig_n, tol = tol)
  M <- aux$M
  N <- aux$N
  Q <- aux$Q
  V <- aux$V
  PP <- aux$PP  # used in initial condition P0 for fkf
  SIG <- V %*% t(V) # used in initial condition P0 for fkf
  
  # transition equation: alpha_t = dt + Tt alpha_{t-1} + Ht eps_t
  # drift
  dt <- matrix(0, nrow = (ord + 1 + (ord * W) + 2))
  
  # transition matrix
  Tt <- diag(W) %x% Q     # %x% is Kronecker product
  Tt <- cbind(Tt, matrix(0, ncol = 2, nrow = (W * ord)))
  Tt <- cbind(matrix(0, ncol = (ord + 1), nrow = (W * ord)), Tt)
  Tt <- rbind(Tt, rep(0, ord + 1 + (ord * W) + 2) )  # row of all zero for e_t
  Tt <- rbind(Tt, c(1 - omega, omega, rep(0, (ord - 1) + (ord * W) + 2)))  #  (1-omega) theta_{t-1}^(0) + omega theta_{t-1}^(1)
  Tt <- rbind(cbind(M, matrix(0, ncol = ((W * ord) + 2), nrow=(ord + 1))), Tt)
  
  # shock variance matrix HHt = E(Ht eta eta' Ht') = Ht Ht'
  Ht <- diag(W) %x% V
  Ht <- cbind(matrix(0, ncol = 2, nrow = (W * ord)), Ht)
  Ht <- rbind(Ht, c(0, 1, rep(0, W)))  # put e_t into current state vector
  Ht <- rbind(Ht, rep(0, 2 + W))  # row of zeros for (1-omega) theta_{t-1}^(0) + omega theta_{t-1}^(1)
  Ht <- rbind(cbind(N, matrix(0, ncol = W, nrow = (ord + 1))), Ht)
  
  HHt <- Ht %*% t(Ht)
  
  # observation equation: y_t = ct + Zt alpha_t + Gt epsilon_t
  # drift
  ct <- matrix(0, nrow = (S + 1))
  
  # observation matrix
  Zt <- matrix(0, nrow = (1 + S), ncol = (ord + 1 + (W * ord) + 2))
  Zt[1, ord + 1 + (ord * W) + 1] <- sig_e  # consensus price is (1-omega) theta_{t-1}^(0) + omega theta_{t-1}^(1) + sig.e e_t
  Zt[1, ord + 1 + (ord * W) + 2] <- 1
  Zt[2:(W + 1), 2] <- 1  # "weak" submitters i=1,...,W submit theta_t^(1) + x_{i,t}^(1) + sig.z z_{i,t}
  for (s in 1:W){
    Zt[(1 + s), (ord + 1 + (s - 1) * ord + 1)] <- 1
  }
  
  if (W < S){
    Zt[(W + 2):(S + 1), 1] <- 1 # "strong" submitters i=W+1,...,S submit theta_t^(0) + sig.z z_{i,t}
  }
  
  # shock variance matrix GGt = E(Gt e e' Gt') = Gt Gt'
  Gt = diag(sig_z, S)
  Gt = rbind(rep(0, S), Gt)
  
  GGt <- Gt %*% t(Gt)
  #GGt <- diag(0,N+1)

  return(list(dt = dt, Tt = Tt, HHt = HHt, ct = ct, Zt = Zt, GGt = GGt, PP = PP, SIG = SIG))
}