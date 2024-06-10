SP_model_demeaned <- function(rho, omega, sig_u, sig_e, sig_n, ord, S, W, sig_z=0, tol=1e-15){

  #########################  CREATE STATE SPACE MODEL FOR FKF PACKAGE ########################################
  
  # first created: 10 Apr 2024
  # last updated: 10 Jun 2024
  # author: Andreas Uthemann
  
  # (i) bring model into state space form for Kalman Filter
  # alpha_t = T alpha_{t-1} + H eps_t, eps_t ~ N(0,I) where eps_t = (u_t, e_t, eta_{1,t},...eta{W,t})'
  # y_t = Z alpha_t + G gamma_t , gamma_t ~ N(0,I)
  # Let i=1,...,W be the "weak" submitters with W <= S (number of submitters)
  # Then alpha_t = (theta_t^(0), ..., theta_t^(ord), x_{1,t}^(1),...,x_{1,t}^(ord), ..., x_{W,t}^(ord), e_t)
  # where x_{i,t}^(l) = theta_{i,t}^(l) - theta_t^(l) and 1 <= l <= ord.
  # (ii) raw submission data is demeaned with cross-sectional mean
  # weak submitters: s_{i,t} = a_t + theta_t^(1) + x_{i,t},
  # strong submitters: s_{i,t} = a_t + theta_t^(0) + sig_z z_{i,t},
  # cross-section mean m_t = a_t + (1-omega) theta_t^(0) + omega theta_t^(1)
  # consensus price (lagged) p_t = m_{t-1} + sig_e e_t
  # demeaned data:
  # strong: y_{i,t} = submission_{i,t} - m_t = omega (theta_t^(0) - theta_t^(1)) + sig_z z_{i,t},
  # weak: y_{i,t} = submission_{i,t} - m_t = x_{i,t} - (1-omega) (theta_t^(0) - theta_t^(1)),
  # consensus price (lagged): q_t = p_t - m_{t-1} = sig_e e_t
  # shape of data for estimation y_t = (consensus price (lagged), weak submitters, strong submitters):
  # y_t = (q_t, y_{1,t}, ...., y_{W,t}, y_{W+1,t}, ...., y_{S,t}}.

  source("createStateSpaceMat2_Lagged_pubpriv_heterogenous.R")

  aux <- StateSpaceMatLag(ord = ord, rho = rho, omega = omega, sig_u = sig_u, sig_e = sig_e, sig_n = sig_n, tol = tol)
  M <- aux$M
  N <- aux$N
  Q <- aux$Q
  V <- aux$V
  PP <- aux$PP  # used in initial condition P0 for Kalman filter
  SIG <- V %*% t(V) # used in initial condition P0 for Kalman filter

  # transition equation: alpha_t = dt + Tt alpha_{t-1} + Ht eps_t
  # drift
  dt <- matrix(0, nrow = (ord + 1 + (ord * W) + 1))

  # transition matrix
  Tt <- diag(W) %x% Q     # %x% is Kronecker product
  Tt <- cbind(Tt, matrix(0, ncol = 1, nrow = (W * ord)))
  Tt <- cbind(matrix(0, ncol = (ord + 1), nrow = (W * ord)), Tt)
  Tt <- rbind(Tt, rep(0, ord + 1 + (ord * W) + 1))  # row of all zero for e_t
  Tt <- rbind(cbind(M, matrix(0, ncol = ((W * ord) + 1), nrow = (ord + 1))), Tt)

  # shock variance matrix HHt = E(Ht eta eta' Ht') = Ht Ht'
  Ht <- diag(W) %x% V
  Ht <- cbind(matrix(0, ncol = 2, nrow = (W * ord)), Ht)
  Ht <- rbind(Ht, c(0, 1, rep(0, W)))  # put e_t into current state vector
  Ht <- rbind(cbind(N, matrix(0, ncol = W, nrow = (ord + 1))), Ht)

  HHt <- Ht %*% t(Ht)

  # observation equation: y_t = ct + Zt alpha_t + Gt epsilon_t
  # drift
  ct <- matrix(0, nrow = (S + 1))

  # observation matrix
  Zt <- matrix(0, nrow = (1 + S), ncol = (ord + 1 + (W * ord) + 1))

  # demeaned consensus price q_t = sig_e e_t
  Zt[1, ord + 1 + (ord * W) + 1] <- sig_e

  # weak submitters i=1,...,W submit x_{i,t}^(1) - (1-omega)(theta_t^(0)-theta_t^(1))
  Zt[2:(W + 1), 1] <- -(1 - omega)
  Zt[2:(W + 1), 2] <- 1 - omega
  for (s in 1:W){
    Zt[(1 + s), (ord + 1 + (s - 1) * ord + 1)] <- 1
  }
  # "strong" submitters i=W+1,...,S submit omega (theta_t^(0) - theta_t^(1)) + sig_z z_{i,t}
  if (W < S){
    Zt[(W + 2):(S + 1), 1] <- omega
    Zt[(W + 2):(S + 1), 2] <- -omega
  }

  # shock variance matrix GGt = E(Gt e e' Gt') = Gt Gt'
  Gt <- diag(c(rep(0, W), rep(sig_z, S - W)), S)
  Gt <- rbind(rep(0, S), Gt)

  GGt <- Gt %*% t(Gt)

  return(list(dt = dt, Tt = Tt, HHt = HHt, ct = ct, Zt = Zt, GGt = GGt, PP = PP, SIG = SIG))
}