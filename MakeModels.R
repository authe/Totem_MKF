# first created: 28 Sep 2023
# last updated: 28 Sep 2023
# author: Andreas Uthemann


MakeModels <- function(y, types, paras, init, ord){
  
  # MakeModels(y, types, paras, init, ord) states space models for type draw types
  # of order "ord" given parameters "paras" in list "models"
  # "init" provides (potentially data-driven) values to set reasonable priors (a0,P0) for the initial state 
  # models[[ind.mod[m]]] : state space model for draw m (vector ind_mod of length M maps draw m to corresponding state space model) 
  # (Return dimension M can be smaller than initial M if type draws are duplicated)

  source("FKFModel2_demeaned_Lagged_pubpriv_heterogenous.R")   # function that creates matrices for state space model
  
  # parameters for state space model
  ord <- ord  # highest order of beliefs
  S <- dim(y)[1] - 1  # number of submitters
  
  # "estimation" parameters
  rho <- paras[1]
  omega <- paras[2]
  sig_u <- paras[3]
  sig_e <- paras[4]
  sig_n <- paras[5]
  sig_z <- paras[6]
  
  
  # data-based parameter values used to set "reasonable" data-driven priors (a0,P0) for the initial state 
  # Important: paras.0 only depend on data y and do not change with type m or "estimation" parameters paras
  # init <- InitialParas(y)
  paras_0 <- init$par  
  rho_0 <- paras_0[1]
  omega_0 <- paras_0[2]
  sig_u_0 <- paras_0[3]
  sig_e_0 <- paras_0[4]
  sig_n_0 <- paras_0[5]
  sig_z_0 <- paras_0[6]
  
  
  ###################################   CREATE STATE SPACE MODELS  ####################################
  
  # create model for W.m weak submitters (unlike for type space (size 2^S), we need at most S different models as the ordering of types does not matter for models since data is rearranged accord to type vector [see below])
  
  Wm <- rowSums(types)  # vector with number of weak submitters for M paths
  size_models <- sort(unique(Wm))  # # figure out distinct values of W (# of weak submitters)
  ind_mod <- sapply(Wm, function(x) match(x, size_models) )  # index of state space for draw m
  
  models <- vector(mode = "list", length = length(size_models))  # list to contain with vectors/matrices of state space models (dt, Tt, HHt, ct, Zt, GGt), i.e. models[[ind.mod[m]]]$Tt
  
  for (i in 1:length(size_models)){
    w <- size_models[i]
    ss <- SP_model_demeaned(rho = rho, omega = omega, sig_u = sig_u, sig_e = sig_e, 
                            sig_n = sig_n, sig_z = sig_z, ord = ord, tol = 1e-15, S = S, W = w)
    
    # initialize prior mean and covariance of state vector alpha_t = (theta_t^(0),..., theta_t^(ord),x_{1,t}^(1),...,x_{W,t}^(ord),e_t, (1-omega) theta_{t-1}^(0)+ omega theta_{t-1}^(1))
    a0 <- c(rep(0, ord+1), rep(0,w * ord), 0, 0)
    
    P0 <- diag(w) %x% matrix(sig_n_0^2, nrow = ord, ncol = ord)
    P0 <- rbind(matrix(0, nrow = (1 + ord), ncol = (w * ord)), P0)
    sig2_theta_0 <- sig_u_0^2 /(1 - rho_0^2)
    aux <- matrix(sig2_theta_0, nrow = (1 + ord), ncol = (1 + ord))
    aux <- rbind(aux, matrix(0, nrow = (w * ord), ncol = (1 + ord)))
    P0 <- cbind(aux, P0)  # covariance matrix for theta^(0:ord) and W submitters' deviations
    P0 <- cbind(P0, rep(0, ord + 1 + (w * ord)))
    P0 <- rbind(P0, c(rep(0, ord + 1 + (w * ord)), sig_e_0)) # add covariance of price error e_t
    aux <- c(rep(rho_0 * sig2_theta_0, (1 + ord)), rep(0, ord * w), 0)
    P0 <- cbind(P0, aux)
    P0 <- rbind(P0, c(aux, sig2_theta_0))  # add covariance of theta_{t-1}^(0)  [approx covariance for (1-omega) theta_{t-1}^(0) + omega theta_{t-1}^(1)]
    
    models[[i]] <- list(dt = ss$dt, Tt = ss$Tt, HHt = ss$HHt, ct = ss$ct,
                        Zt = ss$Zt, GGt = ss$GGt, W = w, a0 = a0, P0 = P0)
  }
  
  mod <- list(models = models, ind_mod = ind_mod)
  
  return(mod)
}