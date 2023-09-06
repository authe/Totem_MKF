# first created: 30 Apr 2023
# last updated: 6 Sep 2023
# author: Andreas Uthemann

# function calculates Kalman updating step for a given prior for state vector (at, Pt) and log sample weight for a given sample path (ln.wt) given new observation (yt)
# output: updated log-likelihood (LLt), log-weight of sample path (ln.wt), posterior of state for t (att, Ptt) and prior for t+1 (anext, Pnext)

KalmanStep <- function(yt, ln_wt, at, Pt, LL, dt, Tt, HHt, ct, Zt, GGt){
  
  # at and Pt are prior mean and covariance matrix for state vector at t: E(alpha_t | y_1,...y_(t-1)), Var(alpha_t | y_1,...y_(t-1))
  # LL is the log-likelihood of data (y_1,...,y_{t-1})
  # ln.wt is the weight of the given sample path using data (y_1,...,y_{t-1})
  
  # Adjusting for missing observation in yt:
  # modify observation equation yt = ct + Zt alpha_t + Ht eta_t with eta_t ~ N(0,I_m) for missing observations in yt
  # WWt is a selection matrix of dim(yt). Starting from an identity matrix, if yt(j) is missing, row j of WWt is deleted. 
  # New observation equation is yt* = Zt* alpha_t + Gt* eta_t with yt* = WWt yt, Zt* = WWt Zt and Gt* = WWt Gt (implying GGt* = WWt GGt WWt')
  
  #library(mvtnorm)  # only used for density of multivariate normal when updating weights
  
  if (anyNA(yt)){
    
    miss <- is.na(yt)
    
    WWt <- diag(length(yt))
    WWt <- WWt[which(!miss),]
    
    yt <- yt[!miss]
    ct <- ct[!miss]
    
    Zt <- WWt %*% Zt
    GGt <- WWt %*% GGt %*% t(WWt)
  }
  
  ythat <- ct + Zt %*% at        # predicted yt
  vt <- yt - ythat               # forecast error
  tmp_Pt_trZt <- Pt %*% t(Zt)    # avoids repeat calculation of Pt * Zt'    # crossprod(Pt, t(Zt))
  Ft <- Zt %*% tmp_Pt_trZt + GGt # covariance matrix of forecast errors vt
  
  # calculate inverse of Ft (via Cholesky decomposition: Ft is symmetric & positive definite)
  # faster than: invFt <- solve(Ft) 
  cholFt <- chol(Ft)
  invFt <- chol2inv(cholFt)
  
  Kt <- tmp_Pt_trZt %*% invFt   # Kalman gain
  att <- at + Kt %*% vt     # updated mean
  Ptt <- Pt - tmp_Pt_trZt %*% t(Kt)    # updated variance
  
  # update log-likelihood: LL(t) = LL(t-1) - 0.5 log det(Ft) - 0.5 * v' * Ft^(-1) * v
  
  # calculate determinant of Ft using its Cholesky decomposition : det(M) = det(L L') = det(L)det(L') = det(L)^2 and det(L) = prod(diag(L))
  # faster than: logDet <- log(det(Ft))
  logDet <- 2 * sum(log(diag(cholFt)))
  MahaDist <- t(vt) %*% invFt %*% vt
  
  aux <- logDet + MahaDist
  LLt <- LL - 0.5 * aux
  
  # update log weight of sample path
  #ln.ut <- dmvnorm(c(vt), mean = rep(0,length(c(vt))), sigma = Ft, log = TRUE, checkSymmetry = FALSE)  # log density of forecast error
  ln_ut <- -0.5 * aux
  ln_wt <- ln_ut + ln_wt  
  
  # (t+1) prediction step
  anext <- dt + Tt %*% att
  Pnext <- Tt %*% Ptt %*% t(Tt) + HHt
  
  out_data <- list(LLt = LLt, ln_wt = ln_wt, att = att, Ptt = Ptt, anext = anext, Pnext = Pnext)
  return(out_data)
}