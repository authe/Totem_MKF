# first created: 30 Apr 2023
# last updated: 19 Aug 2023
# author: Andreas Uthemann

InitialParas <- function(y){

  # IntialParas(y) returns an initial guess for model parameters based on moments of submission data y
  # (can be used to pick initial value for optimization routine and to set "reasonable" priors (a0,P0) for the initial state)
  
  S <- dim(y)[1] - 1  # number of submitters
  TT <- dim(y)[2]     # number of submission periods
  
  # starting guess for weak v. strong submitters based on ts stddev of their submissions (strong > weak because of stronger updating on info)
  sd_quotes <- apply(y[2:(S+1),], MARGIN = 1, FUN = sd)
  groups <- kmeans(sd_quotes, centers = 2, iter.max = 10)
  # weak submitters are group with lower ts stddev of submissions
  if (groups$centers[1]  <= groups$centers[2]){
    weak <- groups$cluster == 1
  } else{
    weak <- groups$cluster == 2
  }
  W_0 <- sum(weak == TRUE)
  
  # reorder yt.sim so that weak submitters appear first, then strong submitters (groups$cluster keeps track of original order)
  aux <- y[2:(S+1),]
  y_weak <- aux[weak,]
  y_strong <- aux[!weak,]
  y_reorder <- rbind(y[1,], y_weak, y_strong)
  
  meanIV <- colMeans(y_reorder[(W_0 + 2):(S + 1),])
  est_AR1 <- arima(meanIV, c(1, 0, 0))
  
  rho_0 <- min(abs(est_AR1$coef[1]), 0.99)
  sig_u_0 <- sqrt(est_AR1$sigma2)
  
  aux <- sweep(y_reorder[2:(W_0 + 1), 1:(TT - 1)], 2, y_reorder[1, 2:TT])
  sig_n_0 <- mean(apply(aux, 2, sd, na.rm=TRUE))
  
  aux <- sweep(y_reorder[(W_0 + 2):(S + 1),], 2, meanIV)
  sig_z_0 <- mean(apply(aux, 2, sd, na.rm=TRUE))
  
  aux <- colMeans(y_reorder[2:(S + 1), 1:(TT - 1)]) - y_reorder[1, 2:(TT)] # e_t = (1-omega) theta^(0)_{t-1} + omega theta^(1)_{t-1} - p_{t-1}   
  sig_e_0 <- sd(aux)
  
  omega_0 <- W_0 / S
  
  
  paras_0 <- c(rho_0, omega_0, sig_u_0, sig_e_0, sig_n_0, sig_z_0)
  ret <- list(paras_0, weak)
  names(ret) <- c("par", "weak")
  
  return(ret)
}