# first created: 30 Apr 2023
# last updated: 28 Aug 2023
# author: Andreas Uthemann

rm(list=ls())

library(mcmc)
library(tictoc)

source('SimulateData2_demeaned_Lagged_pubpriv_heterogenous.R')
source('LogLike_MKF.R')
source("InitialParasGuess.R") 

# -------------------------- Simulate data for testing ------------------------

seed_sim <- 1
set.seed(seed_sim)

S <- 40  # number of submitters
TT <- 200  # submission dates

ord_sim <- 2

rho_sim <- 0.94
omega_sim <- 0.33 # probability of weak submitter
sig_u_sim <- 0.09
sig_e_sim <- 0.01
sig_n_sim <- 0.09  # 0.09 works fine
sig_z_sim <- 0.01

paras_sim <- c(rho_sim, omega_sim, sig_u_sim, sig_e_sim, sig_n_sim, sig_z_sim)

# draw W weak submitters out of S with i.i.d. prob omega_sim
W <- sum(rbinom(S, 1, omega_sim))   
# have a least one weak submitter (code works with 0)
W <- ifelse(W > 0, W, 1)   

# simulate data (consensus price and S submission time series)
y_sim <- SimulateDataLagged(S = S, W = W, TT = TT, rho = paras_sim[1], 
                            omega = paras_sim[2], sig_u = paras_sim[3], 
                            sig_e = paras_sim[4], sig_n = paras_sim[5], 
                            sig_z = paras_sim[6], ord = ord_sim, seed = seed_sim)

save(y_sim, paras_sim, seed_sim, file = "data/simdata_paraset1.RData")

# calculate log-likelihood at true parameters (mostly test)
M_mkf <- 50000 # types draws for MKF
p_class <- 0.9
tic()
LL <- LogLike_MKF(y = y_sim, paras = paras_sim, ord = ord_sim, 
                  p_class = p_class, M = M_mkf)
toc()

# -----------------------------------------------------------------------------

# -------------------------------- MCMC estimation  ---------------------------

# log-likelihood for MCMC
# (fct ll_mcmc accepts unrestricted parameters; original parameters are scaled)
source("LogLike_MKF_MCMC.R")

# parameters for Mixed Kalman Filter
ord <- 2
p_class <- 0.9
M <- 50000  # number of draws from types space ( 2^S possible type combinations)

# parameters for parameter rescaling 
l <- 0  # lower bound for standard deviation parameters (sig.u, sig.e, sig.n, sig.z)
h <- 5 # upper bound

# initial values for parameter estimation
init <- InitialParas(y_sim)
par_0 <- init$par
par_0 <- c(log(par_0[1:2] / (1- par_0[1:2])), 
           log((par_0[3:6] - l) / (h - par_0[3:6]))) 

# MCMC parameters
scale_mcmc <- diag(6) 
scale_mcmc[1,1] <- 0.2   # rho
scale_mcmc[2,2] <- 0.1   # omega
scale_mcmc[3,3] <- 0.02   # sig_u
scale_mcmc[4,4] <- 0.04  # sig_e
scale_mcmc[5,5] <- 0.025 # sig_n
scale_mcmc[6,6] <- 0.005 # sig_z

nbatch_mcmc <- 5000
burnin_mcmc <- 1000

# ----------------- MCMC estimation and processing of results  ----------------

tic()
out_mcmc <- metrop(ll_mcmc, initial = par_0, nbatch = nbatch_mcmc, 
                   scale = scale_mcmc, y = y_sim,ord = ord, 
                   p_class = p_class,  M = M, l = l, h = h)
toc()

# acceptance rate of chain (ideally between 20% to 30%)
acc_prob <- out_mcmc$accept

# results of chain ( size nbatch * 6)
batch_mcmc <- out_mcmc$batch

# chain of rescaled parameters (plot to check burn in phase and mixing of chain) 
rho_mcmc <- exp(batch_mcmc[,1]) / (1 + exp(batch_mcmc[,1]))
omega_mcmc <- exp(batch_mcmc[,2]) / (1 + exp(batch_mcmc[,2]))
sig_u_mcmc <- (l + h * exp(batch_mcmc[,3])) / (1 + exp(batch_mcmc[,3]))
sig_e_mcmc <- (l + h * exp(batch_mcmc[,4])) / (1 + exp(batch_mcmc[,4]))
sig_n_mcmc <- (l + h * exp(batch_mcmc[,5])) / (1 + exp(batch_mcmc[,5]))
sig_z_mcmc <- (l + h * exp(batch_mcmc[,6])) / (1 + exp(batch_mcmc[,6]))

paras_est <- c(mean(rho_mcmc[burnin_mcmc:nbatch_mcmc]), 
               mean(omega_mcmc[burnin_mcmc:nbatch_mcmc]), 
               mean(sig_u_mcmc[burnin_mcmc:nbatch_mcmc]),
               mean(sig_e_mcmc[burnin_mcmc:nbatch_mcmc]), 
               mean(sig_n_mcmc[burnin_mcmc:nbatch_mcmc]), 
               mean(sig_z_mcmc[burnin_mcmc:nbatch_mcmc]))

paras_est_se <- c(sd(rho_mcmc[burnin_mcmc:nbatch_mcmc]), 
                  sd(omega_mcmc[burnin_mcmc:nbatch_mcmc]), 
                  sd(sig_u_mcmc[burnin_mcmc:nbatch_mcmc]),
                  sd(sig_e_mcmc[burnin_mcmc:nbatch_mcmc]), 
                  sd(sig_n_mcmc[burnin_mcmc:nbatch_mcmc]), 
                  sd(sig_z_mcmc[burnin_mcmc:nbatch_mcmc]))

results_mcmc <- list(out_mcmc, paras_sim, paras_est, paras_est_se,
                     seed_sim, scale_mcmc)
names(results_mcmc) <- c("out_mcmc", "paras_sim", "paras_est", "paras_est_se",
                         "seed", "scale_mcmc")

outfile <- paste0("results_mcmc_simdata_paraset1_scale6_seed", seed_sim, ".RData")
path_out <- "results/"
save(results_mcmc, file = paste0(path_out, outfile))
