# first created: 30 Apr 2023
# last updated: 28 Sep 2023
# author: Andreas Uthemann

rm(list=ls())

library(mcmc)
library(tictoc)

source('SimulateData2_demeaned_Lagged_pubpriv_heterogenous.R')
source('LogLike_MKF.R')
source("InitialParasGuess.R")
source("DrawTypes.R")

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

# introduce missing values (optional)
p_miss <- 0.01   # percent of missing observation per submitter
miss <- rbinom(S*TT, 1, p_miss)
miss <- matrix(as.logical(miss), nrow = S)

for (s in 1:40){
  y_sim[(s+1), miss[s, ]] <- NA
}

#save(y_sim, paras_sim, seed_sim, file = "data/simdata_paraset1_NAs.RData")

# data-driven guesses for parameters and submitter types
init <- InitialParas(y_sim)

# calculate log-likelihood at true parameters (mostly test)

M_mkf <- 50000 # types draws for MKF
p_class <- 0.9
crit_eff = 0.2
types_mkf <- DrawTypes(M_mkf, init, p_class)

tic()
LL <- LogLike_MKF(y = y_sim, paras = paras_sim, init = init, 
                  types = types_mkf, ord = ord_sim, crit_eff = crit_eff)
toc()

# -----------------------------------------------------------------------------

# -------------------------------- MCMC estimation  ---------------------------

# log-likelihood for MCMC
# (fct ll_mcmc accepts unrestricted parameters; original parameters are scaled)
source("LogLike_MKF_MCMC.R")

# parameters for Mixed Kalman Filter
ord <- 2
p_class <- 0.9
crit_eff = 0.2
M <- 50000  # number of draws from types space ( 2^S possible type combinations)

# parameters for parameter rescaling 
l <- 0  # lower bound for standard deviation parameters (sig.u, sig.e, sig.n, sig.z)
h <- 5 # upper bound

# draw M types for MKF
types <- DrawTypes(M_mkf, init, p_class)

# initial values for parameter estimation
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

nbatch_mcmc <- 10
burnin_mcmc <- 1

# ----------------- MCMC estimation and processing of results  ----------------

tic()
out_mcmc <- metrop(ll_mcmc, initial = par_0, nbatch = nbatch_mcmc, 
                   scale = scale_mcmc, y = y_sim, init = init, types = types, 
                   ord = ord, crit_eff = crit_eff, l = l, h = h)
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

outfile <- paste0("results_mcmc_simdata_NAs_paraset1_scale6_seed", seed_sim, ".RData")
path_out <- "results/"
save(results_mcmc, file = paste0(path_out, outfile))