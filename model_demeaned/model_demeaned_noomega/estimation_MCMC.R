# first created: 20 Aug 2023
# last updated: 12 Aug 2024
# author: Andreas Uthemann

# main program for MCMC estimation of MKF
rm(list = ls())

library(mcmc)
library(tictoc)

seed <- 1
set.seed(seed)

# ------- contract specification and input & output file paths ----------------

ind <- "SPX"
term <- "12"
strike <- "150"
opttyp <- "CE"

# path to input data and output/results
path_data <- "data/"
data_file <- paste0("estdata_", ind, "Term", term, "Strike", strike, 
                    opttyp, ".RData")

path_out <- "results/"
outfile <- paste0("results_mcmc_", ind, "Term", term, "Strike", strike, opttyp,
                   "_Seed", seed, ".RData")
# output file for posterior weights wT of MCMC
path_types <- paste0(path_out, "types_mcmc_", ind, "Term", term, "Strike",
                     strike, opttyp, "_Seed", seed, ".RData")


# -------------------- load data for estimation -------------------------------                    

load(file = paste0(path_data, data_file))

y_est <- log(ydat)  # panel with lagged consensus IV and submissions

# -------------------- set up log-likelihood (MKF)  ---------------------------

# log-likelihood for MCMC
# (fct ll_mcmc accepts unrestricted parameters; original parameters are scaled)
source("LogLike_MKF_MCMC.R")
source("DrawTypes.R")
source("InitialParasGuess.R")

# parameters for Mixed Kalman Filter
ord <- 2
p_class <- 0.9
crit_eff = 0.2
M <- 50000  # number of draws from types space ( 2^S possible type combinations)

LL_spec <- c(ord, p_class, crit_eff, M)
names(LL_spec) <- c("ord", "p_class", "crit_eff", "M")

# draw types (M*S matrix)
init <- InitialParas(y_est, seed = seed)  # init$weak is initial guess for weak/strong
types <- DrawTypes(M, init, p_class)

# -------------------- set parameters for MCMC estimation ---------------------

# list to keep track of index of surviving types and corresponding weights
types_mcmc <- list()
save(types_mcmc, file = path_types)

# bounds for parameter rescaling 
l <- 0  # lower bound for stddev parameters (sig.u, sig.e, sig.n, sig.z)
h <- 1  # upper bound

# initial values for parameter estimation
par_0 <- init$par
par_0 <- c(log(par_0[1] / (1 - par_0[1])), 
           log((par_0[3:6] - l) / (h - par_0[3:6])))

# MCMC parameters
scale_mcmc <- diag(5) * 0.01   # (h=5) 0.1 -> 0.066, 0.01 -> 0.1, 0.005 -> 0.066, 1 -> 0
scale_mcmc[1,1] <- 0.2   # rho
#scale_mcmc[2,2] <- 0.02   # sig_u
#scale_mcmc[3,3] <- 0.0005  # sig_e 0.0001 bad
#scale_mcmc[4,4] <- 0.002 # sig_n
#scale_mcmc[5,5] <- 0.002 # sig_z

nbatch_mcmc <- 30
burnin_mcmc <- 10

# -------------------------------- MCMC estimation  ---------------------------

tic()
out_mcmc <- metrop(ll_mcmc, initial = par_0, nbatch = nbatch_mcmc, 
                   scale = scale_mcmc, y = y_est, types = types, 
                   init = init, ord = ord, crit_eff = crit_eff, l = l, h = h,
                   path_types = path_types)
toc()

# ---------------------------- process results  -------------------------------

# acceptance rate of chain (ideally between 20% to 30%)
acc_prob <- out_mcmc$accept

# results of chain ( size nbatch * 6)
batch_mcmc <- out_mcmc$batch

# chain of rescaled parameters (plot to check burn in phase and mixing of chain) 
rho_mcmc <- exp(batch_mcmc[, 1]) / (1 + exp(batch_mcmc[, 1]))
sig_u_mcmc <- (l + h * exp(batch_mcmc[, 2])) / (1 + exp(batch_mcmc[, 2]))
sig_e_mcmc <- (l + h * exp(batch_mcmc[, 3])) / (1 + exp(batch_mcmc[, 3]))
sig_n_mcmc <- (l + h * exp(batch_mcmc[, 4])) / (1 + exp(batch_mcmc[, 4]))
sig_z_mcmc <- (l + h * exp(batch_mcmc[, 5])) / (1 + exp(batch_mcmc[, 5]))

paras_est <- c(mean(rho_mcmc[burnin_mcmc:nbatch_mcmc]),
               mean(sig_u_mcmc[burnin_mcmc:nbatch_mcmc]),
               mean(sig_e_mcmc[burnin_mcmc:nbatch_mcmc]),
               mean(sig_n_mcmc[burnin_mcmc:nbatch_mcmc]),
               mean(sig_z_mcmc[burnin_mcmc:nbatch_mcmc]))

paras_est_se <- c(sd(rho_mcmc[burnin_mcmc:nbatch_mcmc]),
                  sd(sig_u_mcmc[burnin_mcmc:nbatch_mcmc]),
                  sd(sig_e_mcmc[burnin_mcmc:nbatch_mcmc]),
                  sd(sig_n_mcmc[burnin_mcmc:nbatch_mcmc]),
                  sd(sig_z_mcmc[burnin_mcmc:nbatch_mcmc]))

results_mcmc <- list(out_mcmc, paras_est, paras_est_se, seed,
                     init, scale_mcmc, LL_spec)

names(results_mcmc) <- c("out_mcmc", "paras_est", "para_est_se",
                         "seed", "init", "scale_mcmc", "LL_spec")

# --------------------------- save results ------------------------------------

save(results_mcmc, file = paste0(path_out, outfile))