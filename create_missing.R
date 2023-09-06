N <- 40*200
miss <- rbinom(N, 1, 0.01)
miss <- matrix(as.logical(miss), nrow = 40)

y_miss <- y_sim

for (s in 1:40){
  y_miss[miss[, s], (s + 1)] <- NA
}


# calculate log-likelihood at true parameters (mostly test)
M_mkf <- 50000 # types draws for MKF
p_class <- 0.9
tic()
LL_miss <- LogLike_MKF(y = y_miss, paras = paras_sim, ord = ord_sim, 
                  p_class = p_class, M = M_mkf)
toc()