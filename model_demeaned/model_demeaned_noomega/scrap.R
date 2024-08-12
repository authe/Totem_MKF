S <- 100         # number of submitters
W <- 66         # number of weak submitters
TT <- 200       # submission dates
rho <- 0.94      # persistence of fundamental theta_t^(0) = rho theta_(t-1)^(0) + u_t
omega <- 0.66    # population share of weak submitters
sig_u <- 0.09  # sd of fundamental shock u_t
sig_e <- 0.05   # sd of noise shock e_t to public signal p_t = (1-omega) theta_{t-1}^(0) + omega theta_{t-1}^(1) + e_t
sig_n <- 0.2   # sd of noise shock n_(j,t) to private signal of weak submitters j: s_(j,t) = theta_t^(0) + n_(j,t)
sig_z <- 0.01      # sd of measurement error for individual submissions
ord <- 2       # order of weak submitters' average belief hierachy theta_t = (theta_t^(0),theta_t^(1),...,theta_t^(ord))
seed=1
tol=1e-15

pc <- prcomp(submissions_demeaned,
             center = TRUE,
             scale. = TRUE)
attributes(pc)

PC1 <- pc$x[, 1]

facreg <- function(x){
  aux <- lm(x ~ PC1)
  b1 <- aux$coefficients[2]
  ar_1 <- ar.ols(aux$residuals, order.max = 3)
  ar_ord <- ar_1$order
  ar_sd <- sqrt(ar_1$var.pred)
  ret <- list(b1 = b1, ar_ord = ar_ord, ar_sd = ar_sd)
  return(ret)
}

reg <- apply(submissions_demeaned, 2, facreg)