# first created: 1 Oct 2023
# last updated: 4 Oct 2023
# author: Andreas Uthemann

# function to calculate covariance of beliefs for weak and strong submitters
# for the counterfactual experiments (1 & 2) and baseline (3):
# 1. no consensus price
# 2. consensus price, but the given submitter does not observe it
# 3. baseline, submitter observes private signal and consensus price

cov_counterfactual <- function(paras, ord, tol=1e-15) {

    source("createStateSpaceMat2_Lagged_pubpriv_heterogenous.R")
    source("CovKal_SteadyState.R")

    rho <- paras[1]
    omega <- paras[2]
    sig_u <- paras[3]
    sig_e <- paras[4]
    sig_n <- paras[5]
    
    # ------------------- 1. no consensus price  ------------------------------

    # calculate stationary Kalman gain and variance of (1st order) belief
    # for weak submitter
    count_max <- 100
    test_convergence <- TRUE
    count <- 0
    sig2 <- sig_u^2 / (1 - rho^2)

    while (test_convergence && count <= count_max) {

        k <- sig2 / (sig2 + sig_n^2)  # Kalman gain
        sig2_next <- rho^2 * (1 - k) * sig2 + sig_u^2 # recursion for variance

        dist <- abs(sig2_next - sig2)
        test_convergence <- (dist > tol)

        sig2 <- sig2_next
        count <- count + 1

    }

    # SS system for weak submitters' average 1st order beliefs
    # theta_t^(0:1) = M theta_{t-1}^(0:1) + N u_t ; u_t ~ N(0,1)
    M <- rbind(c(rho, 0), c(k * rho, (1 - k) * rho))
    N <- rbind(sig_u, sig_u * k)

    # 1a. weak submitter
    # signal: s_{i,t} = theta_t + sig_n eta_{i,t}
    Z <- cbind(1, 0)
    H <- sig_n^2
    P0 <- sig2 * diag(2)

    # covariance of weak's beliefs for theta_t^(0:1) is aux$cov
    aux <- covkal_ss(M = M, N = N, Z = Z, H = H, P0 = P0, tol = tol)

    # average expectation: thetabar_t = (1-omega) theta_t + omega theta_t^1
    # covariance matrix for weak's beliefs for (theta_t, thetabar_t)^T :
    selec <- rbind(c(1, 0), c(1 - omega, omega))
    cov_weak_noprice <- selec %*% aux$cov %*% t(selec)

    # 1b. strong submitter
    # signal: s_{i,t} = theta_t
    Z <- cbind(1, 0)
    H <- 0  # perfect signal
    P0 <- sig2 * diag(2)

    # covariance of strong's beliefs for theta_t^(0:1) is aux$cov
    aux <- covkal_ss(M = M, N = N, Z = Z, H = H, P0 = P0, tol = tol)

    # average expectation: thetabar_t = (1-omega) theta_t + omega theta_t^1
    # covariance matrix for weak's beliefs for (theta_t, thetabar_t)^T :
    selec <- rbind(c(1, 0), c(1 - omega, omega))
    cov_strong_noprice <- selec %*% aux$cov %*% t(selec) # equal zero (CHECK!)


    # ---------------- 2. submitter does not participate  ---------------------

    # get M and N (see above) for SS model with consensus price
    ss <- StateSpaceMatLag(ord = ord, rho = rho, omega = omega, sig_u = sig_u,
                           sig_e = sig_e, sig_n = sig_n, tol = tol)
    M <- ss$M
    N <- ss$N

    # 2a. weak submitter
    # signal: s_{i,t} = theta_t + sig_n eta_{i,t}
    Z <- cbind(1, t(rep(0, ord)))
    H <- sig_n^2
    P0 <- sig2 * diag(ord + 1)

    # covariance of weak's beliefs for theta_t^(0:1) is aux$cov
    aux <- covkal_ss(M = M, N = N, Z = Z, H = H, P0 = P0, tol = tol)

    # average expectation: thetabar_t = (1-omega) theta_t + omega theta_t^1
    # covariance matrix for weak's beliefs for (theta_t, thetabar_t)^T :
    selec <- rbind(c(1, rep(0, ord)), c(1 - omega, omega, rep(0, ord - 1)))
    cov_weak_price <- selec %*% aux$cov %*% t(selec)

    # 2b. strong submitter

    # signal: s_{i,t} = theta_t
    Z <- cbind(1, t(rep(0, ord)))
    H <- 0  # perfect signal
    P0 <- sig2 * diag(ord + 1)

    # covariance of strong's beliefs for theta_t^(0:1) is aux$cov
    aux <- covkal_ss(M = M, N = N, Z = Z, H = H, P0 = P0, tol = tol)

    # average expectation: thetabar_t = (1-omega) theta_t + omega theta_t^1
    # covariance matrix for weak's beliefs for (theta_t, thetabar_t)^T :
    selec <- rbind(c(1, rep(0, ord)), c(1 - omega, omega, rep(0, ord - 1)))
    cov_strong_price <- selec %*% aux$cov %*% t(selec)


    # ------------------------ 3. baseline model  -----------------------------

    # 3a. weak submitter (with private signal & consensus price)
    selec <- rbind(c(1, rep(0, ord - 1)), c(1 - omega, omega, rep(0, ord - 2)))
    cov_weak <- selec %*% ss$PP[1:2, 1:2] %*% t(selec)

    # 3b. strong submitter (with private signal & consensus price)
    cov_strong <- diag(2)
    cov_strong[1, 1] <- 0  # strong know theta_t^0
    # consensus price has same info as theta_t^1 + (sig_e / omega) e_t
    cov_strong[2, 2] <- sig_e^2


    # -------------------------------------------------------------------------

    res <- list()

    res$"cov_weak_noprice" <- cov_weak_noprice
    res$"cov_strong_noprice" <- cov_strong_noprice

    res$"cov_weak_price" <- cov_weak_price
    res$"cov_strong_price" <- cov_strong_price

    res$"cov_weak" <- cov_weak
    res$"cov_strong" <- cov_strong

    return(res)

}