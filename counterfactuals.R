# first created: 1 Oct 2023
# last updated: 29 Mar 2024
# author: Andreas Uthemann

# function to calculate covariance of beliefs for weak and strong submitters
# and Kalman gains of weak and stong for their signals
# for the counterfactual experiments (1 & 2) and baseline (3):
# 1. no consensus price
# 2. consensus price, but the given submitter does not observe it
# 3. baseline, submitter observes private signal and consensus price
# 4. perfect consensus price

covkal_counterfactual <- function(paras, ord, tol=1e-15) {

    source("createStateSpaceMat2_Lagged_pubpriv_heterogenous.R")
    source("CovKal_SteadyState.R")
    source("CovKal_SteadyState_Nimark.R")

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
    kg_weak_noprice <- aux$kg


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
    kg_strong_noprice <- aux$kg


    # ---------------- 2. submitter does not participate  ---------------------

    # get M and N (see above) for SS model with consensus price
    ss <- StateSpaceMatLag(ord = (ord-1), rho = rho, omega = omega, sig_u = sig_u,
                           sig_e = sig_e, sig_n = sig_n, tol = tol)
    M <- ss$M
    N <- ss$N

    # 2a. weak submitter
    # signal: s_{i,t} = theta_t + sig_n eta_{i,t}
    Z <- cbind(1, t(rep(0, ord-1)))
    H <- sig_n^2
    P0 <- sig2 * diag(ord)

    # covariance of weak's beliefs for theta_t^(0:ord) is aux$cov
    aux <- covkal_ss(M = M, N = N, Z = Z, H = H, P0 = P0, tol = tol)

    # average expectation: thetabar_t = (1-omega) theta_t + omega theta_t^1
    # covariance matrix for weak's beliefs for (theta_t, thetabar_t)^T :
    selec <- rbind(c(1, rep(0, ord-1)), c(1 - omega, omega, rep(0, ord - 2)))
    cov_weak_price <- selec %*% aux$cov %*% t(selec)
    kg_weak_price <- aux$kg

    # 2b. strong submitter

    # signal: s_{i,t} = theta_t
    Z <- cbind(1, t(rep(0, ord-1)))
    H <- 0  # perfect signal
    P0 <- sig2 * diag(ord)

    # covariance of strong's beliefs for theta_t^(0:ord) is aux$cov
    aux <- covkal_ss(M = M, N = N, Z = Z, H = H, P0 = P0, tol = tol)

    # average expectation: thetabar_t = (1-omega) theta_t + omega theta_t^1
    # covariance matrix for weak's beliefs for (theta_t, thetabar_t)^T :
    selec <- rbind(c(1, rep(0, ord-1)), c(1 - omega, omega, rep(0, ord - 2)))
    cov_strong_price <- selec %*% aux$cov %*% t(selec)
    kg_strong_price <- aux$kg


    # ------------------------ 3. baseline model  -----------------------------

    # 3a. weak submitter (with private signal & consensus price)
    A <- ss$M
    C <- ss$N
    C <- cbind(C, rep(0, ord))

    D1 <- rbind(c(1, rep(0, ord-1)), rep(0, ord))
    D2 <- rbind(rep(0, ord), c(1 - omega, omega, rep(0, ord - 2)))

    R <- rbind(c(0, 0, sig_n), c(0, sig_e, 0))  # noisy private signal

    P0 <- sig2 * diag(ord)

    # covariance of weak's beliefs for theta_t^(0:ord) is aux$cov
    aux <- covkal_ss_nimark(A = A, C = C, D1 = D1, D2 = D2, R = R,
                             P0 = P0, tol = tol)

    selec <- rbind(c(1, rep(0, ord-1)), c(1 - omega, omega, rep(0, ord - 2)))
    cov_weak <- selec %*% aux$cov %*% t(selec)
    kg_weak <- aux$kg

    # 3b. strong submitter (with private signal & consensus price)
    R <- rbind(c(0, 0, 0), c(0, sig_e, 0))  # perfect private signal

    # covariance of strong's beliefs for theta_t^(0:ord) is aux$cov
    aux <- covkal_ss_nimark(A = A, C = C, D1 = D1, D2 = D2, R = R,
                             P0 = P0, tol = tol)

    selec <- rbind(c(1, rep(0, ord-1)), c(1 - omega, omega, rep(0, ord - 2)))
    cov_strong <- selec %*% aux$cov %*% t(selec)
    kg_strong <- aux$kg


    # ------------------------ 4. perfect consensus price ---------------------

    # perfect consenus price: p_t = theta_t

    # Kalman gain k11 for weak
    k11_weak <- sig_u^2 / (sig_u^2 + sig_n^2)

    # beta is cofficient on u_t in average 1st order expectation:
    # thetabar_t = rho theta_{t-1} + [(1-omega) + omega k11] u_t
    beta <- (1-omega) + k11_weak * omega

    # covariance matrix of X  (see paper appendix for definition)
    Cov_X <- rbind(c(sig_u^2, beta * sig_u^2), 
                c(beta * sig_u^2, beta^2 * sig_u^2)) 
    # covariance of X and signal y_{i,t} = s_{i,t} - rho theta_{t-1}
    Cov_XY <- c(sig_u^2, beta * sig_u^2)
    sig2_y <- sig_u^2 + sig_n^2

    # covariance for 1st and 2nd order beliefs of weak
    cov_weak_perfect <- Cov_X - (Cov_XY %*% t(Cov_XY)) / sig2_y
    kg_weak_perfect <- rbind(c(k11_weak, (1 - k11_weak) * rho), 
                    c(beta * k11_weak, (beta + omega) * (1 - k11_weak) * rho))

    # strong know fundamental and average 1st order expectation
    cov_strong_perfect <- matrix(rep(0,4), nrow = 2)
    kg_strong_perfect <- rbind(c(1, 0), c(beta, omega * (1 - k11_weak) * rho))
    

    # -------------------------------------------------------------------------

    # set covariance entries with absolute val below tol level to 0
    
    cov_weak_noprice[abs(cov_weak_noprice) < tol] <- 0
    cov_strong_noprice[abs(cov_strong_noprice) < tol] <- 0

    cov_weak_price[abs(cov_weak_price) < tol] <- 0
    cov_strong_price[abs(cov_strong_price) < tol] <- 0

    cov_weak[abs(cov_weak) < tol] <- 0
    cov_strong[abs(cov_strong) < tol] <- 0

    # -------------------------------------------------------------------------

    res <- list()

    res$"cov_weak_noprice" <- cov_weak_noprice
    res$"kg_weak_noprice" <- kg_weak_noprice

    res$"cov_strong_noprice" <- cov_strong_noprice
    res$"kg_strong_noprice" <- kg_strong_noprice

    res$"cov_weak_price" <- cov_weak_price
    res$"kg_weak_price" <- kg_weak_price

    res$"cov_strong_price" <- cov_strong_price
    res$"kg_strong_price" <- kg_strong_price

    res$"cov_weak" <- cov_weak
    res$"kg_weak" <- kg_weak

    res$"cov_strong" <- cov_strong
    res$"kg_strong" <- kg_strong

    res$"cov_weak_perfect" <- cov_weak_perfect
    res$"kg_weak_perfect" <- kg_weak_perfect

    res$"cov_strong_perfect" <- cov_strong_perfect
    res$"kg_strong_perfect" <- kg_strong_perfect

    return(res)

}