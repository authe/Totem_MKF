# first created: 30 Sep 2023
# last updated: 30 Sep 2023
# author: Andreas Uthemann

# function returns steady state Kalman gains and covariance matrix
# for a given state space system:
# alpha_{t+1} = M alpha_t + N eta_t ; eta_t ~ N(0, I)
# y_t = Z alpha_t + epsilon_t; epsilon_t ~ N(0, H)
# with priors alpha_0 ~ N(a0, P0)

covkal_ss <- function(M, N, Z, H, P0, tol=1e-15) {

    NN <- N %*% t(N)

    # update Kalman gain Kt and covariance matrix Pt
    # until convergence to stationary values
    count_max <- 1000
    test_convergence <- TRUE
    count <- 0
    Pt <- P0

    while (test_convergence && count <= count_max) {

        tmp_Pt_trZ <- Pt %*% t(Z)    # avoids repeat calculation of Pt * Zt'
        Ft <- Z %*% tmp_Pt_trZ + H # covariance matrix of forecast errors vt

        # calculate inverse of Ft
        # (via Cholesky decomposition: Ft is symmetric & positive definite)
        cholFt <- chol(Ft)
        invFt <- chol2inv(cholFt)

        Kt <- tmp_Pt_trZ %*% invFt   # Kalman gain
        Ptt <- Pt - tmp_Pt_trZ %*% t(Kt)    # updated variance
        Pnext <- M %*% Ptt %*% t(M) + NN

        # matrix distance between PP and P
        dist <- max(abs(Pnext - Pt))
        test_convergence <- (dist > tol)

        Pt <- Pnext
        count <- count + 1

    }

    ss <- list()
    ss$"cov" <- Pt
    ss$"kg" <- Kt

    return(ss)
}