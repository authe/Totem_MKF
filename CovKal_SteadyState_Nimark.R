# first created: 5 Oct 2023
# last updated: 5 Oct 2023
# author: Andreas Uthemann

covkal_ss_nimark <- function(A, C, D1, D2, R, P0, tol=1e-15) {

    # function returns steady state Kalman gains and covariance matrix
    # for state space system:
    # theta_t = A theta_(t-1) + C (w_t eta_{j,t})' with w_t = (u_t, e_t)
    # weak submitters:
    # y_t = D1 theta_t + D2 theta_{t-1} + R (w_t eta_{j,t})'
    # strong submitters: 
    # y_t = D1 theta_t + D2 theta_{t-1} + R (w_t eta_{j,t})'

    # for formulas used for Kalman filter see:
    # Nimark (2015) "A Low Dimensional Kalman Filter", Economic Letters

    #--------------------------------------------------------------------------

    # update Kalman gain Kt and covariance matrix Pt
    # until convergence to stationary values
    count_max <- 1000
    test_convergence <- TRUE
    count <- 0
    P <- P0

    while (test_convergence & count <= count_max ){
        K <- ( A %*% P %*% t(D1 %*% A + D2) + C %*% t(C) %*% t(D1) + C %*% t(R)) %*% 
            solve( (D1 %*% A + D2) %*% P %*% t(D1 %*% A + D2) + (D1 %*% C + R) %*% t(D1 %*% C + R))

        Pnext <- A %*% P %*% t(A) + C %*% t(C) - 
            K %*% ( (D1 %*% A + D2) %*% P %*% t(D1 %*% A + D2) + (D1 %*% C + R) %*% t(D1 %*% C + R)) %*% t(K)

        # matrix distance between PP and P
        dist <- max(abs(Pnext - P))
        test_convergence <- (dist > tol)

        P <- Pnext
        count <- count + 1

    }

    ss <- list()
    ss$"cov" <- P # return Var( alpha_t | I_t)
    ss$"kg" <- K

    return(ss)
}