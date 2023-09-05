# Code for estimating the consensus pricing model via mixture Kalman filter

This repository contains a collection of R programs to simulate and estimate the consensus pricing model developed in "Strategic Uncertainty in OTC Markets: Evidence from a Consensus Pricing Service" (Lerby Ergun and Andreas Uthemann).

The model is estimated by maximum likelihood using MCMC methods with flat priors (metrop algorithm of the "mcmc" R package). The likelihood function is computed using a Mixture Kalman Filter (Chen, R., & Liu, J. S. (2000). Mixture kalman filters. Journal of the Royal Statistical Society: Series B, 62(3), 493-508.).

The logic of the files is as follows:

**test_MKF.R**
Main file to test estimation code and performance
- The code simulates data from the model using **SimulateData2_demeaned_Lagged_pubpriv_heterogenous.R**
- Parameters are then estimated via MCMC (metrop of mcmc package) using the mixture Kalman filter
- Initial conditions for the estimation are set using **InitialParasGuess.R** 

The log-likelihood of the data given parameters is calculated via the mixture Kalman filter (MKF). **LogLike_MKF.R** is the main R function that returns the value of the log-likelihood. **LogLike_MKF_MCMC.R** is an R function that modifies the previous function so that it can be used as an argument of the Metropolis algorithm of the MCMC package. It also uses logit transformations for the parameters so that the metrop algorithm can sample from an unconstrained parameters space.

**LogLike_MKF.R** call a collection of R subroutines to calculate the log-likelihood:
- **SampleTypes.R** samples M draws of submitters types (weak or strong) according to a given binomial prior p_class and returns the resulting matrix of types and the collection of state-space model corresponding to these type draws for given parameters.
- **FKFModel2_demeaned_Lagged_privpub_heterogenous.R** returns the matrices of the state-space model for a given type draw and parameters. **createStateSpaceMat2_Lagged_pubpriv_heterogenous.R** calculates the coefficients of these matrices given parameters. 
- **MixtureKalmanFilter.R** The function MKF returns the log-likelihoods (LL) and posteriors for states (at,Pt) and types (ln.wt) for the type draws (less or equal M because of resampling) given the data (consensus price and submissions). **LogLike_MFK.R** then gives the log-likelihood of the data as a weighted average of the log-likelihoods of the (resampled) draws.
- **KalmanStep.R** Subroutine for the function MKF that updates the likelihood from period t to t+1 given the submissions of period t+1 and the consensus price of t. (Most computationally intensive part of the likelihood)
