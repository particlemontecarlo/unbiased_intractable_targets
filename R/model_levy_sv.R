#'@export
get_levydriven <- function(){
  model <- list(thetadim = 5, ydim = 1)
  # log prior
  model$dprior <- function(thetas){
    if (is.null(dim(thetas))) thetas <- matrix(thetas, nrow = 1)
    density_evals <- dnorm(thetas[,1], mean = 0, sd = sqrt(2), log = TRUE)
    density_evals <- density_evals + dnorm(thetas[,2], mean = 0, sd = sqrt(2), log = TRUE)
    density_evals <- density_evals + dexp(thetas[,3], rate = 0.2, log = TRUE)
    density_evals <- density_evals + dexp(thetas[,4], rate = 0.2, log = TRUE)
    density_evals <- density_evals + dexp(thetas[,5], rate = 1, log = TRUE)
    return(density_evals)
  }
  #
  model$rinit <- function(nparticles, theta){
    x <- matrix(nrow = nparticles, ncol = 2)
    for (i in 1:nparticles){
      x[i,] <- rgamma(2, shape = theta[3] * theta[3]/theta[4], scale = theta[4]/theta[3])
    }
    return(x)
  }
  model$rtransition <- levydriven_rtrans_
  model$dobs <- function(observations, time, xparticles, theta){
    return(dnorm(observations[time], mean = theta[1] + theta[2] * xparticles[,1], sd = sqrt(xparticles[,1]), log = TRUE))
  }
  return(model)
}

#'@export
levydriven_pf <- function(nparticles, model, theta, observations){
  theta2overtheta3 <- theta[3] / theta[4]
  theta4theta2theta2overtheta3 <- theta[5] * theta[3] * theta[3] / theta[4]
  expminustheta4 <- exp(-theta[5])
  oneovertheta4 <- 1/theta[5]
  thetatransform <- c(theta2overtheta3, theta4theta2theta2overtheta3, expminustheta4, oneovertheta4)
  datalength <- nobservations
  # initialization
  xparticles <- model$rinit(nparticles, theta)
  logw <- rep(0, nparticles)
  maxlw <- max(logw)
  w <- exp(logw - maxlw)
  # update log likelihood estimate
  ll <- maxlw + log(mean(w))
  normweights <- w / sum(w)
  # step t > 1
  for (time in 1:datalength){
    if (time > 1){
      ancestors <- systematic_resampling_n(normweights, nparticles, runif(1))
      xparticles <- xparticles[ancestors,,drop=F]
    }
    # xparticles <- model$rtransition(xparticles, theta)
    model$rtransition(xparticles, theta, thetatransform)
    logw <- model$dobs(observations, time, xparticles, theta)
    if (all(is.infinite(logw))){
      return(NA)
    }
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    # update log likelihood estimate
    ll <- ll + maxlw + log(mean(w))
    normweights <- w / sum(w)
  }
  return(ll)
}


#' #'@rdname get_levy_sv
#' #'@title Levy-driven stochastic volatility model as in Barndorff-Nielsen and Shephard (2001)
#' #'@description This function returns a list with objects such as
#' #'* rinit, rinit_rand to sample from the initial distribution
#' #'* rtransition, rtransition_rand to sample from the transition
#' #'* dtransition to evaluate the transition density
#' #'* dmeasurement to evaluate the measurement density
#' #'* dimension, which represents the dimension of the latent process
#' #'@return A list
#' #'@export
#' # SV model
#' # theta = (mu, beta, xi, w2, lambda)
#' # transformed so all are in R...
#' get_model_SVLevy_singlefactor <- function(){
#'   model = list()
#'   # Dimension of parameter, observations, and possibly latent states (int)
#'   model$dimtheta = 5
#'   model$dimY = 1
#'   model$dimX = 2
#'   model$dimension <- 2
#'
#'   #----------------------------------------------------------------------------------------------------
#'   #----------------------------------------------------------------------------------------------------
#'   # Note: if no likelihood nor predictive is provided, the method will be SMC2, which requires
#'   # specifying the transition kernel and the observation density
#'   #----------------------------------------------------------------------------------------------------
#'   #----------------------------------------------------------------------------------------------------
#'   # Sampler from the initial distribution of the latent states
#'   # inputs: theta (single vector), Nx (int)
#'   # outputs: matrix (dimX by Nx) of latent states
#'   model$generate_randomness <- function(nparticles, nobservations){
#'     return(NULL)
#'   }
#'
#'   model$rinitial = function(Nx, theta_transformed, rand, ...){
#'     # xi = theta[3]
#'     # w2 = theta[4]
#'     # lambda = theta[5]
#'     # Note: to avoid numerical issues, we artificially set the variance vt to machine epsilon
#'     # instead of 0 whenever needed
#'     theta <- transformed2params(theta_transformed)
#'
#'     Xs = debiasedpmcmc:::rinitial_SVLevy_cpp(Nx, 1, theta[3], theta[4], theta[5])
#'     if (all(Xs[,1]<.Machine$double.eps)) {Xs[,1] = rep(.Machine$double.eps, Nx)}
#'     if (all(Xs[,2]<.Machine$double.eps)) {Xs[,2] = rep(.Machine$double.eps, Nx)}
#'     return (Xs)
#'   }
#'   # Sampler from the transition distribution of the latent states
#'   # inputs: current states Xs at time (t-1) (dimX by Nx matrix), time t (int), theta (single vector)
#'   # outputs: updated states (dimX by Nx)
#'   model$rtransition = function(Xs, theta_transformed, t, rand, nparticles, ...){
#'     # xi = theta[3]
#'     # w2 = theta[4]
#'     # lambda = theta[5]
#'     # Note: to avoid numerical issues, we artificially set the variance vt to machine epsilon
#'     # instead of 0 whenever needed
#'     theta <- transformed2params(theta_transformed)
#'
#'     new_Xs = debiasedpmcmc:::rtransition_SVLevy_cpp(Xs, t-1, t, theta[3], theta[4], theta[5])
#'     if (all(new_Xs[,1]<.Machine$double.eps)) {new_Xs[,1] = rep(.Machine$double.eps, nparticles)}
#'     if (all(new_Xs[,2]<.Machine$double.eps)) {new_Xs[,2] = rep(.Machine$double.eps, nparticles)}
#'     return (new_Xs)
#'   }
#'   # observation density
#'   # inputs: single observation Yt (dimY by 1), states Xts (dimX by Nx), time t, theta (single vector), log (TRUE by default)
#'   # outputs: observation (log)-densities ("vectorized" with respect to the states Xt)
#'   # model$dobs = function(Xts, theta, Yt){
#'   #   # mu = theta[1]
#'   #   # beta = theta[2]
#'   #   return (debiasedpmcmc:::dobs_SVLevy_cpp(Yt, Xts, theta[1], theta[2], TRUE)[1,])
#'   # }
#'   model$dmeasurement <- function(xparticles, theta_transformed, observation, ...) {
#'     theta <- transformed2params(theta_transformed)
#'
#'     return (debiasedpmcmc:::dobs_SVLevy_cpp(matrix(observation, ncol = 1), xparticles, theta[1], theta[2], TRUE)[1,])
#'   }
#'   #
#'   model$robs = function(Xt,theta_transformed){
#'     theta <- transformed2params(theta_transformed)
#'     mu = theta[1]
#'     beta = theta[2]
#'     return (rnorm(1, mu + beta*Xt[1], sd=sqrt(Xt[1])))
#'   }
#'   return(model)
#' }
#'
#'
#' #'@rdname params2transformed
#' #'@title Transform Levy params
#' #'@description Transform native params to R for Levy-driven stochastic volatility model
#' #'@return [x[1:2],log(x[3:5])]
#' #'@export
#' params2transformed <- function(x){
#'   if(is.vector(x)){
#'     return(c(x[1:2],log(x[3:5])))
#'   }else if(length(dim(x))==2){
#'     return(cbind(x[,1:2],log(x[,3:5])))
#'   }
#' }
#'
#'
#' #'@rdname transformed2params
#' #'@title Revert transformed Levy params
#' #'@description Take transformed params and convert to native params for Levy-driven stochastic volatility model
#' #'@return [x[1:2],exp(x[3:5])]
#' #'@export
#' transformed2params <- function(x){
#'   if(is.vector(x)){
#'     return(c(x[1:2],exp(x[3:5])))
#'   }else if(length(dim(x))==2){
#'     return(cbind(x[,1:2],exp(x[,3:5])))
#'   }
#' }
#'
#'
#'
#'
