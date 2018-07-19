# pf used in pmmh inference
#'@rdname pf
#'@title bootstrap particle filter inference
#'@description runs a particle filter with specified model, model parameters and number of particles
#'
#'@param y observation set
#'@param theta model parameters
#'@param model initial and transition densities for HMM
#'@param nparticles number of particles
#'@export
pf <- function(y, theta, model, nparticles,ess_threshold=1,systematic=FALSE){

  nobservations <- nrow(y)
  dimension <- model$dimension
  Tree <- new(TreeClass, nparticles, 10*nparticles*dimension, dimension)

  # initial step
  randomness <- model$generate_randomness(nparticles, nobservations)
  xparticles <- model$rinit(nparticles, theta, randomness)
  Tree$init(t(xparticles))
  normweights <- rep(1/nparticles, nparticles)
  #
  loglik <- rep(0, nobservations)
  log_W <- rep(log(1/nparticles),nparticles)
  logw <- rep(0,nparticles)
  pf_means <- matrix(0, nrow = nobservations, ncol = dimension)

  for (time in 1:nobservations){
    xparticles <- model$rtransition(xparticles, theta, time, randomness,nparticles)
    logw <- model$dmeasurement(xparticles, theta, y[time,])

    log_weights <- logw + log_W
    maxlw <- max(log_weights)
    w <- exp(log_weights - maxlw)
    loglik[time] <- maxlw + log(sum(w))
    normweights <- w / sum(w)
    # filtering means
    pf_means[time,] <- apply(X = xparticles, MARGIN = 2, FUN = function(x) sum(normweights*x))
    #
    if(((1/sum(normweights**2))/nparticles)<ess_threshold){
      if(!systematic){
        ancestors <- sample(x = 1:nparticles, nparticles, replace = TRUE, prob = normweights)
      }else{
        ancestors <- systematic_resampling_n(normweights, nparticles, runif(1))
      }

      log_W <- rep(log(1/nparticles),nparticles)
      logw <- rep(0,nparticles)
    }else{
      ancestors = 1:nparticles
      log_W <- log(normweights)
    }
    xparticles <- xparticles[ancestors,,drop=F]
    Tree$update(t(xparticles), ancestors - 1)
    #
  }
  index_path <- sample(1:nparticles, 1, prob = normweights)
  path <- Tree$get_path(index_path - 1)
  return(list(pf_means = pf_means, loglik = sum(loglik), path = t(path),loglik_t =cumsum(loglik)))
}


#'@rdname CPF
#'@title Conditional Particle Filter
#'@description Runs a conditional particle filter, with or without ancestor sampling
#'@param nparticles number of particles
#'@param model a list representing a model, for instance as given by \code{\link{get_ar}}.
#'@param theta a parameter to give to the model functions
#'@param observations a matrix of observations of size datalength x dimension(observation)
#'@param ref_trajectory a reference trajectory, of size dimension(process) x datalength; if missing, runs a standard
#'particle filter.
#'@param with_as whether ancestor sampling should be used (TRUE/FALSE)
#'@return A new trajectory.
#'@export
CPF <- function(observations, theta, model, nparticles, ref_trajectory = NULL, with_as = FALSE){
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  xparticles <- model$rinit(nparticles, theta, model$rinit_rand(nparticles, theta), model_precomputed)
  if (!is.null(ref_trajectory)){
    xparticles[,nparticles] <- ref_trajectory[,1]
  }
  Tree$init(xparticles)
  #
  normweights <- rep(1/nparticles, nparticles)
  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    x_last <- xparticles
  }
  # step t > 1
  for (time in 1:datalength){
    ancestors <- multinomial_resampling_n(normweights, nparticles)
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors <- 1:nparticles
    }
    xparticles <- xparticles[,ancestors]
    if (is.null(dim(xparticles))) xparticles <- matrix(xparticles, nrow = dimension)
    xparticles <- model$rtransition(xparticles, theta, time, model$rtransition_rand(nparticles, theta), model_precomputed)
    if (!is.null(ref_trajectory)){
      xparticles[,nparticles] <- ref_trajectory[,time+1]
      if (with_as){
        # Ancestor sampling
        logm <- model$dtransition(ref_trajectory[,time+1], x_last, theta, time, model_precomputed)
        logm <- log(normweights) + logm
        w_as <- exp(logm - max(logm))
        w_as <- w_as / sum(w_as)
        ancestors[nparticles] <- systematic_resampling_n(w_as, 1, runif(1))
        x_last <- xparticles
      } else {
        ancestors[nparticles] <- nparticles
      }
    }
    #
    logw <- model$dmeasurement(xparticles, theta, observations[time,], model_precomputed)
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    normweights <- w / sum(w)
    #
    Tree$update(xparticles, ancestors - 1)
  }
  new_trajectory <- Tree$get_path(systematic_resampling_n(normweights, 1, runif(1))-1)
  return(new_trajectory)
}


#'@rdname CPF_coupled
#'@title Coupled Conditional Particle Filter
#'@description Runs a coupled conditional particle filter, with or without ancestor sampling
#'@param nparticles number of particles
#'@param model a list representing a model, for instance as given by \code{\link{get_ar}}.
#'@param theta a parameter to give to the model functions
#'@param observations a matrix of observations of size datalength x dimension(observation)
#'@param ref_trajectory1 a first reference trajectory, of size dimension(process) x datalength
#'@param ref_trajectory2 a second reference trajectory, of size dimension(process) x datalength
#'@param coupled_resampling a coupled resampling scheme, such as \code{\link{CR_indexmatching}}.
#'@param with_as whether ancestor sampling should be used (TRUE/FALSE)
#'@return A pair of new trajectories.
#'@export
CPF_coupled <- function(nparticles, model, theta, observations, ref_trajectory1, ref_trajectory2,
                        coupled_resampling, with_as = FALSE){
  #
  datalength <- nrow(observations)
  # create tree representation of the trajectories
  Tree1 <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  Tree2 <- new(TreeClass, nparticles, 10*nparticles*model$dimension, model$dimension)
  # initialization
  model_precomputed <- model$precompute(theta)
  init_rand <- model$rinit_rand(nparticles, theta)
  xparticles1 <- model$rinit(nparticles, theta, init_rand, model_precomputed)
  xparticles1[,nparticles] <- ref_trajectory1[,1]
  Tree1$init(xparticles1)
  normweights1 <- rep(1/nparticles, nparticles)
  #
  xparticles2 <- model$rinit(nparticles, theta, init_rand, model_precomputed)
  xparticles2[,nparticles] <- ref_trajectory2[,1]
  Tree2$init(xparticles2)
  normweights2 <- rep(1/nparticles, nparticles)
  #
  # if ancestor sampling, needs to keep the last generation of particles at hand
  if (with_as){
    x_last1 <- xparticles1
    x_last2 <- xparticles2
  }
  # step t > 1
  for (time in 1:datalength){
    ancestors <- coupled_resampling(xparticles1, xparticles2, normweights1, normweights2)
    ancestors1 <- ancestors[,1]
    ancestors2 <- ancestors[,2]
    # if no observation or first time, no resampling
    if (time == 1 || (time > 1 && is.na(observations[time-1,1]))){
      ancestors1 <- 1:nparticles
      ancestors2 <- 1:nparticles
    }
    #
    xparticles1 <- xparticles1[,ancestors1]
    xparticles2 <- xparticles2[,ancestors2]

    if (is.null(dim(xparticles1))) xparticles1 <- matrix(xparticles1, nrow = dimension)
    if (is.null(dim(xparticles2))) xparticles2 <- matrix(xparticles2, nrow = dimension)
    #
    transition_rand <- model$rtransition_rand(nparticles, theta)
    xparticles1 <- model$rtransition(xparticles1, theta, time, transition_rand, model_precomputed)
    xparticles2 <- model$rtransition(xparticles2, theta, time, transition_rand, model_precomputed)
    if (is.null(dim(xparticles1))) xparticles1 <- matrix(xparticles1, nrow = dimension)
    if (is.null(dim(xparticles2))) xparticles2 <- matrix(xparticles2, nrow = dimension)
    #
    xparticles1[,nparticles] <- ref_trajectory1[,time+1]
    xparticles2[,nparticles] <- ref_trajectory2[,time+1]
    if (with_as){
      # % Ancestor sampling
      logm1 <- model$dtransition(ref_trajectory1[,time+1], x_last1, theta, time, model_precomputed)
      logm1 <- log(normweights1) + logm1
      w_as1 <- exp(logm1 - max(logm1))
      w_as1 <- w_as1 / sum(w_as1)
      unif_resampling_as <- runif(1)
      ancestors1[nparticles] = systematic_resampling_n(w_as1, 1, unif_resampling_as)
      x_last1 <- xparticles1
      #
      logm2 <- model$dtransition(ref_trajectory2[,time+1], x_last2, theta, time, model_precomputed)
      logm2 <- log(normweights2) + logm2
      w_as2 <- exp(logm2 - max(logm2))
      w_as2 <- w_as2 / sum(w_as2)
      ancestors2[nparticles] = systematic_resampling_n(w_as2, 1, unif_resampling_as)
      x_last2 <- xparticles2
    } else {
      ancestors1[nparticles] <- nparticles
      ancestors2[nparticles] <- nparticles
    }
    #
    logw1 <- model$dmeasurement(xparticles1, theta, observations[time,], model_precomputed)
    logw2 <- model$dmeasurement(xparticles2, theta, observations[time,], model_precomputed)
    #
    maxlw1 <- max(logw1)
    w1 <- exp(logw1 - maxlw1)
    normweights1 <- w1 / sum(w1)
    #
    maxlw2 <- max(logw2)
    w2 <- exp(logw2 - maxlw2)
    normweights2 <- w2 / sum(w2)
    #
    Tree1$update(xparticles1, ancestors1 - 1)
    Tree2$update(xparticles2, ancestors2 - 1)
  }
  u <- runif(1)
  k_path1 <- systematic_resampling_n(normweights1, 1, u)
  k_path2 <- systematic_resampling_n(normweights2, 1, u)
  ##
  new_trajectory1 <- Tree1$get_path(k_path1 - 1)
  new_trajectory2 <- Tree2$get_path(k_path2 - 1)
  return(list(new_trajectory1 = new_trajectory1, new_trajectory2 = new_trajectory2))
}




# kalman filter
#'@rdname kf
#'@title kalman filter
#'@description runs a kalman filter with specified model, model parameters
#'
#'@param y observation set
#'@param theta model parameters
#'@param mu_0 initial mean
#'@param Sigma_0 initial variance
#'@export
kf <- function(y, theta, mu_0, Sigma_0){
  nobservations <- nrow(y)
  dimension <- ncol(y)
  if(dimension>1){
    kf_means <- matrix(0, ncol = dimension, nrow = nobservations+1)
    m_current <- mu_0
    kf_means[1,] <- m_current
    V_current <- Sigma_0
    Sigma_W <- diag(theta[2]^2, dimension, dimension)
    Sigma_V <- diag(theta[3]^2, dimension, dimension)
    Id <- diag(1, dimension, dimension)
    A <- diag(1, dimension, dimension)
    Phi <- diag(theta[1], dimension, dimension)
    loglik <- rep(0, nobservations)
    for (t in 1:nobservations){
      # prediction step
      m_next <- Phi %*% matrix(m_current, ncol = 1)
      V_next <- Phi %*% V_current %*% t(Phi) + Sigma_W
      # likelihood calculation
      loglik[t] <- fast_dmvnorm(y[t,,drop=F], A %*% m_next, A %*% V_next %*% t(A) + Sigma_V)
      # update step
      K <- V_next %*% t(A) %*% solve(A %*% V_next %*% t(A) + Sigma_V)
      m_current <- m_next + K %*% (y[t,] - A %*% m_next)
      V_current <- (Id - K %*% A) %*% V_next
      kf_means[t+1,] <- m_current
    }
    return(list(kf_means = kf_means[2:(nobservations+1),,drop=F], loglik = sum(loglik)))
  }else{
    return(kf_1d(y, theta, mu_0, Sigma_0))
  }
}


### Kalman filter
# kalman filter
#'@rdname kf_1d
#'@title kalman filter (efficient 1D implementation)
#'@description runs a kalman filter with specified model, model parameters
#'
#'@param y observation set
#'@param theta model parameters
#'@param mu_0 initial mean
#'@param Sigma_0 initial variance
#'@export
kf_1d <- function(y, theta, mu_0, Sigma_0){
  nobservations <- nrow(y)
  dimension <- ncol(y)
  kf_means <- matrix(0, ncol = dimension, nrow = nobservations+1)
  m_current <- mu_0
  kf_means[1,] <- m_current
  V_current <- Sigma_0
  Sigma_W <- theta[2]^2
  Sigma_V <- theta[3]^2
  Id <- 1
  A <- 1
  Phi <- theta[1]
  loglik <- rep(0, nobservations)
  for (t in 1:nobservations){
    # prediction step
    m_next <- Phi * m_current
    V_next <- Phi * V_current * Phi + Sigma_W
    # likelihood calculation
    loglik[t] <- dnorm(y[t,,drop=F], A * m_next,sd=( A * V_next * A + Sigma_V)^0.5,log=T)
    # update step
    K <- V_next * A / (A * V_next * A + Sigma_V)
    m_current <- m_next + K * (y[t,] - A * m_next)
    V_current <- (Id - K * A) * V_next
    kf_means[t+1,] <- m_current
  }
  return(list(kf_means = kf_means[2:(nobservations+1),,drop=F], loglik = sum(loglik)))
}



