#'@rdname get_lgssm_2params
#'@title Linear-gaussian state-space model (2 params)
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
#'@export
# LG-SSM model
# theta = (a, sigma_X, sigma_Y)
get_lgssm_2params <- function(dimension,sigma_y){
  #
  generate_randomness <- function(nparticles, nobservations){
    return(matrix(rnorm((nobservations+1)*nparticles*dimension),nrow=nobservations+1,ncol=nparticles*dimension))
  }
  rinit <- function(nparticles, theta, randomness, ...){
    return(matrix(randomness[1,], nrow = nparticles, ncol = dimension))
  }
  rtransition <- function(xparticles, theta, time, randomness,nparticles, ...){
    # time is going from 1...T
    Phi <- diag(theta[1], dimension, dimension)
    noise <- matrix(randomness[time+1,], nrow = nparticles, ncol = dimension)
    return(xparticles %*% Phi + noise %*% diag(theta[2], dimension, dimension))
  }
  dtransition <- function(next_x, xparticles, theta, time, ...){
    Phi <- diag(theta[1], dimension, dimension)
    return(fast_dmvnorm_chol_inverse(xparticles %*% Phi, next_x, diag(1/theta[2], dimension, dimension)))
  }
  #
  dmeasurement <- function(xparticles, theta, observation, ...){
    return(fast_dmvnorm_chol_inverse(xparticles, observation, diag(1/sigma_y, dimension, dimension)))
    # return(fast_dmvnorm(xparticles, observation, diag(theta[3]^2, dimension, dimension)))
  }

  #
  #  also add rprior, dprior
  model <- list(rinit = rinit, rtransition = rtransition,dtransition =dtransition,
                dmeasurement = dmeasurement, generate_randomness = generate_randomness,
                dimension = dimension)

  return(model)
}
