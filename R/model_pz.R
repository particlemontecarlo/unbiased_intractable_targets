#'@rdname get_pz
#'@title Phytoplankton-zooplankton model as in Jones, Parslow, Murray 2010
#'@description This function returns a list with objects such as
#'* rinit, rinit_rand to sample from the initial distribution
#'* rtransition, rtransition_rand to sample from the transition
#'* dtransition to evaluate the transition density
#'* dmeasurement to evaluate the measurement density
#'* dimension, which represents the dimension of the latent process
#'@return A list
#'@export
# PZ model
# theta = (mu_alpha, sd_alpha, c, e, ml, mq)
# transformed so that all parameters are in R
# theta = (mu_alpha, log sd_alpha, logit c, logit e, logit ml, logit mq)
get_pz <- function(){
  # logit <- function(z) log(z / (1 - z))
  expit <- function(z) 1 / (1 + exp(-z))

  generate_randomness <- function(nparticles, nobservations){
    return(matrix(rnorm((nobservations+2)*nparticles),nrow=nobservations+2,ncol=nparticles))
  }

  rinit <- function(nparticles, theta, rand, ...){
    return(exp(matrix(log(2) + as.numeric(rand[1:2,]), nrow = nparticles)))
  }

  rtransition <- function(xparticles, theta, time, rand, nparticles, ...){
    alphas <- theta[1] + exp(theta[2]) * rand[time+2,]
    xparticles <- pz_transition(xparticles, alphas, time-1, expit(theta[3:6]))
    return(xparticles)
  }
  #
  dmeasurement <- function(xparticles, theta, observation, ...) {
    return(dnorm(x = observation, mean = log(xparticles[,1]), sd = 0.2, log = TRUE))
  }
  #
  pz_model <- list(rinit = rinit, rtransition = rtransition,
                   generate_randomness = generate_randomness,
                   dmeasurement = dmeasurement, dimension = 2)
  # model <- list(rinit = rinit, rtransition = rtransition,dtransition =dtransition,
  #               dmeasurement = dmeasurement, generate_randomness = generate_randomness,
  #               dimension = dimension)
  return(pz_model)
}





#'@rdname logit
#'@title logit
#'@description logit
#'@return A vector or matrix
#'@export
logit <- function(z) log(z / (1 - z))

#'@rdname expit
#'@title expit
#'@description expit
#'@return A vector or matrix
#'@export
expit <- function(z) 1 / (1 + exp(-z))


#'@rdname pz_params2transformed
#'@title PZ-parameter conversion
#'@description Converts params to R^D
#'@return A vector or matrix
#'@export
pz_params2transformed <- function(x){
  if(is.vector(x)){
    return(c(x[1],log(x[2]),logit(x[3:4])))
  }else{
    return(cbind(x[,1],log(x[,2]),logit(x[,3:6])))
  }
}

#'@rdname pz_transformed2params
#'@title PZ-parameter conversion
#'@description Converts from R^D to the model
#'@return A vector or matrix
#'@export
pz_transformed2params <- function(x){
  if(is.vector(x)){
    return(c(x[1],exp(x[2]),expit(x[3:4])))
  }else{
    return(cbind(x[,1],exp(x[,2]),expit(x[,3:6])))
  }
}




#'@rdname get_pz_naturalparams
#'@title Phytoplankton-zooplankton model as in Jones, Parslow, Murray 2010
#'@description The same as get_pz though implemented with natural params (without the transformation)
#'@return A list
#'@export
# PZ model
# theta = (mu_alpha, sd_alpha, c, e, ml, mq)
get_pz_naturalparams <- function(){
  # logit <- function(z) log(z / (1 - z))
  expit <- function(z) 1 / (1 + exp(-z))

  generate_randomness <- function(nparticles, nobservations){
    return(matrix(rnorm((nobservations+2)*nparticles),nrow=nobservations+2,ncol=nparticles))
  }

  rinit <- function(nparticles, theta, rand, ...){
    return(exp(matrix(log(2) + as.numeric(rand[1:2,]), nrow = nparticles)))
  }

  rtransition <- function(xparticles, theta, time, rand, nparticles, ...){
    alphas <- theta[1] + (theta[2]) * rand[time+2,]
    xparticles <- pz_transition(xparticles, alphas, time-1, (theta[3:6]))
    return(xparticles)
  }
  #
  dmeasurement <- function(xparticles, theta, observation, ...) {
    return(dnorm(x = observation, mean = log(xparticles[,1]), sd = 0.2, log = TRUE))
  }
  #
  pz_model <- list(rinit = rinit, rtransition = rtransition,
                   generate_randomness = generate_randomness,
                   dmeasurement = dmeasurement, dimension = 2)
  # model <- list(rinit = rinit, rtransition = rtransition,dtransition =dtransition,
  #               dmeasurement = dmeasurement, generate_randomness = generate_randomness,
  #               dimension = dimension)
  return(pz_model)
}


