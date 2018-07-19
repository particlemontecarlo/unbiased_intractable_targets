
# contains experiment settings for LG
lgssm_model_2params <- function(){

  nmcmc <- 500000
  nobservations <- 200
  nobservations_lowT <- 50
  dimension <- 1
  mu_0 <- rep(0, dimension)
  Sigma_0 <- diag(1, dimension, dimension)
  phi <- 0.5
  sigma_x <- 1
  sigma_y <- 1
  theta <- c(phi, sigma_x)

  D_theta <- length(theta)


  # mh_prop <- diag(c(0.01,0.01,0.01),D_theta,D_theta)

  N_arr <- floor(seq(100,500,length.out=20))
  # N_arr_Hkm <- floor(seq(N_arr[1],N_arr[length(N_arr)],length.out = 5))
  N_arr_Hkm <- floor(seq(50,300,50))
  N_arr_serial <- floor(rep(seq(100,300,length.out = 10),5))

  exp_name <- 'lgssm_lowT_calibrated'
  exp_name_serial <- 'lgssm_lowT_calibrated_serial'
  data_file <- "inst/lg_ssm/2params/data/lgssm_model_lowTcalibrated.RData"
  serial_data_file <- "inst/lg_ssm/2params/data/lgssm_serial.RData"
  coupled_data_file <- "inst/lg_ssm/2params/data/lgssm_coupled.RData"
  lowTcalibration_file <-  "inst/lg_ssm/2params/data/lgssm_model_lowTcalibration_calibration_Data.RData"


  # prior function
  logprior <- function(theta){
    # want weakly informative prior on standard dev params
    if (theta[1]>1 || theta[1]< 0 || theta[2]<0){return(-Inf)}
    else{
      theta1_p <- log(0.5)

      shape_theta2 <- 2
      scale_theta2 <- 2
      theta2_p <- dgamma(theta[2], shape_theta2, scale = scale_theta2, log = TRUE)

      return(theta1_p+theta2_p)
    }
  }

  # set init
  rinit <- function(){
    # sample according to the multivariate normal, truncated so the initial conditions
    # are positive for variance parameters
    theta_init <- -Inf
    while(logprior(theta_init)==-Inf){
      theta_init <- mvrnorm(1,mean_est,5*cov_est)
    }
    return(theta_init)
  }



  return(list(nobservations=nobservations,
              nobservations_lowT=nobservations_lowT,
              dimension=dimension,
              mu_0=mu_0,
              Sigma_0=Sigma_0,
              theta=theta,
              sigma_y=sigma_y,
              D_theta=D_theta,
              exp_name=exp_name,
              exp_name_serial=exp_name_serial,
              logprior=logprior,
              N_arr=N_arr,
              N_arr_Hkm=N_arr_Hkm,
              N_arr_serial=N_arr_serial,
              serial_data_file=serial_data_file,
              coupled_data_file=coupled_data_file,
              data_file=data_file,
              lowTcalibration_file=lowTcalibration_file,
              rinit=rinit,
              nmcmc=nmcmc))

}

lgssm_model_2params_priorinit <- function(){
  settings <- lgssm_model_2params()
  # set init
  rinit <- function(){
    return(c(runif(1),rgamma(1,shape=2,scale=2)))
  }

  settings$nmcmc <- 250000
  settings$rinit <- rinit
  settings$N_arr_mt <- seq(75,300,25)
  settings$N_arr_Hkm <- settings$N_arr_mt
  settings$N_arr_serial <- rep(settings$N_arr_mt,10)
  settings$n_rep_Hkm <- 500

  return(settings)
}



lgssm_model_2params_100obs <- function(){
  settings <- lgssm_model_2params_priorinit()
  # set init
  settings$N_arr_mt <- round(seq(50,225,25))
  settings$N_arr_Hkm <- round(seq(50,225,25))
  settings$N_arr_serial <- rep(settings$N_arr_mt,10)
  settings$nobservations <- 100
  settings$nobservations_lowT <- 25
  settings$nmcmc <- 500000
  settings$rinit <- function() {return(c(runif(1),5*runif(1)))}
  settings$coupledres_folder <- 'inst/reproduce/Section 3.1 lgssm/data/'
  settings$nrep_mt <- 20000

  settings$exp_name <- 'lgssm_2params_100obs'
  settings$exp_name_serial <- 'lgssm_2params_100obs_serial_calibrated'
  settings$data_file <- "inst/reproduce/Section 3.1 lgssm/data/lgssm_2params_100obs.RData"
  settings$serial_data_file <- "inst/reproduce/Section 3.1 lgssm/data/lgssm_2params_100obs_serial.RData"
  settings$mt_data_file <- "inst/reproduce/Section 3.1 lgssm/data/lgssm_2params_100obs_mt.RData"

  return(settings)
}






# contains experiment settings for LG
lgssm_model_2params_scalingT <- function(){

  nobservations_highT <- 1000
  nobservations_lowT <- 50
  dimension <- 1
  mu_0 <- rep(0, dimension)
  Sigma_0 <- diag(1, dimension, dimension)
  phi <- 0.5
  sigma_x <- 1
  sigma_y <- 1
  theta <- c(phi, sigma_x)

  D_theta <- length(theta)

  cov_0_prop <- 4*diag(1,D_theta)
  cov_0_init <- 50*diag(1,D_theta)

  # mh_prop <- diag(c(0.01,0.01,0.01),D_theta,D_theta)
  nobs_arr <- seq(nobservations_lowT,nobservations_highT,50)
  N_arr <- 1*nobs_arr

  exp_name <- 'lgssm_scalingT'
  exp_name_serial <- 'lgssm_scalingT_serial'
  data_file <- "inst/lg_ssm/2params/data/lgssm_model_scalingT_data.RData"
  serial_data_file <- "inst/lg_ssm/2params/data/lgssm_model_scalingT_serial.RData"
  mt_resfile <- "inst/lg_ssm/2params/data/lgssm_model_scalingTres.RData"
  lowTcalibration_file <-  "inst/lg_ssm/2params/data/lgssm_model_scalingT_calibration.RData"


  # prior function
  logprior <- function(theta){
    # want weakly informative prior on standard dev params
    if (theta[1]>1 || theta[1]< 0 || theta[2]<0){return(-Inf)}
    else{
      theta1_p <- log(0.5)

      shape_theta2 <- 2
      scale_theta2 <- 2
      theta2_p <- dgamma(theta[2], shape_theta2, scale = scale_theta2, log = TRUE)

      return(theta1_p+theta2_p)
    }
  }


  return(list(nobs_arr=nobs_arr,
              nobservations_lowT=nobservations_lowT,
              nobservations=nobservations_highT,
              dimension=dimension,
              mu_0=mu_0,
              Sigma_0=Sigma_0,
              theta=theta,
              sigma_y=sigma_y,
              D_theta=D_theta,
              exp_name=exp_name,
              exp_name_serial=exp_name_serial,
              logprior=logprior,
              N_arr=N_arr,
              serial_data_file=serial_data_file,
              mt_resfile=mt_resfile,
              data_file=data_file,
              lowTcalibration_file=lowTcalibration_file,
              cov_0_prop=cov_0_prop,
              cov_0_init=cov_0_init))

}


#
#
#
# # contains experiment settings for LG
# lgssm_model_2params_1000obs <- function(){
#
#   nobservations <- 1000
#   nobservations_lowT <- 100
#   dimension <- 1
#   mu_0 <- rep(0, dimension)
#   Sigma_0 <- diag(1, dimension, dimension)
#   phi <- 0.5
#   sigma_x <- 1
#   sigma_y <- 1
#   theta <- c(phi, sigma_x)
#
#   D_theta <- length(theta)
#
#
#   # mh_prop <- diag(c(0.01,0.01,0.01),D_theta,D_theta)
#
#   N_arr <- floor(seq(100,300,length.out=20))
#
#   exp_name <- 'lgssm_lowT_calibrated'
#   exp_name_serial <- 'lgssm_lowT_calibrated_serial'
#   data_file <- "inst/lg_ssm/2params/data/lgssm_model_1000obs.RData"
#   lowTcalibration_file <-  "inst/lg_ssm/2params/data/lgssm_model_1000obs_calibration.RData"
#
#
#   # prior function
#   logprior <- function(theta){
#     # want weakly informative prior on standard dev params
#     if (theta[1]>1 || theta[1]< 0 || theta[2]<0){return(-Inf)}
#     else{
#       theta1_p <- log(0.5)
#
#       shape_theta2 <- 2
#       scale_theta2 <- 2
#       theta2_p <- dgamma(theta[2], shape_theta2, scale = scale_theta2, log = TRUE)
#
#       return(theta1_p+theta2_p)
#     }
#   }
#
#
#   return(list(nobservations=nobservations,
#               nobservations_lowT=nobservations_lowT,
#               dimension=dimension,
#               mu_0=mu_0,
#               Sigma_0=Sigma_0,
#               theta=theta,
#               sigma_y=sigma_y,
#               D_theta=D_theta,
#               exp_name=exp_name,
#               exp_name_serial=exp_name_serial,
#               logprior=logprior,
#               N_arr=N_arr,
#               data_file=data_file,
#               lowTcalibration_file=lowTcalibration_file))
#
# }
#
#
#
#
#
#
# # contains experiment settings for LG
# lgssm_model_flat_prior <- function(){
#
#   nobservations <- 20
#   dimension <- 1
#   mu_0 <- rep(0, dimension)
#   Sigma_0 <- diag(1, dimension, dimension)
#   phi <- 0.9
#   sigma_x <- 1
#   sigma_y <- 0.5
#   theta <- c(phi, sigma_x, sigma_y)
#
#   D_theta <- length(theta)
#
#   posterior_cov <- matrix(c(  0.013893593, -0.0081610438, -0.0031266948,
#                              -0.008161044,  0.0661270457, -0.0002974635,
#                              -0.003126695, -0.0002974635,  0.0870091739),3,3)
#
#
#   #mh_prop<-diag(c(0.1,0.1,0.1),D_theta,D_theta)
#   mh_prop <- 3*posterior_cov
#
#   nparticles_arr <- seq(from=100,to=500,length.out=11)#seq(100,1000,by=100)
#
#   exp_name <- 'lgssm_unbiased_eff_flat_prior'
#   exp_name_serial <- 'lgssm_serial_eff_flat_prior'
#
#
#
#   # prior function
#   logprior <- function(theta){
#     # want weakly informative prior on standard dev params
#     if (theta[1]>1 || theta[1]< -1 || theta[2]<0 || theta[3] <0){return(-Inf)}
#     else{
#       theta1_p <- log(0.5)
#
#       shape_theta2 <- 2
#       scale_theta2 <- 1
#       theta2_p <- dgamma(theta[2], shape_theta2, scale = scale_theta2, log = TRUE)
#
#       shape_theta3 <- 2
#       scale_theta3 <- 1
#       theta3_p <- dgamma(theta[3], shape_theta3, scale = scale_theta3, log = TRUE)
#
#       return(theta1_p+theta2_p+theta3_p)
#     }
#   }
#
#   rinit <- function() c(runif(1,min=0,max=1),abs(rnorm(1)),abs(rnorm(1)))
#
#   return(list(nobservations=nobservations,
#               dimension=dimension,
#               mu_0=mu_0,
#               Sigma_0=Sigma_0,
#               theta=theta,
#               N_arr=nparticles_arr,
#               mh_prop=mh_prop,
#               D_theta=D_theta,
#               exp_name=exp_name,
#               exp_name_serial=exp_name_serial,
#               rinit=rinit,
#               logprior=logprior))
#
# }
#
#
#
#
# pz_settings <- function(){
#
#   nobservations <- 25
#   nobservations_lowT <- 25
#   # untransformed parameters
#   theta_dgp <- c(0.7, 0.5, 0.25, 0.3, 0.1, 0.1)
#   D_theta <- length(theta_dgp)
#
#
#
#   nrep <- 500
#   nrep_initial <- 500
#
#   N_arr <- c(20,40,60,80,100)
#   m <- 150000
#   max_iterations <- 1.01*m
#
#   data_file <- 'inst/pz/pz_obs.RData'
#   lowTcalibration_file <- 'inst/pz/pz_lowT_calibration.RData'
#
#   # pmmh settings
#   nmcmc <- 500000#50000
#
#   # cov_scaling <- 0.5#0.2 (0.2 looks about optimal for serial algorithm for 50 obs)
#   # cov_vec <- c(   0.026690880,  0.002021423,  0.035507876, -0.05882410, -0.123796322,  -0.05294267,
#   #                 0.002021423, 0.020095446 ,-0.001717783,  0.01226273,  0.005164332,   0.03490882,
#   #                 0.035507876, -0.001717783,  1.350985413, -1.02098043,  0.028033578,   0.66215747,
#   #                -0.058824102,  0.012262728, -1.020980432,  1.48925376,  0.564696914,  -0.08484970,
#   #                -0.123796322,  0.005164332,  0.028033578,  0.56469691,  2.600985500,   0.49796597,
#   #                -0.052942667,  0.034908817,  0.662157472, -0.08484970,  0.497965974,   2.26202997)
#   # proposal_covariance <- cov_scaling * structure(cov_vec, .Dim = c(6L, 6L))
#
#   # function to generate data
#   generate_data <- function(nobservations,theta){
#     state <- matrix(exp(log(2) + rnorm(2, mean = 0, sd = 1)), ncol = 2)
#     states <- matrix(ncol = 2, nrow = nobservations+1)
#     states[1,] <- state
#     log_obs <- rep(0, nobservations)
#     for (t in 1:nobservations){
#       alpha <- rnorm(n = 1, mean = theta[1], sd = theta[2])
#       state <- pz_transition(state, alpha, t-1, theta[3:6])
#       states[t+1,] <- state
#       log_obs[t] <- rnorm(1, mean = log(state[1,1]), sd = 0.2)
#     }
#     return(list(states=states,log_obs=log_obs))
#   }
#
#
#   #prior
#   logprior <- function(transformed_theta){
#     ## evaluate prior density on the transformed parameter
#     # normal prior on theta1
#     density_eval <- dnorm(transformed_theta[1], mean = 0, sd = 10, log = TRUE)
#     # exponential prior on theta2
#     density_eval <- density_eval + dexp(exp(transformed_theta[2]), rate = 1, log = TRUE)
#     ## add Jacobian term to account for the transformation
#     # log transform for theta2
#     density_eval <- density_eval + transformed_theta[2]
#     # uniform prior on the other parameters
#     # logit transform
#     for (iparam in 3:6){
#       density_eval <- density_eval - transformed_theta[iparam] - 2 * log(1 + exp(-transformed_theta[iparam]))
#     }
#     return(density_eval)
#   }
#
#   #prior
#   logprior_naturalparams <- function(theta){
#     if(theta[2]<0 | any(theta[3:6]>1) | any(theta[3:6]<0)){
#       return(-Inf)
#     }else{
#       ## evaluate prior density on the transformed parameter
#       # normal prior on theta1
#       density_eval <- dnorm(theta[1], mean = 0, sd = 10, log = TRUE)
#       # exponential prior on theta2
#       density_eval <- density_eval + dexp(theta[2], rate = 1, log = TRUE)
#       # uniform prior on the other parameters
#       # logit transform
#       return(density_eval)
#     }
#   }
#
#
#
#   # posterior density function up to normalising constant
#   estimate_pf_target <- function(transformed_theta,nparticles){
#     log_prior<-logprior(transformed_theta)
#     if(log_prior==-Inf){
#       return(list(log_target=-Inf,path=NA))
#     }else{
#       pf_results <- pf(y, transformed_theta, pz_model,nparticles,ess_threshold = 0.75,systematic = TRUE)
#       log_target <- pf_results$loglik + log_prior
#       path <- pf_results$path
#       return(list(log_target=log_target,path=path))
#     }
#   }
#
# #
# #   rinit <- function(){
# #
# #     init_means <- c(0.44435692, -0.69450502, -0.17233041, -0.08403124,  0.01082659,  0.59587534)
# #     init_vars <-  c(0.01726431, 0.01265800, 0.46045588, 0.86885887, 1.88446711, 1.77972311)
# #     init_std <- 2*init_vars^0.5
# #     return(init_means + init_std*rnorm(D_theta))
# #   }
#
#
#
#   return(list(nobservations=nobservations,
#               theta_dgp=theta_dgp,
#               nmcmc=nmcmc,
#               generate_data=generate_data,
#               logprior=logprior,
#               logprior_naturalparams=logprior_naturalparams,
#               estimate_pf_target=estimate_pf_target,
#               nobservations_lowT=nobservations_lowT,
#               max_iterations=max_iterations,
#               nrep=nrep,
#               nrep_initial=nrep_initial,
#               data_file=data_file,
#               lowTcalibration_file=lowTcalibration_file,
#               N_arr=N_arr,
#               m=m))
#
#
#
# }
#
#
#
# pz_large_data_settings <- function(){
#
#
#
#   nobservations <- 100
#   nobservations_lowT <- 50
#   # untransformed parameters
#   theta_dgp <- c(0.7, 0.5, 0.25, 0.3, 0.1, 0.1)
#   D_theta <- length(theta_dgp)
#
#
#
#   nrep <- 500
#   nrep_initial <- 1000
#
#   N_arr <- c(20,40,60,80,100)
#   m <- 150000
#   max_iterations <- 1.01*m
#
#   data_file <- 'inst/pz/pz_obs_large_data.RData'
#   lowTcalibration_file <- 'inst/pz/pz_lowT_calibration_large_data.RData'
#
#   # pmmh settings
#   nmcmc <- 500000#50000
#
#   # cov_scaling <- 0.5#0.2 (0.2 looks about optimal for serial algorithm for 50 obs)
#   # cov_vec <- c(   0.026690880,  0.002021423,  0.035507876, -0.05882410, -0.123796322,  -0.05294267,
#   #                 0.002021423, 0.020095446 ,-0.001717783,  0.01226273,  0.005164332,   0.03490882,
#   #                 0.035507876, -0.001717783,  1.350985413, -1.02098043,  0.028033578,   0.66215747,
#   #                -0.058824102,  0.012262728, -1.020980432,  1.48925376,  0.564696914,  -0.08484970,
#   #                -0.123796322,  0.005164332,  0.028033578,  0.56469691,  2.600985500,   0.49796597,
#   #                -0.052942667,  0.034908817,  0.662157472, -0.08484970,  0.497965974,   2.26202997)
#   # proposal_covariance <- cov_scaling * structure(cov_vec, .Dim = c(6L, 6L))
#
#   # function to generate data
#   generate_data <- function(nobservations,theta){
#     state <- matrix(exp(log(2) + rnorm(2, mean = 0, sd = 1)), ncol = 2)
#     states <- matrix(ncol = 2, nrow = nobservations+1)
#     states[1,] <- state
#     log_obs <- rep(0, nobservations)
#     for (t in 1:nobservations){
#       alpha <- rnorm(n = 1, mean = theta[1], sd = theta[2])
#       state <- pz_transition(state, alpha, t-1, theta[3:6])
#       states[t+1,] <- state
#       log_obs[t] <- rnorm(1, mean = log(state[1,1]), sd = 0.2)
#     }
#     return(list(states=states,log_obs=log_obs))
#   }
#
#
#   #prior
#   logprior <- function(transformed_theta){
#     ## evaluate prior density on the transformed parameter
#     # normal prior on theta1
#     density_eval <- dnorm(transformed_theta[1], mean = 0, sd = 10, log = TRUE)
#     # exponential prior on theta2
#     density_eval <- density_eval + dexp(exp(transformed_theta[2]), rate = 1, log = TRUE)
#     ## add Jacobian term to account for the transformation
#     # log transform for theta2
#     density_eval <- density_eval + transformed_theta[2]
#     # uniform prior on the other parameters
#     # logit transform
#     for (iparam in 3:6){
#       density_eval <- density_eval - transformed_theta[iparam] - 2 * log(1 + exp(-transformed_theta[iparam]))
#     }
#     return(density_eval)
#   }
#
#   #prior
#   logprior_naturalparams <- function(theta){
#     if(theta[2]<0 | any(theta[3:6]>1) | any(theta[3:6]<0)){
#       return(-Inf)
#     }else{
#       ## evaluate prior density on the transformed parameter
#       # normal prior on theta1
#       density_eval <- dnorm(theta[1], mean = 0, sd = 10, log = TRUE)
#       # exponential prior on theta2
#       density_eval <- density_eval + dexp(theta[2], rate = 1, log = TRUE)
#       # uniform prior on the other parameters
#       # logit transform
#       return(density_eval)
#     }
#   }
#
#
#
#   # posterior density function up to normalising constant
#   estimate_pf_target <- function(transformed_theta,nparticles){
#     log_prior<-logprior(transformed_theta)
#     if(log_prior==-Inf){
#       return(list(log_target=-Inf,path=NA))
#     }else{
#       pf_results <- pf(y, transformed_theta, pz_model,nparticles,ess_threshold = 0.75,systematic = TRUE)
#       log_target <- pf_results$loglik + log_prior
#       path <- pf_results$path
#       return(list(log_target=log_target,path=path))
#     }
#   }
#
#   #
#   #   rinit <- function(){
#   #
#   #     init_means <- c(0.44435692, -0.69450502, -0.17233041, -0.08403124,  0.01082659,  0.59587534)
#   #     init_vars <-  c(0.01726431, 0.01265800, 0.46045588, 0.86885887, 1.88446711, 1.77972311)
#   #     init_std <- 2*init_vars^0.5
#   #     return(init_means + init_std*rnorm(D_theta))
#   #   }
#
#
#
#   return(list(nobservations=nobservations,
#               theta_dgp=theta_dgp,
#               nmcmc=nmcmc,
#               generate_data=generate_data,
#               logprior=logprior,
#               logprior_naturalparams=logprior_naturalparams,
#               estimate_pf_target=estimate_pf_target,
#               nobservations_lowT=nobservations_lowT,
#               max_iterations=max_iterations,
#               nrep=nrep,
#               nrep_initial=nrep_initial,
#               data_file=data_file,
#               lowTcalibration_file=lowTcalibration_file,
#               N_arr=N_arr,
#               m=m))
#
#
#
# }
#
# sv_settings <- function(){
#
#   nobservations <- 1000
#   # untransformed parameters
#   theta_dgp <- c(0, 0, 0.5, 0.0625, 0.01)
#   D_theta <- length(theta_dgp)
#
#   # number particles for pmcmc
#   nparticles <- 20000
#
#   max_iterations <- 1000
#   nrep <- 100
#   nrep_initial <- 100
#
#
#   # pmmh settings
#   nmcmc <- 10000
#
#   cov_vec <- c(
#       0.0017669711, -0.0030540913, -0.0001128849, -0.001476659, -0.002448329,
#      -0.0030540913,  0.0070933927,  0.0005334601,  0.004451518,  0.003743204,
#      -0.0001128849,  0.0005334601,  0.2086905483,  0.493183121, -0.275014887,
#      -0.0014766586,  0.0044515177,  0.4931831211,  1.617976568, -1.117019053,
#      -0.0024483295,  0.0037432042, -0.2750148870, -1.117019053,  1.791992498
#   )
#   proposal_covariance <- 1*(2.2^2/D_theta)* structure(cov_vec, .Dim = c(5L, 5L))
#
#   # function to generate data
#   generate_data <- function(nobservations,transformed_theta){
#     x <- sv_model$rinitial(1, transformed_theta, NULL)
#     states <- matrix(nrow = nobservations+1, ncol = 2)
#     states[1,] <- x
#     y <- matrix(nrow = nobservations, ncol = 1)
#     for (time in 1:nobservations){
#       x <- sv_model$rtransition(x, transformed_theta, time, NULL, 1)
#       states[time+1,] <- x
#       y[time,] <- sv_model$robs(states[time+1,], transformed_theta)
#     }
#     return(list(states=states,y=y))
#   }
#
#
#   sv_hyperparams <- function(){
#     mu0mu <- 0
#     sigma02mu <- 10
#     mu0beta <- 0
#     sigma02beta <- 10
#     r0xi <- 1/5
#     r0w2 <- 1/5
#     r0lambda <- 1
#     return(list(mu0mu=mu0mu,
#                 sigma02mu=sigma02mu,
#                 mu0beta=mu0beta,
#                 sigma02beta=sigma02beta,
#                 r0xi=r0xi,
#                 r0w2=r0w2,
#                 r0lambda=r0lambda))
#   }
#
#   #prior
#   # Sampler from the prior distribution on parameters
#   # inputs: Ntheta (int)
#   # outputs: matrix (dimtheta by Ntheta) of prior draws
#   rprior = function(Ntheta){
#     hp <- sv_hyperparams()
#     mu0mu <- hp$mu0mu
#     sigma02mu <- hp$sigma02mu
#     mu0beta <- hp$mu0beta
#     sigma02beta <- hp$sigma02beta
#     r0xi <- hp$r0xi
#     r0w2 <- hp$r0w2
#     r0lambda <- hp$r0lambda
#
#     mu = rnorm(Ntheta, mu0mu, sqrt(sigma02mu))
#     beta = rnorm(Ntheta, mu0beta, sqrt(sigma02beta))
#     xi = rexp(Ntheta, r0xi)
#     w2 = rexp(Ntheta, r0w2)
#     lambda = rexp(Ntheta, r0lambda)
#     return (rbind(mu, beta, xi, w2, lambda))
#   }
#
#   # prior density on parameters
#   # inputs: theta (single vector), log (TRUE by default)
#   # outputs: prior (log)-density theta (double)
#   logprior <- function(theta_transformed){
#     hp <- sv_hyperparams()
#     mu0mu <- hp$mu0mu
#     sigma02mu <- hp$sigma02mu
#     mu0beta <- hp$mu0beta
#     sigma02beta <- hp$sigma02beta
#     r0xi <- hp$r0xi
#     r0w2 <- hp$r0w2
#     r0lambda <- hp$r0lambda
#
#     theta <- transformed2params(theta_transformed)
#     lmu <- dnorm(theta[1], mu0mu, sqrt(sigma02mu), log = TRUE)
#     lbeta <- dnorm(theta[2], mu0beta, sqrt(sigma02beta), log = TRUE)
#     lxi <- dexp(theta[3], r0xi, log = TRUE) + theta_transformed[3]
#     lw2 <- dexp(theta[4], r0w2, log = TRUE) + theta_transformed[4]
#     llambda <- dexp(theta[5], r0lambda, log = TRUE) + theta_transformed[5]
#     lp = lmu + lbeta + lxi + lw2 + llambda
#     return (lp)
#   }
#
#
#
#   # posterior density function up to normalising constant
#   estimate_pf_target <- function(transformed_theta,nparticles){
#     log_prior<-logprior(transformed_theta)
#     if(log_prior==-Inf){
#       return(list(log_target=-Inf,path=NA))
#     }else{
#       pf_results <- pf(y, transformed_theta, sv_model,nparticles,ess_threshold = 0.75,systematic = TRUE)
#       log_target <- pf_results$loglik + log_prior
#       path <- pf_results$path
#       return(list(log_target=log_target,path=path))
#     }
#   }
#
#
#   rinit <- function(){
#     mu0 <- c(0.0000000,  0.0000000, -0.6931472, -2.7725887, -4.6051702)
#     return(c(0.1,0.1,0.1,0.1,0.1)*rnorm(D_theta)+mu0)
#   }
#
#
#
#   return(list(nobservations=nobservations,
#               theta_dgp=theta_dgp,
#               nmcmc=nmcmc,
#               nparticles=nparticles,
#               proposal_covariance=proposal_covariance,
#               generate_data=generate_data,
#               logprior=logprior,
#               estimate_pf_target=estimate_pf_target,
#               rinit=rinit,
#               max_iterations=max_iterations,
#               nrep=nrep,
#               nrep_initial=nrep_initial))
#
#
#
# }
#
#
#
# # contains experiment settings for cpm
# correlated_lg_iid_test_model <- function(){
#   ##
#   ## model settings
#   sigma2 <- 1
#   var_prior <- 100
#   thetastar <- 0.5
#   nobs <- 1000
#   coupled_state <- FALSE
#
#   # PMCMC settings
#   N_block_arr <- rep(c(1,seq(5,40,5)),5)
#   N_pm_arr <- rep(seq(500,1200,100),5)
#   nmcmc_serial_block <- 250000
#   nmcmc_serial_pm <- 250000
#   # CPM settings
#   rho <- 0
#   rho_component <- 0.0
#   nrep <- 1000
#   nrep_initial <- 10000
#
#   mh_prop <- 0.01
#
#   rinit <- function() 2*runif(1)-1
#
#   data_fname <- 'inst/blockpm_toy/data/cpm_data.RData'
#   serial_fname <- 'inst/blockpm_toy/data/cpm_iid_serial.RData'
#   couple_fname <- 'inst/blockpm_toy/data/cpm_iid_coupled.RData'
#   differentT_fname <- 'inst/blockpm_toy/cpm_iid_differentT.RData'
#
#   return(list(sigma2=sigma2,
#               var_prior=var_prior,
#               thetastar=thetastar,
#               nobs=nobs,
#               coupled_state=coupled_state,
#               N_block_arr=N_block_arr,
#               N_pm_arr=N_pm_arr,
#               nmcmc_serial_block=nmcmc_serial_block,
#               nmcmc_serial_pm=nmcmc_serial_pm,
#               rho=rho,
#               rho_component=rho_component,
#               nrep=nrep,
#               nrep_initial=nrep_initial,
#               data_fname=data_fname,
#               rinit=rinit,
#               serial_fname=serial_fname,
#               couple_fname=couple_fname,
#               differentT_fname=differentT_fname,
#               mh_prop=mh_prop))
#
# }
#
#
#
#
#
# # contains experiment settings for LG
# coupled_like_est_lg_iid_model <- function(){
#   ##
#   ## model settings
#   sigma2 <- 1
#   var_prior <- 100^2
#   thetastar <- 0.5
#   nobs <- 200
#   mh_prop <- 0.07
#   nparticles <- 500
#
#   return(list(sigma2=sigma2,
#               var_prior=var_prior,
#               thetastar=thetastar,
#               nobs=nobs,
#               mh_prop=mh_prop,
#               nparticles=nparticles))
#
# }
#
#
#
#
#
# ising_model <- function(){
#
#   M <- 50
#   K <- 10000
#   max_iterations <- 1.1*K
#   mh_prop <- 1e-4*diag(1,2,2)#2e-5*matrix(c(5,-7,-7,12),2)
#   nmcmc <- 200000
#   nrep <- 100
#   nrep_initial <- 100
#
#   datafile <- 'inst/doubly_intractable/ising_data.RData'
#   serial_datafile <- 'inst/doubly_intractable/ising_serial.RData'
#   coupled_datafile <- 'inst/doubly_intractable/ising_coupled.RData'
#
#   J_true <- 0.3
#   H_true <- 0.1
#   theta_true <- c(J_true,H_true)
#
#   rinit <- function(){
#     return(c(runif(1),2*runif(1)-1))
#   }
#
#   return(list(
#     M=M,
#     K=K,
#     max_iterations=max_iterations,
#     mh_prop=mh_prop,
#     nmcmc=nmcmc,
#     nrep=nrep,
#     nrep_initial=nrep_initial,
#     J_true=J_true,
#     H_true=H_true,
#     datafile=datafile,
#     theta_true=theta_true,
#     serial_datafile=serial_datafile,
#     coupled_datafile=coupled_datafile,
#     rinit=rinit
#
#   ))
# }
#
#
# ising_model_no_external <- function(){
#
#   M <- 80
#   K <- 1100
#   max_iterations <- 1.1*K
#   mh_prop <- 1e-4*diag(1,1,1)#2e-5*matrix(c(5,-7,-7,12),2)
#   nmcmc <- 200000
#   nrep <- 1000
#   nrep_initial <- 1000
#   nrep_serial <- 10
#   external_field <- F
#
#   datafile <- 'inst/doubly_intractable/ising_data_noext.RData'
#   serial_datafile <- 'inst/doubly_intractable/ising_serial_noext.RData'
#   coupled_datafile <- 'inst/doubly_intractable/ising_coupled_noext.RData'
#
#   image_res_folder <- 'inst/doubly_intractable/noext_res/'
#
#   beta_critical <- 0.5*log(1+sqrt(2))
#
#   J_true <- beta_critical*(1/2)
#   theta_true <- J_true
#
#   rinit <- function(){
#     return(beta_critical*runif(1))
#   }
#
#   J_lim <- c(0,beta_critical)
#
#   return(list(
#     M=M,
#     K=K,
#     max_iterations=max_iterations,
#     mh_prop=mh_prop,
#     nmcmc=nmcmc,
#     nrep=nrep,
#     nrep_initial=nrep_initial,
#     nrep_serial=nrep_serial,
#     external_field=external_field,
#     J_true=J_true,
#     J_lim=J_lim,
#     datafile=datafile,
#     theta_true=theta_true,
#     serial_datafile=serial_datafile,
#     coupled_datafile=coupled_datafile,
#     image_res_folder=image_res_folder,
#     rinit=rinit
#
#   ))
# }
#
#
#
# random_effects_settings <- function(){
#
#
#   alpha_true <- 0.1814
#   beta_true <- c(0.4193,0.02289,0.1737,0.09443)
#   gamma_true <- c(0.1028,0.4033,0.1360,0.05086,0.2353)
#   delta_true <- -0.0288
#   lambda_true <- 2.1892
#   sigma_true <- 0.3522
#   sigma_u_true <- 0.3204
#   sigma_v_true <- 0.1463
#   sigma_w_true <- 0.04
#
#   pmcmc_nparticles <- 60 # seemed optimal for pmcmc
#   nparticles <- 10
#   coupled_state <- F
#
#   N_arr_block <- seq(1,10)
#   N_arr_pm <- c(1,seq(5,40,5))
#   N_arr_block_coupled <- c(1,5,10,20,30,40)
#   N_arr_pm_coupled <- c(1,5,seq(10,40,10))
#   nmcmc_serial <- 500000
#
#   max_iterations_block <- 200000
#   max_iterations_pm <- 100000
#   nrep <- 500
#
#   theta_true <- c(alpha_true,beta_true,gamma_true,delta_true,lambda_true,sigma_true,sigma_w_true)
#   D_theta <- length(theta_true)
#
#
#   data_fname <- 'inst/random_effects/banks.csv'
#   cov_fname <- 'inst/random_effects/prop_var.RData'
#   save_fname <- 'inst/random_effects/mt_save.RData'
#   save_serial_fname <- 'inst/random_effects/serial_save.RData'
#   save_mt_fname <- 'inst/random_effects/block_mt.RData'
#
#   return(list(data_fname=data_fname,
#               cov_fname=cov_fname,
#               save_fname=save_fname,
#               save_mt_fname=save_mt_fname,
#               theta_true=theta_true,
#               D_theta=D_theta,
#               nparticles=nparticles,
#               pmcmc_nparticles=pmcmc_nparticles,
#               coupled_state=coupled_state,
#               N_arr_block=N_arr_block,
#               N_arr_pm=N_arr_pm,
#               N_arr_block_coupled=N_arr_block_coupled,
#               N_arr_pm_coupled=N_arr_pm_coupled,
#               nmcmc_serial=nmcmc_serial,
#               max_iterations_block=max_iterations_block,
#               max_iterations_pm=max_iterations_pm,
#               save_serial_fname=save_serial_fname,
#               nrep=nrep))
#
# }
#
#
#
# dynamic_probit_settings <- function(){
#
#
#   cov_fname <- 'inst/dynamic_probit/prop_var.RData'
#
#
#   return(list(cov_fname=cov_fname))
#
# }
#
#
# election_model_settings <- function(){
#
#
#   cov_file <- 'inst/election_panel_data/cov_data.RData'
#   election_datafile <- 'inst/election_panel_data/elect_data.RData'
#   save_mt_fname <- 'inst/election_panel_data/data/bpm_pm_mt.RData'
#   save_Hkm_fname <- 'inst/election_panel_data/data/bpm_pm_Hkm.RData'
#   save_serial_fname <- 'inst/election_panel_data/data/bpm_serial.RData'
#   save_serial_PMMH_fname <- 'inst/election_panel_data/data/pm_PMMH_serial.RData'
#
#   nobs_selection <- 2000
#   N_cov_est <- 10*3
#   nmcmc_serial_cov_est <- 500000
#
#   # meeting time settings
#   N_arr_block_coupled <- c(5,10,20)*3
#   nrep_max <- 10000
#   max_iterations_block <- 50000
#
#   N_arr_pm_coupled <- c(600,700,800)*3#c(550,700,850)*3
#   max_iterations_pm <-50000
#
#   N_block_Hkm <- 10*3
#   N_pm_Hkm <- 700*3
#
#   nrep_Hkm <- 2*cores_requested
#   k <- 500
#   K <- 5000
#   K_large <- 20000
#
#   nmcmc_serial <- 250000
#   nmcmc_serial_PMMH <- 5000
#   nrep_serial <- 20
#   nrep_serial_PMMH <- 40
#
#   mean_est <- c(0.01464665, -0.08091029,  0.99063885,  0.95055937,  0.96655334 )
#   cov_est <- matrix(c(  6.964512e-05, -2.323131e-04,  1.292314e-06, -5.004271e-07, 3.558656e-07,
#                        -2.323131e-04,  1.432540e-03, -5.199375e-06,  4.784542e-06, 1.553754e-06,
#                         1.292314e-06, -5.199375e-06,  3.948509e-06,  4.837124e-06, 1.348633e-06,
#                        -5.004271e-07,  4.784542e-06,  4.837124e-06,  4.056879e-05, 2.256164e-05,
#                         3.558656e-07,  1.553754e-06,  1.348633e-06,  2.256164e-05, 2.352386e-05),nrow=length(mean_est))
#   return(list(cov_file=cov_file,
#               election_datafile=election_datafile,
#               save_serial_fname=save_serial_fname,
#               save_Hkm_fname=save_Hkm_fname,
#               nobs_selection=nobs_selection,
#               N_cov_est=N_cov_est,
#               nmcmc_serial_cov_est=nmcmc_serial_cov_est,
#               save_mt_fname=save_mt_fname,
#               N_arr_block_coupled=N_arr_block_coupled,
#               nrep_max=nrep_max,
#               max_iterations_block=max_iterations_block,
#               N_arr_pm_coupled=N_arr_pm_coupled,
#               max_iterations_pm=max_iterations_pm,
#               nrep_Hkm=nrep_Hkm,
#               nmcmc_serial_PMMH=nmcmc_serial_PMMH,
#               k=k,
#               K=K,
#               nmcmc_serial=nmcmc_serial,
#               nrep_serial=nrep_serial,
#               nrep_serial_PMMH=nrep_serial_PMMH,
#               mean_est=mean_est,
#               cov_est=cov_est,
#               N_block_Hkm=N_block_Hkm,
#               N_pm_Hkm=N_pm_Hkm))
#
#
# }
