
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


