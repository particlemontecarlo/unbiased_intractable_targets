if (!require(debiasedpmcmc)){
  library(devtools)
  devtools::document()
}
library(debiasedpmcmc)
library(coda)
library(doMC)
library(MASS)
library(ggplot2)
library(latex2exp)
rm(list = ls())
set.seed(17)
setmytheme()

source("inst/reproduce/Section 3.1 lgssm/exp_settings.R")


##
## Try a possible coupling with correlated latent variables
##
## Settings
settings <- lgssm_model_2params_100obs()

nobservations<- settings$nobservations
dimension<- settings$dimension
mu_0<- settings$mu_0
Sigma_0<- settings$Sigma_0
theta<- settings$theta
N_arr_mt<- settings$N_arr_mt
D_theta <- settings$D_theta
exp_name <- settings$exp_name
exp_name_serial <- settings$exp_name_serial
logprior <- settings$logprior
sigma_y <- settings$sigma_y
rinit <- settings$rinit
mt_data_file <- settings$mt_data_file

run_mh_beforehand <- T

n_N_arr <- length(N_arr_mt)

data_file <- settings$data_file
lowTcalibration_file <- settings$lowTcalibration_file

# iteration settings
cores_requested <- min(100,detectCores()-1)

# coupling settings
nrep <- settings$nrep_mt
nrep_pmmh <- settings$nrep_mt

# specifying the maximum iterations is helpful for diagnostic purposes
# N.B. it should be checked that no chains have exceed this.
max_iterations <- 20000


load(data_file)
x <- as.matrix(x[1:(nobservations+1),],ncol=1)
y <- as.matrix(y[1:nobservations,],ncol=1)
print(sprintf('number of obs : %i',nobservations))




# set proposal covariance
mh_prop <- diag(c(0.04,0.04))#mh_prop_chosen


# request cores
registerDoMC(cores = cores_requested)

ar_model <- get_lgssm_2params(dimension,sigma_y)

module_tree <<- Module("module_tree", PACKAGE = "debiasedpmcmc")
TreeClass <<- module_tree$Tree


mh_loglikelihood <- function(theta){
  kf_results <- kf(y, c(theta,sigma_y), mu_0, Sigma_0)
  return(kf_results$loglik)
}


# posterior density function up to normalising constant
mh_logtarget <- function(theta) mh_loglikelihood(theta) + logprior(theta)

# posterior density function up to normalising constant (we use an unnormalised prior)
estimate_pftarget <- function(theta,nparticles){
  log_prior<-logprior(theta)
  if(log_prior==-Inf){
    return(list(log_target=-Inf,path=NA))
  }else{
    pf_results <- pf(y, theta, ar_model,nparticles)
    log_target <- pf_results$loglik + log_prior
    path <- pf_results$path
    return(list(log_target=log_target,path=path))
  }
}

pmmh_init <- function(nparticles){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  pf_target1 <- estimate_pftarget(chain_state1,nparticles)
  pf_target2 <- estimate_pftarget(chain_state2,nparticles)
  log_pdf_state1 <- pf_target1$log_target
  log_pdf_state2 <- pf_target2$log_target
  path_state1 <- pf_target1$path
  path_state2 <- pf_target2$path

  return(list(chain_state1=chain_state1,chain_state2=chain_state2,
              log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2,
              path_state1=path_state1,path_state2=path_state2))
}



# get MH and PMMH kernels
mh_kernels <- get_mh_kernel(logtarget = mh_logtarget, Sigma_proposal = mh_prop, dimension = D_theta)
single_mh_kernel <- mh_kernels$kernel
coupled_mh_kernel <- mh_kernels$coupled_kernel

pmmh_kernels <- get_pmmh_kernel(estimate_pftarget, Sigma_proposal = mh_prop, dimension = D_theta)
single_pmmh_kernel <- pmmh_kernels$kernel
coupled_pmmh_kernel <- pmmh_kernels$coupled_kernel

# run mh before hand
mh_mt_batch <- foreach(i = 1:nrep) %dopar% {
  coupled_mcmc_chains(single_mh_kernel, coupled_mh_kernel, rinit)
}

mh_mts <- matrix(sapply(mh_mt_batch,function(x) x$meetingtime),nrow=1)
# hist(mh_mts)



# run coupled pmmh
pmmh_mts <- matrix(NA,n_N_arr,nrep_pmmh)
pmmh_times <-list()
pmmh_mt_batches <- list()
for(i in 1:n_N_arr){
  print(sprintf('pm batch %i of %i',i,n_N_arr))

  nparticles <- N_arr_mt[i]
  pmmh_batch_timed <- foreach(i = 1:nrep_pmmh) %dopar% {
    t1<-Sys.time()
    couple_res <- coupled_pmmh_chains(single_pmmh_kernel, coupled_pmmh_kernel, pmmh_init,nparticles,nobs,max_iterations=max_iterations)
    couple_time <- Sys.time()-t1
    return(list(couple_res=couple_res,couple_time=couple_time))
  }

  pmmh_batch <- lapply(pmmh_batch_timed,function(x) x$couple_res)
  pmmh_batch_times <- lapply(pmmh_batch_timed,function(x) x$couple_time)

  pmmh_mts[i,] <- sapply(pmmh_batch,function(x) x$meetingtime)
  pmmh_times[[i]] <- pmmh_batch_times
  pmmh_mt_batches[[i]] <- pmmh_batch
}


# save mt results
save(pmmh_mt_batches,mh_mt_batch,file=mt_data_file)
print(sprintf('data written to %s',mt_data_file))




