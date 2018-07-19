if (!require(debiasedpmcmc)){
  library(devtools)
  devtools::document()
}
library(debiasedpmcmc)
library(coda)
library(ggplot2)
library(doMC)
library(MASS)
library(foreign)
library(TruncatedNormal)
rm(list = ls())
set.seed(1)
setmytheme()




cores_requested <- min(100,detectCores()-1)
registerDoMC(cores = cores_requested)



source("inst/reproduce/Section 4.1 blockpm/exp_settings.R")
settings <- election_model_settings()
# settings
cov_file <- settings$cov_file
election_datafile <- settings$election_datafile
save_mt_fname <- settings$save_mt_fname
save_serial_fname <- settings$save_serial_fname
nobs_selection <- settings$nobs_selection
N_block_Hkm <- settings$N_block_Hkm

mean_est <- settings$mean_est
cov_est <- settings$cov_est

nmcmc_serial <- settings$nmcmc_serial
nrep_serial <- settings$nrep_serial


load(file=election_datafile)

X_save <- X_save[1:nobs_selection,1:3]
y_save <- y_save[1:nobs_selection,1:3]


X <- matrix(c(t(X_save)))

y <- c(t(y_save))
y <- 1*(y==1)


T_length <- 3
nobs <- length(y)/T_length

# parameter ests
theta <- c(-0.01266113,  0,0.95,0.95,0.95)
D_theta <- length(theta)



y_obs <- matrix(NA,ncol=T_length, nrow=nobs)
for(i in 1:T_length){
  y_obs[,i] <- y[seq(i,T_length*nobs,T_length)]
}



init_cov <- 0.01^2*diag(1,D_theta)


coupled_state <- F


re_model <- get_mv_probit_model(mean_est,init_cov)
block_pm_loglikelihood <- re_model$block_pm_loglikelihood
block_pm_logtarget <- re_model$block_pm_logtarget
state_crn_sample <- re_model$state_crn_sample
latent_state <- re_model$latent_state
logprior <- re_model$logprior
rinit <- re_model$rinit
block_coupled_init <- re_model$block_coupled_init
pm_logtarget <- re_model$pm_logtarget
pm_coupled_init <- re_model$pm_coupled_init



mh_prop <- (2.38^2/D_theta)*cov_est
D_theta <- length(theta)
# run block pm



t1 <- Sys.time()



# compare the distribution of meeting times
couple_block_res <- list()

N <- N_block_Hkm

print(sprintf('getting serial block pm, %i particles',N))

# marginal kernel with joint updates only
cpm_kernels <- get_cpm_T_blocking(logtarget = block_pm_logtarget, Sigma_proposal=mh_prop, dimension=D_theta,
                                  rho_component = 0,joint_update=F,rho = 0,
                                  single_kernel_only = F)
single_cpm_kernel_component <- cpm_kernels$kernel
coupled_cpm_kernel_component <- cpm_kernels$coupled_kernel


run_bpm <- function(){
  initial_conditions <- block_coupled_init( N,nobs,coupled_state=coupled_state)

  chain <- matrix(NA,nmcmc_serial,D_theta)

  t1 <- Sys.time()

  chain_state1 <- initial_conditions$chain_state1
  log_pdf_state1 <- initial_conditions$log_pdf_state1
  state_crn1 <- initial_conditions$state_crn1
  loglik_t1 <- initial_conditions$loglik_t1
  accepts <- 0
  for(iter in 1:nmcmc_serial){
    if(iter%%1000==0){
      print(sprintf('progress %.4f, accepts=%.4f',iter/nmcmc_serial,accepts/iter))
    }

    chain[iter,] <- chain_state1

    cpm_step <- single_cpm_kernel_component(chain_state1,state_crn1,log_pdf_state1,loglik_t1,iter)

    if(any(chain_state1!=cpm_step$chain_state)){
      accepts <- accepts+1
    }

    chain_state1 <- cpm_step$chain_state
    log_pdf_state1 <- cpm_step$log_pdf_state
    state_crn1 <- cpm_step$state_crn
    loglik_t1 <- cpm_step$loglik_t

  }
  res_time <- Sys.time()-t1

  chain <- chain[1:iter,]
  n_chain <- dim(chain)[1]

  return(list(chain=chain,res_time=res_time))
}



# run a few block pseudo-marginal in parallel and estimate the inefficiency
bpm_batch <- foreach(pp = 1:nrep_serial) %dopar% {

  bpm_res <- run_bpm()
  return(bpm_res)
}

# save.image(file=save_serial_fname)


print('Serial estimation complete!')

load(file=save_serial_fname)



bpm_chains <- lapply(bpm_batch,function(x) x$chain)
bpm_burnin_chains <- lapply(bpm_chains,function(x) x[(0.1*nmcmc_serial):nmcmc_serial,])
bpm_burnin_chains_correct_frequency <- lapply(bpm_burnin_chains,function(x) x[seq(1,dim(x)[1],2),])
bpm_acf_est <- sapply(bpm_burnin_chains_correct_frequency,function(x) spectrum0.ar(rowSums(x)+rowSums(x^2))$spec) * (nmcmc_serial/(nmcmc_serial-0.1*nmcmc_serial))




