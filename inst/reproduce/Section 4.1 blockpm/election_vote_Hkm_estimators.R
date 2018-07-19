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
settings <- election_model_settings()
cov_file <- settings$cov_file
election_datafile <- settings$election_datafile
save_Hkm_fname <- settings$save_Hkm_fname
nobs_selection <- settings$nobs_selection


N_block_Hkm <- settings$N_block_Hkm
N_pm_Hkm <- settings$N_pm_Hkm

mean_est <- settings$mean_est
cov_est <- settings$cov_est
k <- settings$k
K <- settings$K
K_large <- settings$K_large

nrep_Hkm <- settings$nrep_Hkm

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


#


init_cov <- 0.01^2*diag(1,D_theta)

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
coupled_state <- F




# compare the distribution of meeting times
N <- N_block_Hkm

print(sprintf('getting Hkm block pm, %i particles',N))

# marginal kernel with joint updates only
cpm_kernels <- get_cpm_T_blocking(logtarget = block_pm_logtarget, Sigma_proposal=mh_prop, dimension=D_theta,
                                  rho_component = 0,joint_update=F,rho = 0,
                                  single_kernel_only = F)
single_cpm_kernel_component <- cpm_kernels$kernel
coupled_cpm_kernel_component <- cpm_kernels$coupled_kernel


couple_block <- foreach(pp = 1:nrep_Hkm) %dopar% {
  ran_experiment <- F
  t1 <- Sys.time()

  couple_res <- tryCatch({
    res <- coupled_cpm_chains_Tblocking(single_cpm_kernel_component, coupled_cpm_kernel_component, block_coupled_init,N,nobs,iters_per_cycle=1,K = 2*K, max_iterations = 4*K,verbose=F)
  }, error = function(err) {
    return(err)
  })

  couple_res_time <- Sys.time() - t1
  ran_experiment <- T
  return(list(couple_res=couple_res,couple_res_time=couple_res_time,ran_experiment=ran_experiment))
}

save.image(file=save_Hkm_fname)
print('Hkm estimators complete!')
print(sprintf('saved to %s',save_Hkm_fname))




# estimate for pm


# compare the distribution of meeting times
N <- N_pm_Hkm

print(sprintf('getting Hkm for standard pm, %i particles',N))

# run coupled chain for pm
pmmh_kernels <- get_pm_kernel(logtarget = pm_logtarget,state_crn_sample, latent_state, Sigma_proposal = mh_prop, dimension = D_theta)
single_pmmh_kernel <- pmmh_kernels$kernel
coupled_pmmh_kernel <- pmmh_kernels$coupled_kernel


couple_pm <- foreach(pp = 1:nrep_Hkm) %dopar% {
  ran_experiment <- F
  t1 <- Sys.time()

  couple_res <- tryCatch({
    res <- coupled_pm_chains(single_pmmh_kernel, coupled_pmmh_kernel, pm_coupled_init,N,nobs,K=K,max_iterations=2*K,coupled_state = coupled_state)
  }, error = function(err) {
    return(err)
  })

  ran_experiment <- T
  couple_res_time <- Sys.time() - t1

  return(list(couple_res=couple_res,couple_res_time=couple_res_time,ran_experiment=ran_experiment))
}

save.image(file=save_Hkm_fname)
print('Hkm estimators pm complete!')
print(sprintf('saved to %s',save_Hkm_fname))



# perform analysis



load(file=save_Hkm_fname)


k <- 500
K <- 5000

estimators_mt <- sapply(couple_block,function(x) x$couple_res$meetingtime)
estimators <- sapply(couple_block,function(x) H_bar(x$couple_res,h=function(x) sum(x)+sum(x^2),k=2*k,K=2*K))
block_timings <- sapply(couple_block,function(x) as.numeric(x$couple_res_time,units='secs'))
estimate_mean <- mean(estimators)

estimate_expected_cost <- mean(2*(estimators_mt/2) + max(1,K+1-(estimators_mt/2)))
estimate_ineff <- (var((estimators)))*(estimate_expected_cost)

estimate_sd_err <- (var((estimators)))^0.5/sqrt(nrep_Hkm)



pm_estimators_mt <- sapply(couple_pm,function(x) x$couple_res$meetingtime)
pm_estimators <- sapply(couple_pm,function(x) H_bar(x$couple_res,h=function(x) sum(x)+sum(x^2),k=k,K=K))

pm_estimate_mean <- mean(pm_estimators)
pm_estimate_expected_cost <- mean(2*(pm_estimators_mt) + max(1,K+1-(pm_estimators_mt)))
pm_estimate_ineff <- (var((pm_estimators)))*pm_estimate_expected_cost
pm_acf_est <- sapply(couple_pm,function(x) spectrum0.ar(rowSums(x$couple_res$samples1[k:K,]))$spec)
pm_timings <- sapply(couple_pm,function(x) as.numeric(x$couple_res_time,units='secs'))

(pm_estimate_ineff*N_pm_Hkm)/(2*estimate_ineff*N_block_Hkm)






