if (!require(debiasedpmcmc)){
  library(devtools)
  devtools::document()
}
library(debiasedpmcmc)
library(coda)
library(doMC)
library(MASS)
library(latex2exp)
library(ggplot2)
rm(list = ls())
set.seed(19)
setmytheme()

source("inst/reproduce/Section 3.1 lgssm/exp_settings.R")



##
## Settings
settings <- lgssm_model_2params_100obs()



nobservations<- settings$nobservations
dimension<- settings$dimension
mu_0<- settings$mu_0
Sigma_0<- settings$Sigma_0
theta<- settings$theta
N_arr_Hkm<- settings$N_arr_Hkm
D_theta <- settings$D_theta
exp_name <- settings$exp_name
exp_name_serial <- settings$exp_name_serial
logprior <- settings$logprior
sigma_y <- settings$sigma_y
rinit <- settings$rinit
n_rep_Hkm <- settings$n_rep_Hkm
coupledres_folder <- settings$coupledres_folder


run_mh_beforehand <- F

n_N_arr_Hkm <- length(N_arr_Hkm)

data_file <- settings$data_file
lowTcalibration_file <- settings$lowTcalibration_file

# iteration settings
cores_requested <- min(100,floor((detectCores()-1)))

# coupling settings
n_rep_Hkm <- 10000
K <- 5000
max_iterations <- K*1.1



load(data_file)
x <- as.matrix(x[1:(nobservations+1),],ncol=1)
y <- as.matrix(y[1:nobservations,],ncol=1)
print(sprintf('number of obs : %i',nobservations))

# load(lowTcalibration_file)



# set proposal covariance
# mh_prop <- diag(c(0.03,0.03))#mh_prop_chosen
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


# function to help strip out the meeting times
get_finished_meeting_times <- function(batch_res){
  finished_meeting_times <- sapply(batch_res, function(x){ if(x$finished){return(x$meetingtime)}})
  null_mask <- sapply(finished_meeting_times, is.null)
  dnf <- sum(null_mask)
  print(paste(dnf,' failed to terminate'))
  finished_meeting_times<-unlist(finished_meeting_times[!null_mask ])
  return(list(finished_meeting_times=finished_meeting_times,dnf=dnf))
}


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


# get estimators for metropolis hastings
mh_batch <- foreach(i = 1:n_rep_Hkm) %dopar% {
  coupled_mcmc_chains(single_mh_kernel, coupled_mh_kernel, rinit,K=K,max_iterations = max_iterations)
}

save(mh_batch,file=sprintf('%slgssm_Hkm_mh_batch.RData',coupledres_folder))

# get estimators for coupled PMMH
for(i in 1:n_N_arr_Hkm){
  print(sprintf('pm batch %i of %i',i,n_N_arr_Hkm))

  nparticles <- N_arr_Hkm[i]
  timed_pmmh_batch <- foreach(i = 1:n_rep_Hkm) %dopar% {
    t1 <- Sys.time()
    coupling_res <- coupled_pmmh_chains(single_pmmh_kernel, coupled_pmmh_kernel, pmmh_init,nparticles,nobs,K=K,max_iterations=max_iterations,verbose=F)
    exec_time <- Sys.time()-t1

    return(list(coupling_res=coupling_res,exec_time=exec_time))
  }
  pmmh_batch <- lapply(timed_pmmh_batch,function(x) x$coupling_res)
  batch_timings <- lapply(timed_pmmh_batch,function(x) x$exec_time)

  save(pmmh_batch,batch_timings,settings,mh_prop,file=sprintf('%slgssm_Hkm_%i.RData',coupledres_folder,i))
}


# this is just used to save everything - it can be very large watch out!
save.image(sprintf('%slgssm_Hkm_allsave.RData',coupledres_folder))
# save.image(sprintf('inst/lg_ssm/2params/lgssm_Hkm_allsave.RData'))



