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
set.seed(17)
setmytheme()

source("inst/reproduce/Section 3.1 lgssm/exp_settings.R")


##
## Settings
settings <- lgssm_model_2params_scalingT()

nobs_arr<- settings$nobs_arr
nobservations_lowT <- settings$nobservations_lowT
dimension<- settings$dimension
mu_0<- settings$mu_0
Sigma_0<- settings$Sigma_0
theta<- settings$theta
N_arr<- settings$N_arr
D_theta <- settings$D_theta
exp_name <- settings$exp_name
exp_name_serial <- settings$exp_name_serial
logprior <- settings$logprior
sigma_y <- settings$sigma_y
mt_resfile <- settings$mt_resfile
cov_0_prop <- settings$cov_0_prop
cov_0_init <- settings$cov_0_init

n_N_arr <- length(N_arr)
n_nobs_arr <-  length(nobs_arr)



# iteration settings
cores_requested <- min(100,detectCores()-1)

# coupling settings - this specifies the number of runs for both MH and PM
nrep <- 1000
nrep_mh <- 10000
max_iterations <- 10000



# generate a long run of the data that we will use sequentially with increasing number of particles
nobs_gendata <- 2*max(nobs_arr)
x_all <- matrix(0, nrow = nobs_gendata+1, ncol = dimension)
y_all <- matrix(0, nrow = nobs_gendata, ncol = dimension)
x_all[1,] <- fast_rmvnorm(1, mu_0, Sigma_0)
for (t in 1:nobs_gendata){
  x_all[t+1,] <- x_all[t,,drop=F] %*% diag(theta[1], dimension, dimension)+ fast_rmvnorm_chol(1, rep(0, dimension), diag(theta[2], dimension, dimension))
  y_all[t,]  <- x_all[t+1,] + fast_rmvnorm_chol(1, rep(0, dimension), diag(sigma_y, dimension, dimension))
}



# request cores
registerDoMC(cores = cores_requested)

# model
ar_model <- get_lgssm_2params(dimension,sigma_y)

# used for particle filter
module_tree <<- Module("module_tree", PACKAGE = "debiasedpmcmc")
TreeClass <<- module_tree$Tree


#### define targets for MH
# loglikelihood for mcmc
mh_loglikelihood <- function(theta){
  kf_results <- kf(y, c(theta,sigma_y), mu_0, Sigma_0)
  return(kf_results$loglik)
}


# posterior density function up to normalising constant
mh_logtarget <- function(theta) mh_loglikelihood(theta) + logprior(theta)

#### definte targets for pmmh
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

# initialisation for pseudo marginal
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




###### get meeting times under the different scalings
# perform scaling 1
# set init
rinit <- function(){
  # sample according to the multivariate normal, truncated so the initial conditions
  # are positive for variance parameters
  theta_init <- -Inf
  while(logprior(theta_init)==-Inf){
    theta_init <- mvrnorm(1,theta,Sigma_n_init)#5*Sigma_n)
  }
  return(theta_init)
}

# loop through each number of observations. need to set:
# 1) the initial distribution to pi_0=N(mustar,5*Sigma_n)
# 2) the propsal covariance to Sigma_n
# 3) make sure only n observations are being used
# where Sigma_n = n_0/n*Sigma_est
t12 <- Sys.time()

# mh results
mh_times_scaling1 <- matrix(NA,n_nobs_arr,nrep_mh)
mh_mts_scaling1 <- matrix(NA,n_nobs_arr,nrep_mh)
mh_ar1_scaling1 <- matrix(NA,n_nobs_arr,nrep_mh)
mh_ar2_scaling1 <- matrix(NA,n_nobs_arr,nrep_mh)


for(i in 1:n_nobs_arr){

  nobs <- nobs_arr[i]

  Sigma_n_prop <- cov_0_prop/nobs
  Sigma_n_init <- cov_0_init/nobs
  mh_prop <- Sigma_n_prop

  y <- as.matrix(y_all[1:nobs,],ncol=1)

  # get MH and PMMH kernels
  mh_kernels <- get_mh_kernel(logtarget = mh_logtarget, Sigma_proposal = mh_prop, dimension = D_theta)
  single_mh_kernel <- mh_kernels$kernel
  coupled_mh_kernel <- mh_kernels$coupled_kernel

  print(sprintf('Running MH with nobs=%i, nrep_mh=%i',nobs,nrep_mh))
  mh_t1 <- Sys.time()
  mh_batch_timed <- foreach(i = 1:nrep_mh) %dopar% {
    mh_t1_i <- Sys.time()
    mh_res <- coupled_mcmc_chains(single_mh_kernel, coupled_mh_kernel, rinit)
    return(list(mh_res=mh_res,mh_res_time=Sys.time()-mh_t1_i))
  }
  mh_time <- Sys.time() - mh_t1
  mh_batch <- lapply(mh_batch_timed,function(x) x$mh_res)

  mh_times_scaling1[i,] <- sapply(mh_batch_timed,function(x) as.numeric(x$mh_res_time,unit='secs'))
  mh_mts_scaling1[i,] <- sapply(mh_batch,function(x) x$meetingtime)
  mh_ar1_scaling1[i,] <- sapply(mh_batch,function(x) acc_rate(x$samples1))
  mh_ar2_scaling1[i,] <- sapply(mh_batch,function(x) acc_rate(x$samples2))



}


# pmmh results
pmmh_times_scaling1 <- matrix(NA,n_nobs_arr,nrep)
pmmh_mts_scaling1 <- matrix(NA,n_nobs_arr,nrep)
pmmh_ar1_scaling1 <- matrix(NA,n_nobs_arr,nrep)
pmmh_ar2_scaling1 <- matrix(NA,n_nobs_arr,nrep)


for(i in 1:n_nobs_arr){

  nobs <- nobs_arr[i]

  Sigma_n_prop <- cov_0_prop/nobs
  Sigma_n_init <- cov_0_init/nobs
  mh_prop <- Sigma_n_prop

  y <- as.matrix(y_all[1:nobs,],ncol=1)

  pmmh_kernels <- get_pmmh_kernel(estimate_pftarget, Sigma_proposal = mh_prop, dimension = D_theta)
  single_pmmh_kernel <- pmmh_kernels$kernel
  coupled_pmmh_kernel <- pmmh_kernels$coupled_kernel

  nparticles <- N_arr[i]
  print(sprintf('Running PMMH with nobs=%i, nparticles=%i, nrep=%i',nobs,nparticles,nrep))
  pmmh_t1 <- Sys.time()
  pmmh_batch_timed <- foreach(i = 1:nrep) %dopar% {
    pmmh_t1_i <- Sys.time()
    pmmh_res<-coupled_pmmh_chains(single_pmmh_kernel, coupled_pmmh_kernel, pmmh_init,nparticles,nobs,max_iterations=max_iterations)
    return(list(pmmh_res=pmmh_res,pmmh_res_time=Sys.time()-pmmh_t1_i))
  }
  pmmh_time <- Sys.time() - pmmh_t1
  print('PMMH time taken:')
  print(pmmh_time)

  pmmh_batch <- lapply(pmmh_batch_timed,function(x) x$pmmh_res)

  pmmh_times_scaling1[i,] <- sapply(pmmh_batch_timed,function(x) as.numeric(x$pmmh_res_time,unit='secs'))

  pmmh_mts_scaling1[i,] <- sapply(pmmh_batch,function(x) x$meetingtime)
  pmmh_ar1_scaling1[i,] <- sapply(pmmh_batch,function(x) acc_rate(x$samples1))
  pmmh_ar2_scaling1[i,] <- sapply(pmmh_batch,function(x) acc_rate(x$samples2))


}
t22 <- Sys.time()

print('Total time taken:')
print(t22-t12)


save.image(file=mt_resfile)
print(sprintf('saved to %s',mt_resfile))



#### perform second scalings
rinit <- function() {
  return(c(runif(1),5*runif(1)))
}

# loop through each number of observations. need to set:
# 1) the initial distribution to pi_0=N(mustar,5*Sigma_n)
# 2) the propsal covariance to Sigma_n
# 3) make sure only n observations are being used
# where Sigma_n = n_0/n*Sigma_est
t1 <- Sys.time()

# mh results
mh_times_scaling2 <- matrix(NA,n_nobs_arr,nrep_mh)
mh_mts_scaling2 <- matrix(NA,n_nobs_arr,nrep_mh)
mh_ar1_scaling2 <- matrix(NA,n_nobs_arr,nrep_mh)
mh_ar2_scaling2 <- matrix(NA,n_nobs_arr,nrep_mh)


for(i in 1:n_nobs_arr){

  nobs <- nobs_arr[i]

  Sigma_n_prop <- cov_0_prop/nobs
  Sigma_n_init <- cov_0_init/nobs
  mh_prop <- Sigma_n_prop

  y <- as.matrix(y_all[1:nobs,],ncol=1)

  # get MH and PMMH kernels
  mh_kernels <- get_mh_kernel(logtarget = mh_logtarget, Sigma_proposal = mh_prop, dimension = D_theta)
  single_mh_kernel <- mh_kernels$kernel
  coupled_mh_kernel <- mh_kernels$coupled_kernel

  print(sprintf('Running MH with nobs=%i, nrep_mh=%i',nobs,nrep_mh))
  mh_t1 <- Sys.time()
  mh_batch_timed <- foreach(i = 1:nrep_mh) %dopar% {
    mh_t1_i <- Sys.time()
    mh_res <- coupled_mcmc_chains(single_mh_kernel, coupled_mh_kernel, rinit)
    return(list(mh_res=mh_res,mh_res_time=Sys.time()-mh_t1_i))
  }
  mh_time <- Sys.time() - mh_t1
  mh_batch <- lapply(mh_batch_timed,function(x) x$mh_res)

  mh_times_scaling2[i,] <- sapply(mh_batch_timed,function(x) as.numeric(x$mh_res_time,unit='secs'))
  mh_mts_scaling2[i,] <- sapply(mh_batch,function(x) x$meetingtime)
  mh_ar1_scaling2[i,] <- sapply(mh_batch,function(x) acc_rate(x$samples1))
  mh_ar2_scaling2[i,] <- sapply(mh_batch,function(x) acc_rate(x$samples2))
}

# pmmh results
pmmh_times_scaling2 <- matrix(NA,n_nobs_arr,nrep)
pmmh_mts_scaling2 <- matrix(NA,n_nobs_arr,nrep)
pmmh_ar1_scaling2 <- matrix(NA,n_nobs_arr,nrep)
pmmh_ar2_scaling2 <- matrix(NA,n_nobs_arr,nrep)
for(i in 1:n_nobs_arr){

  nobs <- nobs_arr[i]

  Sigma_n_prop <- cov_0_prop/nobs
  Sigma_n_init <- cov_0_init/nobs
  mh_prop <- Sigma_n_prop

  y <- as.matrix(y_all[1:nobs,],ncol=1)

  pmmh_kernels <- get_pmmh_kernel(estimate_pftarget, Sigma_proposal = mh_prop, dimension = D_theta)
  single_pmmh_kernel <- pmmh_kernels$kernel
  coupled_pmmh_kernel <- pmmh_kernels$coupled_kernel

  nparticles <- N_arr[i]
  print(sprintf('Running PMMH with nobs=%i, nparticles=%i, nrep=%i',nobs,nparticles,nrep))
  pmmh_t1 <- Sys.time()
  pmmh_batch_timed <- foreach(i = 1:nrep) %dopar% {
    pmmh_t1_i <- Sys.time()
    pmmh_res<-coupled_pmmh_chains(single_pmmh_kernel, coupled_pmmh_kernel, pmmh_init,nparticles,nobs,max_iterations=max_iterations)
    return(list(pmmh_res=pmmh_res,pmmh_res_time=Sys.time()-pmmh_t1_i))
  }
  pmmh_time <- Sys.time() - pmmh_t1
  print('PMMH time taken:')
  print(pmmh_time)

  pmmh_batch <- lapply(pmmh_batch_timed,function(x) x$pmmh_res)

  pmmh_times_scaling2[i,] <- sapply(pmmh_batch_timed,function(x) as.numeric(x$pmmh_res_time,unit='secs'))

  pmmh_mts_scaling2[i,] <- sapply(pmmh_batch,function(x) x$meetingtime)
  pmmh_ar1_scaling2[i,] <- sapply(pmmh_batch,function(x) acc_rate(x$samples1))
  pmmh_ar2_scaling2[i,] <- sapply(pmmh_batch,function(x) acc_rate(x$samples2))


}
t2 <- Sys.time()

print('Total time taken:')
print(t2-t1)


save.image(file=mt_resfile)
print(sprintf('saved to %s',mt_resfile))

#### Plot results
# the following gets quantiles of the meeting times for each value of T and
# each scaling then plots then produces the results

# plot graphs
load(mt_resfile)

# specify the quantiles to plot.
q1 <- 0.8
q4 <- 0.99



q_arr <- c(q1,q4)

q_res_mh1 <- matrix(NA,n_nobs_arr,length(q_arr))
for(i in 1:n_nobs_arr){
  q_res_mh1[i,] <- as.vector(quantile(mh_mts_scaling1[i,],q_arr))
}

q_res_pmmh1 <- matrix(NA,n_nobs_arr,length(q_arr))
for(i in 1:n_nobs_arr){
  q_res_pmmh1[i,] <- as.vector(quantile(pmmh_mts_scaling1[i,],q_arr))
}


q_res_mh2 <- matrix(NA,n_nobs_arr,length(q_arr))
for(i in 1:n_nobs_arr){
  q_res_mh2[i,] <- as.vector(quantile(mh_mts_scaling2[i,],q_arr))
}

q_res_pmmh2 <- matrix(NA,n_nobs_arr,length(q_arr))
for(i in 1:n_nobs_arr){
  q_res_pmmh2[i,] <- as.vector(quantile(pmmh_mts_scaling2[i,],q_arr))
}

nobs_min <- 100



ggplot_mt_compare <- function(q_res_mh,q_res_pmmh,nobs_min){
  # function used to plot the meeting times using ggplot
  nobs_mask <- nobs_arr>=nobs_min
  x <- nobs_arr[nobs_mask]

  # mh
  y_mh <- q_res_mh[nobs_mask,]
  y_pmmh <- q_res_pmmh[nobs_mask,]

  y <- cbind(y_mh,y_pmmh)


  label_mh <- rep(NA,length(q_arr))
  for(i in 1:length(q_arr)){ label_mh[i] <- sprintf('%s (MH)',prettyNum(q_arr[i])) }
  label_pmmh <- rep(NA,length(q_arr))
  for(i in 1:length(q_arr)){ label_pmmh[i] <- sprintf('%s (PMMH)',prettyNum(q_arr[i])) }

  colnames(y) <-c(label_mh,label_pmmh) #q_arr#label_mh

  ym <- reshape2::melt(y)
  ym$x<-x[ym$Var1]


  g <- ggplot( ym,aes(x=x,y=value))+
    geom_line(aes(linetype=Var2)) +
    labs(x='T',y='meeting time') +
    theme(legend.title=element_blank()) +
    guides(linetype = guide_legend(nrow = 2))
  g

}


setmytheme()

name <- 'pmmh_scalingT_scaling1.pdf'
ggplot_mt_compare(q_res_mh1,q_res_pmmh1,nobs_min)
ggsave(name,width=7,height=5)
print(sprintf('written to %s',name))

name <- 'pmmh_scalingT_scaling2.pdf'
ggplot_mt_compare(q_res_mh2,q_res_pmmh2,nobs_min)
ggsave(name,width=7,height=5)
print(sprintf('written to %s',name))

