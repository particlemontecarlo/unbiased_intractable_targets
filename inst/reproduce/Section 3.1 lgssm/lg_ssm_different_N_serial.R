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


## experiment settings
settings <- lgssm_model_2params_100obs()

nobservations<- settings$nobservations
dimension<- settings$dimension
mu_0<- settings$mu_0
Sigma_0<- settings$Sigma_0
theta<- settings$theta
N_arr_serial<- settings$N_arr_serial
D_theta <- settings$D_theta
exp_name <- settings$exp_name
exp_name_serial <- settings$exp_name_serial
serial_data_file <- settings$serial_data_file
logprior <- settings$logprior
sigma_y <- settings$sigma_y
nmcmc <- settings$nmcmc
rinit <- settings$rinit

n_serial <- length(N_arr_serial)
run_mh_beforehand <- T



data_file <- settings$data_file
lowTcalibration_file <- settings$lowTcalibration_file

# iteration settings
cores_requested <- min(100,2*detectCores()-1)



# IO settings
exp_save_file <- "inst/exps/lg/optimal_N_coupling_different_N.RData"
exp_name <- "different_N"
run_on_server <- TRUE

# data file
load(data_file)
x <- as.matrix(x[1:(nobservations+1),],ncol=1)
y <- as.matrix(y[1:nobservations,],ncol=1)
print(sprintf('number of obs : %i',nobservations))




# set proposal covariance
mh_prop <- diag(c(0.04,0.04))#mh_prop_chosen


# request cores
registerDoMC(cores = cores_requested)

# model
ar_model <- get_lgssm_2params(dimension,sigma_y)

# memory for the particle filter
module_tree <<- Module("module_tree", PACKAGE = "debiasedpmcmc")
TreeClass <<- module_tree$Tree

# MH targets
mh_loglikelihood <- function(theta){
  kf_results <- kf(y, c(theta,sigma_y), mu_0, Sigma_0)
  return(kf_results$loglik)
}


# posterior density function up to normalising constant
mh_logtarget <- function(theta) mh_loglikelihood(theta) + logprior(theta)

# PF targets
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

# pm initialisation
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


# function to perform a run of PMMH, this will be called a few times in parallel
run_pmmh <- function(nparticles){

  pmmh_kernels <- get_pmmh_kernel(estimate_pftarget, Sigma_proposal = mh_prop, dimension = D_theta)
  pmmh_kernel <- pmmh_kernels$kernel
  coupled_pmmh_kernel <- pmmh_kernels$coupled_kernel

  chain_state <- rinit()
  pmcmc_chain <- matrix(0, nmcmc,D_theta)
  log_pdf_chain <- matrix(0,nmcmc)

  smc_res <- estimate_pftarget(chain_state,nparticles)
  log_pdf_state <- smc_res$log_target
  path_state <- smc_res$path

  accepts <- 0
  for (i in 1:nmcmc){

    print(sprintf('Progress : %.4f      AR : %.4f',i/nmcmc,accepts/i))
    pmmh_step <- pmmh_kernel(chain_state,log_pdf_state,path_state,i,nparticles)
    chain_state <- pmmh_step$chain_state
    log_pdf_state <-  pmmh_step$log_pdf_state
    path_state <- pmmh_step$path_state

    pmcmc_chain[i,] <- chain_state
    log_pdf_chain[i,] <- log_pdf_state

    if(any(pmcmc_chain[i,]!=pmcmc_chain[i-1,])){
      accepts <- accepts + 1
    }
  }

  return(list(pmcmc_chain=pmcmc_chain,log_pdf_chain=log_pdf_chain))
}


# perform long runs of PMMH
pmmh_runs <- foreach(i = 1:n_serial) %dopar% {
  nparticles <- N_arr_serial[i]
  pmmh_res <- run_pmmh(nparticles)
  return(pmmh_res)
}

# save the results
save(pmmh_runs,file=serial_data_file)

# load the results and plot some things of interest
load(serial_data_file)


library(coda)
acfs <- foreach(i = 1:n_serial) %dopar% {
    pmmh_run <- pmmh_runs[[i]]$pmcmc_chain
    nmcmc <- dim(pmmh_run)[1]
    return(spectrum0.ar(pmmh_run[10000:nmcmc,])$spec)
}

acf1 <- unlist(lapply(acfs,function(x) x[1]))
acf2 <- unlist(lapply(acfs,function(x) x[2]))

eff1 <- acf1*N_arr_serial


N_arr_unique <- unique(N_arr_serial)
eff_serial_mat <- matrix(NA,length(N_arr_unique),n_serial/length(N_arr_unique))
for(i in 1:length(N_arr_unique)){
    eff_serial_mat[i,] <- eff1[N_arr_serial==N_arr_unique[i]]
}
eff_unique <- rowMeans(eff_serial_mat)




eff1 <- acf1*N_arr_serial
name <- 'inst/lg_ssm/2params/serial_ineff.pdf'
pdf(name)
plot(N_arr_unique,eff_unique,ylab='Inefficiency',xlab='N',frame.plot=F)
grid(NULL,NULL)
dev.off()
print(sprintf('written to %s ',name))


