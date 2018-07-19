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


source("inst/reproduce/Section 4.2 coupledexchange/exp_settings.R")
# settings <- ising_model()
settings <- ising_model_no_external()

M <- settings$M
K <- settings$K
max_iterations <- settings$max_iterations
mh_prop <- settings$mh_prop
nmcmc <- settings$nmcmc
nrep <- settings$nrep
nrep_initial <- settings$nrep_initial
nrep_serial <- settings$nrep_serial
J_true <- settings$J_true
J_lim <- settings$J_lim
external_field <- settings$external_field
datafile <- settings$datafile
serial_datafile <- settings$serial_datafile
theta_true <- settings$theta_true
D_theta <- length(theta_true)
rinit <- settings$rinit


ising_model <- get_ising(M,J_lim,external_field=external_field,test_model=T)
log_prior <- ising_model$log_prior
d_log_lik <- ising_model$d_log_lik
r_lik <- ising_model$r_lik


# load data
load(datafile)
y <- data

# parallel settings
cores_requested <- min(100,detectCores()-1)
registerDoMC(cores = cores_requested)



Sigma_proposal <- mh_prop
kernels <- get_exchange_kernel(log_prior,d_log_lik,r_log_lik, Sigma_proposal, dimension=D_theta)
kernel <- kernels$kernel
coupled_kernel <- kernels$coupled_kernel

# test single kernel
run_exchange <- function(){

  # run exchange algorithm
  chain <- matrix(NA,nmcmc,D_theta)
  chain_state <- rinit()

  chain[1,] <- chain_state
  accepts <- 0

  t1 <- Sys.time()
  for(i in 2:nmcmc){
    print(sprintf('iteration %i, accept rate %.4f',i,accepts/i))
    chain[i,] <- kernel(chain[i-1,],i)$chain_state
    if(any(chain[i,]!=chain[i-1,])){
      accepts <- accepts+1
    }

  }
  res_time <- Sys.time()-t1

  return(list(chain=chain,res_time=res_time))
}


# run coupling algorithm to look at meeting times
mh_serial_batches <- foreach(i = 1:nrep_serial) %dopar% {
  run_exchange()
}



save(mh_serial_batches,file=serial_datafile)
print(sprintf('written to %s',serial_datafile))


# plot results
load(serial_datafile)


acfs <- sapply(mh_serial_batches,function(x) spectrum0.ar(x$chain)$spec)



chain <- mh_serial_batches[[1]]$chain
hist(chain[1:i,1],probability=T,breaks=200)
abline(v=J_true,col='blue')
hist(chain[1:i,2],probability=T,breaks=200)
abline(v=H_true,col='blue')

acf(chain[,1],300)

library(coda)
serial_ineff <- spectrum0.ar(chain)$spec







