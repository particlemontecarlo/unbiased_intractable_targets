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
J_true <- settings$J_true
J_lim <- settings$J_lim
external_field <- settings$external_field
datafile <- settings$datafile
coupled_datafile <- settings$coupled_datafile
theta_true <- settings$theta_true
D_theta <- length(theta_true)
rinit <- settings$rinit


ising_model <- get_ising(M,J_lim,external_field=external_field,test_model=T)
log_prior <- ising_model$log_prior
d_log_lik <- ising_model$d_log_lik
r_lik <- ising_model$r_lik

run_initial_kernels <- F


# load data
load(datafile)
y <- data

# parallel settings
cores_requested <- min(detectCores()-1)
registerDoMC(cores = cores_requested)



Sigma_proposal <- mh_prop
kernels <- get_exchange_kernel(log_prior,d_log_lik,r_log_lik, Sigma_proposal, dimension=D_theta)
kernel <- kernels$kernel
coupled_kernel <- kernels$coupled_kernel
#
# if(run_initial_kernels){
#   # test single kernel
#
#   # run exchange algorithm
#   nmcmc <- 10000
#   chain <- matrix(NA,nmcmc,D_theta)
#   chain_state <- rinit()
#
#   chain[1,] <- chain_state
#   accepts <- 0
#
#   t1 <- Sys.time()
#   for(i in 2:nmcmc){
#     print(sprintf('iteration %i, accept rate %.4f',i,accepts/i))
#     chain[i,] <- kernel(chain[i-1,],i)$chain_state
#     if(any(chain[i,]!=chain[i-1,])){
#       accepts <- accepts+1
#     }
#
#   }
#   res_time <- Sys.time()-t1
#   hist(chain[1:i,1],probability=T,breaks=200)
#   abline(v=J_true,col='blue')
#   hist(chain[1:i,2],probability=T,breaks=200)
#   abline(v=H_true,col='blue')
#
#
#
#   # test coupled kernel
#
#
#   # run exchange algorithm
#   chain1 <- matrix(NA,nmcmc,D_theta)
#   chain2 <- matrix(NA,nmcmc,D_theta)
#   chain_state1 <- rinit()
#   chain_state2 <- rinit()
#
#   chain1[1,] <- chain_state1
#   chain2[1,] <- chain_state2
#
#   accepts1 <- 0
#   accepts2 <- 0
#   for(i in 2:nmcmc){
#     print(sprintf('iteration %i, accept rate 1 %.4f,  accept rate 1 %.4f',i,accepts1/i,accepts2/i))
#
#     mh_step <- coupled_kernel(chain_state1,chain_state2,i)
#
#     chain_state1 <- mh_step$chain_state1
#     chain_state2 <- mh_step$chain_state2
#
#     chain1[i,] <- chain_state1
#     chain2[i,] <- chain_state2
#
#     if(any(chain1[i,]!=chain1[i-1,])){
#       accepts1 <- accepts1+1
#     }
#
#     if(any(chain2[i,]!=chain2[i-1,])){
#       accepts2 <- accepts2+1
#     }
#   }
#
#   hist(chain1[,1],probability=T,breaks=200,main='chain1')
#   hist(chain2[,1],probability=T,breaks=200,main='chain2',add=T,col=rgb(0,0,1,0.1))
#
#   hist(chain1[,2],probability=T,breaks=200,main='chain1')
#   hist(chain2[,2],probability=T,breaks=200,main='chain2',add=T,col=rgb(0,0,1,0.1))
# }





# run coupling algorithm to look at meeting times
mh_batch_initial_timed <- foreach(i = 1:nrep_initial) %dopar% {
  t1 <- Sys.time()
  coupled_res <- coupled_mcmc_chains(kernel, coupled_kernel, rinit,K=1,max_iterations = max_iterations)
  res_time <- Sys.time()-t1
  return(list(coupled_res=coupled_res,res_time=res_time))
}

mh_batch_initial <- lapply(mh_batch_initial_timed,function(x) x$coupled_res)
batch_times_initial <- lapply(mh_batch_initial_timed,function(x) x$res_time)

mt <- sapply(mh_batch_initial,function(x) x$meetingtime)
mt_secs <- sapply(batch_times_initial,function(x) as.numeric(x,units='secs'))

name <- 'ising_mt.pdf'
pdf(name)
hist(mt,breaks=25,xlab='Meeting time',main='',cex=2)
dev.off()
print(sprintf('written to %s',name))

name <- 'ising_mt_secs.pdf'
pdf(name)
hist(mt_secs,breaks=25,xlab='Meeting time',main='',cex=2)
dev.off()
print(sprintf('written to %s',name))


# run coupling algorithm to look at meeting times
mh_batch_timed <- foreach(i = 1:nrep) %dopar% {
  t1 <- Sys.time()
  coupled_res <- coupled_mcmc_chains(kernel, coupled_kernel, rinit,K=K,max_iterations = max_iterations)
  res_time <- Sys.time()-t1
  return(list(coupled_res=coupled_res,res_time=res_time))
}

mh_batch <- lapply(mh_batch_timed,function(x) x$coupled_res)
batch_times <- lapply(mh_batch_timed,function(x) x$res_time)
batch_times_secs <- sapply(mh_batch_timed,function(x) as.numeric(x$res_time,units='secs'))


save(mh_batch,mh_batch_initial,batch_times,file=coupled_datafile)




# # get estimators from batch
# ests <- foreach(i = 1:nrep) %dopar% {
#   return(H_bar(mh_batch[[i]],k=0.1*K,K=K))
# }
#
# ests_mat <- do.call(rbind,ests)
#
# name <-  'ising_J.pdf'
# pdf(name,width=10,height=10)
# hist(ests_mat[,1])
# dev.off()
# print(sprintf('written to %s',name))
#
#
#
#
# h <- histogram_c_chains(mh_batch,1,k=0.1*K,K=K)
# g <- plot_histogram(h, with_bar = T)
# g <- g + labs(x = TeX('$J$'),y='Density')
#
# ggsave(filename = "unbiased_hist_ising.pdf", plot = g, width = 10, height = 10)
#
#
# name <-  'ising_mt.pdf'
# pdf(name,width=5,height=5)
# hist(mt,breaks=25,xlab='Meeting time',main='',cex=2)
# dev.off()
# print(sprintf('written to %s',name))
#
#
# name <-  'ising_img.pdf'
# pdf(name,width=5,height=5)
# image(y)
# dev.off()
# print(sprintf('written to %s',name))
#
#




