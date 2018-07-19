# runs BPF and controlled SMC to see if costs are comparable

# remove all objects from R environment
rm(list = ls())
# load package
library(debiasedpmcmc)
### fix the random seed
set.seed(17)
setwd("")
source("binomialmodel.bpf.R")


# compare BFP cost with 4096 particles
theta_mle <- c(0.99, 0.11)
library(doParallel)
library(doRNG)
registerDoParallel(cores = 5)
res.bpf <- foreach(irep = 1:100, .combine = rbind) %dorng% {
  ptm <- proc.time()
  ll <- particle_filter(nparticles, binomialmodel, theta_mle, observations)$ll
  elapsed <- as.numeric((proc.time() - ptm)[3])
  c(ll, elapsed)
}

sd(res.bpf[,1])
mean(res.bpf[,2])

source("binomialmodel.R")

res.csmc <- foreach(irep = 1:100, .combine = rbind) %dorng% {
  ptm <- proc.time()
  ll <- csmc(nparticles, theta_mle, observations, niter)$ll
  elapsed <- as.numeric((proc.time() - ptm)[3])
  c(ll, elapsed)
}

sd(res.csmc[,1])
mean(res.csmc[,2])



