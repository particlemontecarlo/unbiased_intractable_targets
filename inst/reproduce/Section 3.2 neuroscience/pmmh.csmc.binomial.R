## generate long run of PMMH with controlled SMC,
## with a larger std deviation than used in the paper by Heng et al.
## takes a few days to complete (~ 1 second per iteration)

# remove all objects from R environment
rm(list = ls())
# load package
library(debiasedpmcmc)
### fix the random seed
set.seed(17)
setwd("")
source("binomialmodel.R")

current <- rinit()
current$chain_state <- current$theta
current$current_target <- current$target_pdf
## number of MCMC steps
nmcmc <- 250000
# save file
savef <- paste0("pmmh.csmc.T", nmcmc, ".N", nparticles, ".niter", niter, ".RData")

### run and store chain in matrix
# chain <- matrix(nrow = nmcmc, ncol = 2)
# target_pdfs <- rep(0, nmcmc)
# ## number of accepts
# naccepts <- 0
# ptm <- proc.time()
# elapsed <- 0
# ## loop over MCMC iteration
# for (imcmc in 1:nmcmc){
#   if (imcmc %% 100 ==1){
#     print(imcmc)
#     cat("accept rate:", 100*naccepts/imcmc, "%\n")
#     elapsed <- as.numeric((proc.time() - ptm)[3])
#     print(elapsed)
#     save(nmcmc, nparticles, chain, target_pdfs, naccepts, imcmc, current, elapsed, file = savef)
#   }
#   current <- single_kernel(current$chain_state, current$current_target)
#   chain[imcmc,] <- current$chain_state
#   target_pdfs[imcmc] <- current$current_target
#   naccepts <- naccepts + current$accept
# }
# save(nmcmc, nparticles, chain, target_pdfs, naccepts, current, elapsed, file = savef)

load(savef)
naccepts / nmcmc * 100
## Compute inefficiency
burnin <- 10000
testchain <- apply(chain[burnin:(imcmc-1),], MARGIN = 1, function(x) sum(x) + sum(x^2))
mean(testchain)
library(coda)
var(testchain)
# inefficiency parametrized by time
spectrum0(testchain)$spec

