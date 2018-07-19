## generate long run of PMMH with bootstrap particle filter,
## takes a few days to complete (~ 1 second per iteration)

# remove all objects from R environment
rm(list = ls())
# load package
library(debiasedpmcmc)
### fix the random seed
set.seed(17)
setwd("")
source("binomialmodel.bpf.R")

current <- rinit()
## number of MCMC steps
nmcmc <- 250000
# save file
savef <- paste0("pmmh.bpf.T", nmcmc, ".N", nparticles, ".RData")

## store chain in matrix
chain <- matrix(nrow = nmcmc, ncol = 2)
target_pdfs <- rep(0, nmcmc)
## number of accepts
naccepts <- 0
ptm <- proc.time()
elapsed <- 0
## loop over MCMC iteration
for (imcmc in 1:nmcmc){
  if (imcmc %% 100 ==1){
    print(imcmc)
    cat("accept rate:", 100*naccepts/imcmc, "%\n")
    elapsed <- as.numeric((proc.time() - ptm)[3])
    print(elapsed)
    save(nmcmc, nparticles, chain, target_pdfs, naccepts, imcmc, current, elapsed, file = savef)
  }
  current <- single_kernel(current$chain_state, current$current_target)
  chain[imcmc,] <- current$chain_state
  target_pdfs[imcmc] <- current$current_target
  naccepts <- naccepts + current$accept
}
save(nmcmc, nparticles, chain, target_pdfs, naccepts, current, elapsed, file = savef)

#
load(savef)
naccepts / imcmc * 100
matplot(chain[1000:imcmc,1], type = "l")
matplot(chain[1000:imcmc,2], type = "l")

matplot(chain[1000:nmcmc,1], type = "l")
matplot(chain[1000:nmcmc,2], type = "l")

hist(chain[1000:imcmc,1])
hist(chain[1000:imcmc,2])

# print(naccepts/nmcmc)
# matplot(chain[100:nmcmc,1], type = "l")
# matplot(chain[100:nmcmc,2], type = "l")
# hist(chain[300:nmcmc,1])
# hist(chain[300:nmcmc,2])

