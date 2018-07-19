## Example of file to be run in batch mode on the cluster
## This one runs coupled PMMH with BPF for 23 hours (soft deadline), with m = 0 (to produce meeting times)

# time allocated to this
TIME <- 23*3600
# m
m <- 0

# change file paths!
library(debiasedpmcmc)
source("binomialmodel.bpf.R")

library(parallel)
# parallel RNG using L'Ecuyer et al (2002)
RNGkind("L'Ecuyer-CMRG") # L'Ecuyer CMRG required for multiple streams
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

ptm_irep <- proc.time()
last_start <- ptm_irep
nsamples_ <- 0
durations_ <- c()
elapsedtime <- 0
timeleft <- TIME
cchains_ <- list()
while(elapsedtime < TIME){
  if (nsamples_ == 0){ # then produce a sample
    res <- try(coupled_pmmh_(single_kernel, coupled_kernel, rinit, m = m, max_iterations = Inf, verbose = T))
    if (inherits(res, "try-error")) res <- list(finished = FALSE)
  } else { # then produce a sample if time permits
    res <- try(coupled_pmmh_(single_kernel, coupled_kernel, rinit, m = m, max_iterations = Inf, totalduration = timeleft))
    if (inherits(res, "try-error")) res <- list(finished = FALSE)
  }
  elapsedtime <- as.numeric((proc.time() - ptm_irep)[3])
  timeleft <- TIME - elapsedtime
  duration_ <- as.numeric((proc.time() - last_start)[3])
  last_start <- proc.time()
  if (res$finished){
    nsamples_ <- nsamples_ + 1
    cchains_[[nsamples_]] <- res
    durations_ <- c(durations_, duration_)
  }
  save(cchains_, nsamples_, durations_, file = paste0("output/neuro.bpf.N", nparticles, ".m", m, ".ID", igrid, ".RData"))
}

