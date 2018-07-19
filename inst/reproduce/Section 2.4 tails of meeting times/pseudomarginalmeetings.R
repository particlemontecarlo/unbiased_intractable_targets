# This script can be used to
# reproduce figure 1 of the draft
# "Unbiased Markov chain Monte Carlo for intractable target distributions"
# by Middleton, Deligiannidis, Doucet & Jacob, 2018

# the figure illustrates the behavior of the meeting times,
# then the marginal algorithms transition from an exact MH algorithm
# to a pseudo-marginal MH algorithm, with increasing amounts of noise

# loack packages
library(debiasedpmcmc)
setmytheme()
library(doParallel)
registerDoParallel(cores = detectCores()-2)
library(doRNG)
library(dplyr)
library(ggthemes)
library(ggplot2)
rm(list = ls())
set.seed(1)
# set path
# setwd("~/Dropbox/UPMCMC/code/debiasedpmcmc/inst/pierre/")

### target distribution is bivariate Normal
dimension <- 2
post_mean <- 1:dimension
post_cov <- diag(1, dimension, dimension)

# covariance of the proposal
Sigma <- diag(1, dimension, dimension)
Sigma_chol <- chol(Sigma)
Sigma_chol_inv <- solve(chol(Sigma))
# sample from the target
rtarget <- function(n){
  return(fast_rmvnorm(n, post_mean, post_cov))
}
# evaluate target log pdf
dtarget <- function(theta){
  fast_dmvnorm(matrix(theta, nrow = 1), post_mean, post_cov)
}
# get initial distribution of chains, pseudo-marginal kernels (kern) and coupled kernels (ckern)
get_kernels <- function(noiselevel){
  dnoisytarget <- function(thetas){
    dtarget(thetas) + rnorm(1, mean = -0.5 * (noiselevel^2), sd = noiselevel)
  }
  # initial distribution of the chains
  rinit <- function(){
    valid <- FALSE
    while(!valid){
      chain_state <- rnorm(dimension, mean = 0, sd = 1)
      valid <- TRUE
    }
    current_target <- dnoisytarget(chain_state)
    return(list(chain_state = chain_state, current_target = current_target))
  }
  # kernel of pseudo-marginal MCMC
  kern <- function(chain_state, current_target){
    proposal <- chain_state + fast_rmvnorm(1, rep(0, dimension), post_cov)
    proposal_pdf <- dnoisytarget(proposal)
    if (log(runif(1)) < (proposal_pdf - current_target)){
      return(list(chain_state = proposal, current_target = proposal_pdf, accept = 1))
    } else {
      return(list(chain_state = chain_state, current_target = current_target, accept = 0))
    }
  }
  # coupled kernel
  ckern <- function(chain_state1, chain_state2, current_target1, current_target2){
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma_chol, Sigma_chol, Sigma_chol_inv, Sigma_chol_inv)
    proposal1 <- t(proposal_value[,1,drop=F])
    proposal2 <- t(proposal_value[,2,drop=F])
    logU <- log(runif(1))
    coupled_prop <- FALSE
    if (all(proposal1 == proposal2)){
      coupled_prop <- TRUE
      target_proposal1 <-  dnoisytarget(proposal1)
      target_proposal2 <- target_proposal1
    } else { ## proposals are different
      target_proposal1 <- dnoisytarget(proposal1)
      target_proposal2 <- dnoisytarget(proposal2)
    }
    accept <- 0
    if (logU < (target_proposal1 - current_target1)){
      chain_state1 <- proposal1
      current_target1 <- target_proposal1
      accept <- accept + 1
    }
    if (logU < (target_proposal2 - current_target2)){
      chain_state2 <- proposal2
      current_target2 <- target_proposal2
      accept <- accept + 1
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                current_target1 = current_target1, current_target2 = current_target2, coupled_prop = coupled_prop, accept = accept))
  }
  return(list(rinit = rinit, kern = kern, ckern = ckern))
}

# function to run coupled chains until max(meeting time, m)
# the "totalduration" argument refers to a number of seconds after which the execution is stopped
coupled_pmmh_ <- function(single_kernel, coupled_kernel, rinit, m = 1, max_iterations = Inf, preallocate = 10, verbose=FALSE, totalduration = Inf){
  ptm <- proc.time()
  initial_condition1 <- rinit()
  initial_condition2 <- rinit()
  chain_state1 <- initial_condition1$chain_state
  chain_state2 <- initial_condition2$chain_state
  log_pdf_state1 <- initial_condition1$current_target
  log_pdf_state2 <- initial_condition2$current_target
  p <- length(chain_state1)
  samples1 <- matrix(nrow = m+preallocate+1, ncol = p)
  nrowsamples1 <- m+preallocate+1
  samples2 <- matrix(nrow = m+preallocate, ncol = p)
  log_pdf1 <- matrix(nrow = m+preallocate+1, ncol = 1)
  log_pdf2 <- matrix(nrow = m+preallocate, ncol = 1)
  coupled_prop_arr <- matrix(nrow = m+preallocate, ncol = 1)

  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  log_pdf1[1,] <- log_pdf_state1
  log_pdf2[1,] <- log_pdf_state2

  current_nsamples1 <- 1
  iter <- 1
  mh_step <- single_kernel(chain_state1, log_pdf_state1)
  chain_state1 <- mh_step$chain_state
  log_pdf_state1 <- mh_step$current_target
  # Check if time is up already
  elapsedtime <- as.numeric((proc.time() - ptm)[3])
  if (elapsedtime > totalduration){
    # time's up
    return(list(finished = FALSE, message = "interrupted because time's up"))
  }

  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1
  log_pdf1[current_nsamples1,] <- log_pdf_state1

  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    elapsedtime <- as.numeric((proc.time() - ptm)[3])
    if (elapsedtime > totalduration){
      # time's up
      return(list(finished = FALSE, message = "interrupted because time's up"))
    }
    if(verbose){
      print(iter)
      print(elapsedtime)
    }
    iter <- iter + 1
    if (meet){
      mh_step <- single_kernel(chain_state1, log_pdf_state1)
      chain_state1 = mh_step$chain_state
      log_pdf_state1 = mh_step$current_target
      chain_state2 = mh_step$chain_state
      log_pdf_state2 = mh_step$current_target
    } else {
      mh_step <- coupled_kernel(chain_state1, chain_state2, log_pdf_state1, log_pdf_state2)
      chain_state1 <- mh_step$chain_state1
      chain_state2 <- mh_step$chain_state2
      log_pdf_state1 <- mh_step$current_target1
      log_pdf_state2 <- mh_step$current_target2
      coupled_prop <- mh_step$coupled_prop
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1-1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
      log_pdf1 <- rbind(log_pdf1, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf1)))
      log_pdf2 <- rbind(log_pdf2, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf2)))
      coupled_prop_arr <- rbind(coupled_prop_arr, matrix(NA, nrow = new_rows, ncol = ncol(coupled_prop_arr)))
    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    log_pdf1[current_nsamples1+1] <- log_pdf_state1
    log_pdf2[current_nsamples1] <- log_pdf_state2
    coupled_prop_arr[current_nsamples1] <- coupled_prop

    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  log_pdf1 <- log_pdf1[1:current_nsamples1,,drop=F]
  log_pdf2 <- log_pdf2[1:(current_nsamples1-1),,drop=F]
  coupled_prop_arr <- coupled_prop_arr[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2, log_pdf1=log_pdf1, log_pdf2=log_pdf2,
              meetingtime = meetingtime, iteration = iter, finished = finished, coupled_prop_arr=coupled_prop_arr))
}

## This creates meeting times (uncomment and run)
# nrep <- 10000
# noiselevels <- c(0, 0.5, 1, 1.5, 2)
# results.df <- data.frame()
# for (noiselevel in noiselevels){
#   print(noiselevel)
#   kernels <- get_kernels(noiselevel = noiselevel)
#   ccmeetings <- foreach(irep = 1:nrep, .combine = c) %dorng% {
#     res_ <- coupled_pmmh_(kernels$kern, kernels$ckern, kernels$rinit)
#     res_$meetingtime
#   }
#   ns <- 1:as.numeric(quantile(ccmeetings, probs = 0.99))
#   survprobs <- rep(0, length(ns))
#   for (i in seq_along(ns)){
#     survprobs[i] <- mean(ccmeetings > i)
#   }
#   df <- data.frame(ns = ns, survprobs = survprobs, noiselevel = noiselevel)
#   results.df <- rbind(results.df, df)
# }
## this saves meeting times, change the file path
# save(results.df, noiselevels, nrep, file = "survprobs.RData")

## this load meeting times (if the above code has been previously run, change the file path)
load(file = "survprobs.RData")

## Create plots as in the manuscript
g <- ggplot(results.df, aes(x = ns, y = survprobs, linetype = factor(noiselevel))) + geom_line() + scale_y_log10() +
  scale_linetype(name = "noise std dev") + xlab("n") + ylab("survival probability")
g
# ggsave(filename = "mvnormplusnoise.survprobs.pdf", plot = g, width = 7, height = 5)

## second plot
g <- ggplot(results.df %>% filter(ns > 30), aes(x = ns, y = survprobs, linetype = factor(noiselevel))) + geom_line() + scale_y_log10() +
  scale_linetype(name = "noise std dev") + xlab(expression(n)) + ylab("survival probability")
g <- g + scale_x_log10(breaks = c(30, 100, 200))
g
# ggsave(filename = "mvnormplusnoise.survprobs.loglog.pdf", plot = g, width = 7, height = 5)
