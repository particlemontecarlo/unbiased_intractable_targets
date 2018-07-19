### load data as in controlled SMC paper
re <- read.csv("thaldata.csv", header = FALSE, col.names = FALSE)
observations <- matrix(re[1:3000,1], ncol = 1)
datalength <- nrow(observations)
### plot data
# g <- qplot(x = 1:datalength, y = observations, geom = "line") + xlab("time") + ylab("response")
# g

### logistic function
logistic <- function(x) 1/(1+exp(-x))

### Original model
## x_0 ~ Normal(0,1)
## x_t | x_{t-1} ~ Normal(alpha x_{t-1}, sigma^2)
## y_t | x_t ~ Binomial(50, logistic(x_t))
### comments
## the latent process is AR(1) with parameters theta = (alpha, sigma^2)
## given the latent process, the logistic function maps to (0,1)
## and the observation y_t are Binomial with 50 trials and probability logistic(x_t)

binomialmodel <- list(
  rinit = function(nparticles, theta, ...){
    return(matrix(rnorm(nparticles), ncol = nparticles))
  },
  rtransition = function(xparticles, theta, time, rand, precomputed, ...){
    return(matrix(theta[1] * xparticles + rnorm(nparticles, mean = 0, sd = sqrt(theta[2])), nrow = 1))
  },
  dmeasurement = function(xparticles, theta, observation, precomputed, ...){
    return(dbinom(observation, 50, prob = logistic(xparticles[1,]), log = TRUE))
  },
  precompute = function(theta){
    return(NULL)
  },
  dimension = 1)

particle_filter <- function(nparticles, model, theta, observations){
  datalength <- nrow(observations)
  precomputed <- try(model$precompute(theta))
  if (inherits(precomputed, "try-error")){
    precomputed <- NULL
  }
  # initialization
  xparticles <- model$rinit(nparticles, theta)
  normweights <- rep(1/nparticles, nparticles)
  ll <- 0
  #
  # xmeans <- matrix(nrow = datalength + 1, ncol = nrow(xparticles))
  # xmeans[1,] <- apply(xparticles, 1, function(v) sum(v * normweights))
  # step t > 1
  for (time in 1:datalength){
    ancestors <- systematic_resampling_n(normweights, nparticles, runif(1))
    xparticles <- xparticles[,ancestors]
    if (is.null(dim(xparticles))) xparticles <- matrix(xparticles, nrow = model$dimension)
    xparticles <- model$rtransition(xparticles, theta, time, precomputed)
    logw <- model$dmeasurement(xparticles, theta, observations[time,], precomputed)
    maxlw <- max(logw)
    w <- exp(logw - maxlw)
    # update log likelihood estimate
    ll <- ll + maxlw + log(mean(w))
    normweights <- w / sum(w)
    #
    # xmeans[time+1,] <- apply(xparticles, 1, function(v) sum(v * normweights))
  }
  return(list(ll = ll))
}

## prior on theta
## For alpha in (0,1) we put a Uniform prior,
## and for sigma2 we put an inverse gamma prior with parameters (1, 0.2)

## log pdf of inverse Gamma distribution
digamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

## log pdf of prior distribution on theta
dprior <- function(theta){
  if (theta[1] < 0 || theta[1] > 1){
    return(-Inf)
  } else {
    if (theta[2] < 0){
      return(-Inf)
    } else {
      return(digamma(theta[2], 1, 0.1))
    }
  }
}
## with this many particles
nparticles <- 2^12


rinit <- function(){
  theta <- c(runif(1, min = 0, max = 1), runif(1, min = 0, max = 1)) # uniform on [0,1] for both parameters
  target_pdf <- particle_filter(nparticles, binomialmodel, theta, observations)$ll
  return(list(chain_state = theta, current_target = target_pdf))
}

sd_proposal <- 5 * c(0.002, 0.01)
Sigma_chol <- diag(sd_proposal)
Sigma_chol_inv <- diag(1/sd_proposal)

single_kernel <- function(chain_state, current_target){
  ## proposal
  proposal <- chain_state + rnorm(2, mean = 0, sd = sd_proposal)
  ## evaluate prior log pdf
  target_proposal <- dprior(proposal)
  if (is.finite(target_proposal)){
    ## only run particle filter if parameter is in the acceptable range
    target_proposal <- target_proposal + particle_filter(nparticles, binomialmodel, proposal, observations)$ll
  }
  ## accept or not
  if (log(runif(1)) < (target_proposal - current_target)){
    return(list(chain_state = proposal, current_target = target_proposal, accept = TRUE))
  } else {
    return(list(chain_state = chain_state, current_target = current_target, accept = FALSE))
  }
}

coupled_kernel <- function(chain_state1, chain_state2, current_target1, current_target2){
  proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                     Sigma_chol, Sigma_chol, Sigma_chol_inv, Sigma_chol_inv)
  proposal1 <- proposal_value[,1]
  proposal2 <- proposal_value[,2]
  logU <- log(runif(1))
  coupled_prop <- FALSE
  if (all(proposal1 == proposal2)){
    coupled_prop <- TRUE
    target_proposal1 <- dprior(proposal1)
    if (is.finite(target_proposal1)){
      target_proposal1 <- target_proposal1 + particle_filter(nparticles, binomialmodel, proposal1, observations)$ll
    }
    target_proposal2 <- target_proposal1
  } else { ## proposals are different
    target_proposal1 <- dprior(proposal1)
    if (is.finite(target_proposal1)){
      target_proposal1 <- target_proposal1 + particle_filter(nparticles, binomialmodel, proposal1, observations)$ll
    }
    target_proposal2 <- dprior(proposal2)
    if (is.finite(target_proposal2)){
      target_proposal2 <- target_proposal2 + particle_filter(nparticles, binomialmodel, proposal2, observations)$ll
    }
  }
  if (logU < (target_proposal1 - current_target1)){
    chain_state1 <- proposal1
    current_target1 <- target_proposal1
  }
  if (logU < (target_proposal2 - current_target2)){
    chain_state2 <- proposal2
    current_target2 <- target_proposal2
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
              current_target1 = current_target1, current_target2 = current_target2, coupled_prop = coupled_prop))
}

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
