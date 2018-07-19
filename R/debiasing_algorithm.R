# This file consists of the general-purpose functions coupled_chains, continue_coupled_chains,
# and H_bar, which implement our debiased MCMC algorithm for general kernels and functions h(.)

# from coupled_chains -----------------------------------------------------
# Run coupled chains until max(tau, K) where tau is the meeting time and K specified by user
#'@rdname coupled_pmmh_chains
#'@title Coupled PMCMC chains
#'@description Sample two PMMH chains, each following \code{single_kernel} marginally,
#' and \code{coupled_kernel} jointly, until min(max(tau, K), max_iterations), where tau
#' is the first time at which the two chains meet (i.e. take the same value exactly).
#' Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
#' are initialized from the distribution provided in \code{rinit}.
#'
#'  See \code{\link{get_hmc_kernel}}
#' for an example of function returning the appropriate kernels.
#'
#'@param single_kernel function taking a state (in a vector) and an iteration, and returning
#' a list with a key named \code{chain_state} and containing the next state.
#'@param coupled_kernel function taking two states (in two vectors) and an iteration,
#'and returning a list with keys \code{chain_state1} and \code{chain_state2}.
#'@param rinit function taking no arguments are returning an initial state for a Markov chain.
#'@param K number of iterations desired (will be proportional to the computing cost if meeting occurs before \code{K},
#' default to 1).
#'@param max_iterations number of iterations at which the function stops if it is still running  (default to Inf).
#'@param preallocate  expected number of iterations, used to pre-allocate memory (default to 10).
#'@export
coupled_pmmh_chains <- function(single_kernel, coupled_kernel, pmmh_init,N,nobs,K = 1, max_iterations = Inf, preallocate = 10,verbose=FALSE){

  initial_conditions <- pmmh_init(N)
  chain_state1 <- initial_conditions$chain_state1
  chain_state2 <- initial_conditions$chain_state2
  log_pdf_state1 <- initial_conditions$log_pdf_state1
  log_pdf_state2 <- initial_conditions$log_pdf_state2
  path_state1 <- initial_conditions$path_state1
  path_state2 <- initial_conditions$path_state2

  p <- length(chain_state1)
  samples1 <- matrix(nrow = K+preallocate+1, ncol = p)
  nrowsamples1 <- K+preallocate+1
  samples2 <- matrix(nrow = K+preallocate, ncol = p)
  log_pdf1 <- matrix(nrow = K+preallocate+1, ncol = 1)
  log_pdf2 <- matrix(nrow = K+preallocate, ncol = 1)
  coupled_prop_arr <- matrix(nrow = K+preallocate, ncol = 1)

  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  log_pdf1[1,] <- log_pdf_state1
  log_pdf2[1,] <- log_pdf_state2


  current_nsamples1 <- 1
  iter <- 1
  mh_step <- single_kernel(chain_state1,log_pdf_state1,path_state1, iter,N)
  chain_state1 <- mh_step$chain_state
  log_pdf_state1 <- mh_step$log_pdf_state
  path_state1 <- mh_step$path_state

  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1
  log_pdf1[current_nsamples1,] <- log_pdf_state1

  accepts1 <- 0
  accepts2 <- 0
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){

    iter <- iter + 1
    if (meet){
      mh_step <- single_kernel(chain_state1,log_pdf_state1,path_state1, iter,N)
      chain_state1 = mh_step$chain_state
      log_pdf_state1 = mh_step$log_pdf_state
      chain_state2 = mh_step$chain_state
      log_pdf_state2 = mh_step$log_pdf_state
      path_state1 <- mh_step$path_state
      path_state2 <- mh_step$path_state
    } else {
      mh_step <- coupled_kernel(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2, path_state1, path_state2,iter,N)
      chain_state1 <- mh_step$chain_state1
      chain_state2 <- mh_step$chain_state2
      log_pdf_state1 <- mh_step$log_pdf_state1
      log_pdf_state2 <- mh_step$log_pdf_state2
      path_state1 <- mh_step$path_state1
      path_state2 <- mh_step$path_state2

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


    if(any(samples1[current_nsamples1+1,]!=samples1[current_nsamples1,])){
      accepts1 <- accepts1 + 1
    }
    if(any(samples2[current_nsamples1,]!=samples2[current_nsamples1-1,])){
      accepts2 <- accepts2 + 1
    }


    if(verbose){
      print(sprintf('iter : %i   Progress : %.4f   AR_1 : %.4f  AR_2 : %.4f',iter,iter/max_iterations,accepts1/iter,accepts2/iter))
    }

    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(K, tau) steps
    if (iter >= max(meetingtime, K)){
      finished <- TRUE
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  log_pdf1 <- log_pdf1[1:current_nsamples1,,drop=F]
  log_pdf2 <- log_pdf2[1:(current_nsamples1-1),,drop=F]
  coupled_prop_arr <- coupled_prop_arr[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2, log_pdf1=log_pdf1, log_pdf2=log_pdf2,
              meetingtime = meetingtime, iteration = iter, finished = finished,coupled_prop_arr=coupled_prop_arr,
              accepts1=accepts1,accepts2=accepts2))
}






# Run coupled chains until max(tau, K) where tau is the meeting time and K specified by user
#'@rdname coupled_pm_chains
#'@title Coupled pseudo-marginal chains
#'@description Sample two pseudo-marginal MCMC chains, each following \code{single_kernel} marginally,
#' and \code{coupled_kernel} jointly, until min(max(tau, K), max_iterations), where tau
#' is the first time at which the two chains meet (i.e. take the same value exactly).
#' Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
#' are initialized from the distribution provided in \code{rinit}.
#'
#'
#'@param single_kernel function taking a state (in a vector) and an iteration, and returning
#' a list with a key named \code{chain_state} and containing the next state.
#'@param coupled_kernel function taking two states (in two vectors) and an iteration,
#'and returning a list with keys \code{chain_state1} and \code{chain_state2}.
#'@param rinit function taking no arguments are returning an initial state for a Markov chain.
#'@param K number of iterations desired (will be proportional to the computing cost if meeting occurs before \code{K},
#' default to 1).
#'@param max_iterations number of iterations at which the function stops if it is still running  (default to Inf).
#'@param preallocate  expected number of iterations, used to pre-allocate memory (default to 10).
#'@export
coupled_pm_chains <- function(single_kernel, coupled_kernel, coupled_init,N,nobs,K = 1, max_iterations = Inf, preallocate = 10,coupled_state=TRUE,finish_time=NULL){

  early_stopping=F

  initial_conditions <- coupled_init(N,nobs,coupled_state=coupled_state)
  chain_state1 <- initial_conditions$chain_state1
  chain_state2 <- initial_conditions$chain_state2
  log_pdf_state1 <- initial_conditions$log_pdf_state1
  log_pdf_state2 <- initial_conditions$log_pdf_state2


  p <- length(chain_state1)
  samples1 <- matrix(nrow = K+preallocate+1, ncol = p)
  nrowsamples1 <- K+preallocate+1
  samples2 <- matrix(nrow = K+preallocate, ncol = p)
  log_pdf1 <- matrix(nrow = K+preallocate+1, ncol = p)
  log_pdf2 <- matrix(nrow = K+preallocate, ncol = p)

  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  log_pdf1[1,] <- log_pdf_state1
  log_pdf2[1,] <- log_pdf_state2


  current_nsamples1 <- 1
  iter <- 1
  mh_step <- single_kernel(chain_state1,log_pdf_state1, iter,N,nobs)
  chain_state1 = mh_step$chain_state
  log_pdf_state1 = mh_step$log_pdf_state

  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1
  log_pdf1[current_nsamples1,] <- log_pdf_state1

  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations && !early_stopping){
    iter <- iter + 1
    if (meet){
      mh_step <- single_kernel(chain_state1,log_pdf_state1, iter,N,nobs)
      chain_state1 = mh_step$chain_state
      log_pdf_state1 = mh_step$log_pdf_state
      chain_state2 = mh_step$chain_state
      log_pdf_state2 = mh_step$log_pdf_state
    } else {
      mh_step <- coupled_kernel(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2, iter,N,nobs,coupled_state=coupled_state)
      chain_state1 <- mh_step$chain_state1
      chain_state2 <- mh_step$chain_state2
      log_pdf_state1 <- mh_step$log_pdf_state1
      log_pdf_state2 <- mh_step$log_pdf_state2

      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1 - 1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
      log_pdf1 <- rbind(log_pdf1, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf1)))
      log_pdf2 <- rbind(log_pdf2, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf2)))
    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    log_pdf1[current_nsamples1+1,] <- log_pdf_state1
    log_pdf2[current_nsamples1,] <- log_pdf_state2

    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(K, tau) steps
    if (iter >= max(meetingtime, K)){
      finished <- TRUE
    }
    if(!is.null(finish_time)){
      if(Sys.time()>finish_time){
        early_stopping <- T
      }
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  log_pdf1 <- log_pdf1[1:current_nsamples1,,drop=F]
  log_pdf2 <- log_pdf2[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2, log_pdf1=log_pdf1, log_pdf2=log_pdf2,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}





# from coupled_chains -----------------------------------------------------
# Run coupled chains until max(tau, K) where tau is the meeting time and K specified by user
#'@rdname coupled_cpm_chains
#'@title Coupled correlated pseudo-marginal chains
#'@description Sample two CPM chains, each following \code{single_kernel} marginally,
#' and \code{coupled_kernel} jointly, until min(max(tau, K), max_iterations), where tau
#' is the first time at which the two chains meet (i.e. take the same value exactly).
#' Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
#' are initialized from the distribution provided in \code{rinit}.
#'
#'
#'@param single_kernel function taking a state (in a vector) and an iteration, and returning
#' a list with a key named \code{chain_state} and containing the next state.
#'@param coupled_kernel function taking two states (in two vectors) and an iteration,
#'and returning a list with keys \code{chain_state1} and \code{chain_state2}.
#'@param rinit function taking no arguments are returning an initial state for a Markov chain.
#'@param K number of iterations desired (will be proportional to the computing cost if meeting occurs before \code{K},
#' default to 1).
#'@param max_iterations number of iterations at which the function stops if it is still running  (default to Inf).
#'@param preallocate  expected number of iterations, used to pre-allocate memory (default to 10).
#'@export
coupled_cpm_chains <- function(single_kernel, coupled_kernel, coupled_init,nparticles,nobs,iters_per_cycle=1,K = 1, max_iterations = Inf, preallocate = 10,coupled_state=F){
  if((max_iterations%%iters_per_cycle)!=0){
    stop('Error max_iterations must be a multiple of iters_per_cycle')
  }

  initial_conditions <- coupled_init( nparticles,nobs,coupled_state=F)
  chain_state1 <- initial_conditions$chain_state1
  log_pdf_state1 <- initial_conditions$log_pdf_state1
  state_crn1 <- initial_conditions$state_crn1

  chain_state2 <- initial_conditions$chain_state2
  log_pdf_state2 <- initial_conditions$log_pdf_state2
  state_crn2 <- initial_conditions$state_crn2


  p <- length(chain_state1)
  samples1 <- matrix(nrow = K+preallocate+1, ncol = p)
  nrowsamples1 <- K+preallocate+1
  samples2 <- matrix(nrow = K+preallocate, ncol = p)
  log_pdf1 <- matrix(nrow = K+preallocate+1, ncol = p)
  log_pdf2 <- matrix(nrow = K+preallocate, ncol = p)

  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  log_pdf1[1,] <- log_pdf_state1
  log_pdf2[1,] <- log_pdf_state2


  current_nsamples1 <- 1
  iter <- 0

  # update one step
  for(i in 1:iters_per_cycle){
    cpm_step <- single_kernel(chain_state1,state_crn1,log_pdf_state1,i)
    chain_state1 <- cpm_step$chain_state
    log_pdf_state1 <- cpm_step$log_pdf_state
    state_crn1 <- cpm_step$state_crn
  }


  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1
  log_pdf1[current_nsamples1,] <- log_pdf_state1

  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){

      # update one step
      cpm_step <- single_kernel(chain_state1,state_crn1,log_pdf_state1,iter)

      chain_state1 <- cpm_step$chain_state
      log_pdf_state1 <- cpm_step$log_pdf_state
      state_crn1 <- cpm_step$state_crn

      chain_state2 <- cpm_step$chain_state
      log_pdf_state2 <- cpm_step$log_pdf_state
      state_crn2 <- cpm_step$state_crn
    } else {
      cpm_step <- coupled_kernel(chain_state1,chain_state2,
                                 state_crn1, state_crn2,
                                 log_pdf_state1,log_pdf_state2,
                                 iter)

      chain_state1 = cpm_step$chain_state1
      log_pdf_state1 = cpm_step$log_pdf_state1
      state_crn1 = cpm_step$state_crn1

      chain_state2 = cpm_step$chain_state2
      log_pdf_state2 = cpm_step$log_pdf_state2
      state_crn2 = cpm_step$state_crn2

      if (all(chain_state1 == chain_state2) && all(state_crn1==state_crn2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1 - 1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
      log_pdf1 <- rbind(log_pdf1, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf1)))
      log_pdf2 <- rbind(log_pdf2, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf2)))
    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    log_pdf1[current_nsamples1+1,] <- log_pdf_state1
    log_pdf2[current_nsamples1,] <- log_pdf_state2

    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(K, tau) steps
    if (iter >= max(meetingtime, K)){
      finished <- TRUE
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  log_pdf1 <- log_pdf1[1:current_nsamples1,,drop=F]
  log_pdf2 <- log_pdf2[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2, log_pdf1=log_pdf1, log_pdf2=log_pdf2,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}




# from coupled_chains -----------------------------------------------------
# Run coupled chains until max(tau, K) where tau is the meeting time and K specified by user
#'@rdname coupled_cpm_chains_Tblocking
#'@title Coupled correlated pseudo-marginal chains with blocking
#'@description Sample two CPM chains, each following \code{single_kernel} marginally,
#' and \code{coupled_kernel} jointly, until min(max(tau, K), max_iterations), where tau
#' is the first time at which the two chains meet (i.e. take the same value exactly).
#' Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}. The chains
#' are initialized from the distribution provided in \code{rinit}.
#'
#'
#'@param single_kernel function taking a state (in a vector) and an iteration, and returning
#' a list with a key named \code{chain_state} and containing the next state.
#'@param coupled_kernel function taking two states (in two vectors) and an iteration,
#'and returning a list with keys \code{chain_state1} and \code{chain_state2}.
#'@param coupled_init function taking no arguments are returning an initial state for a Markov chain.
#'@param K number of iterations desired (will be proportional to the computing cost if meeting occurs before \code{K},
#' default to 1).
#'@param max_iterations number of iterations at which the function stops if it is still running  (default to Inf).
#'@param preallocate  expected number of iterations, used to pre-allocate memory (default to 10).
#'@export
coupled_cpm_chains_Tblocking <- function(single_kernel, coupled_kernel, coupled_init,nparticles,nobs,iters_per_cycle=1,K = 1, max_iterations = Inf, preallocate = 10,coupled_state=F,verbose=F, finish_time=NULL){
  if((max_iterations%%iters_per_cycle)!=0){
    stop('Error max_iterations must be a multiple of iters_per_cycle')
  }

  early_stopping <- F

  initial_conditions <- coupled_init( nparticles,nobs,coupled_state=F)
  chain_state1 <- initial_conditions$chain_state1
  log_pdf_state1 <- initial_conditions$log_pdf_state1
  state_crn1 <- initial_conditions$state_crn1
  loglik_t1 <- initial_conditions$loglik_t1
  chain_state2 <- initial_conditions$chain_state2
  log_pdf_state2 <- initial_conditions$log_pdf_state2
  state_crn2 <- initial_conditions$state_crn2
  loglik_t2 <- initial_conditions$loglik_t2


  p <- length(chain_state1)
  samples1 <- matrix(nrow = K+preallocate+1, ncol = p)
  nrowsamples1 <- K+preallocate+1
  samples2 <- matrix(nrow = K+preallocate, ncol = p)
  log_pdf1 <- matrix(nrow = K+preallocate+1, ncol = p)
  log_pdf2 <- matrix(nrow = K+preallocate, ncol = p)
  prop_coupled <- matrix(nrow = K+preallocate,ncol=1)
  theta_coupled <- matrix(nrow = K+preallocate,ncol=1)


  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  log_pdf1[1,] <- log_pdf_state1
  log_pdf2[1,] <- log_pdf_state2


  current_nsamples1 <- 1
  iter <- 0

  # update one step
  for(i in 1:iters_per_cycle){
    cpm_step <- single_kernel(chain_state1,state_crn1,log_pdf_state1,loglik_t1,i)
    chain_state1 <- cpm_step$chain_state
    log_pdf_state1 <- cpm_step$log_pdf_state
    state_crn1 <- cpm_step$state_crn
    loglik_t1 <- cpm_step$loglik_t
  }


  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1
  log_pdf1[current_nsamples1,] <- log_pdf_state1

  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations && !early_stopping){
    iter <- iter + 1
    if(verbose){
      print(iter)
    }
    if (meet){

      # update one step
      cpm_step <- single_kernel(chain_state1,state_crn1,log_pdf_state1,loglik_t1,iter)

      chain_state1 <- cpm_step$chain_state
      log_pdf_state1 <- cpm_step$log_pdf_state
      state_crn1 <- cpm_step$state_crn
      loglik_t1 <- cpm_step$loglik_t

      chain_state2 <- cpm_step$chain_state
      log_pdf_state2 <- cpm_step$log_pdf_state
      state_crn2 <- cpm_step$state_crn
      loglik_t2 <- cpm_step$loglik_t

      theta_prop_coupled <- T
    } else {
      cpm_step <- coupled_kernel(chain_state1,chain_state2,
                                 state_crn1, state_crn2,
                                 log_pdf_state1,log_pdf_state2,
                                 loglik_t1,loglik_t2,
                                 iter)

      chain_state1 = cpm_step$chain_state1
      log_pdf_state1 = cpm_step$log_pdf_state1
      state_crn1 = cpm_step$state_crn1
      loglik_t1 <- cpm_step$loglik_t1
      theta_prop1 <- cpm_step$theta_prop1

      chain_state2 = cpm_step$chain_state2
      log_pdf_state2 = cpm_step$log_pdf_state2
      state_crn2 = cpm_step$state_crn2
      loglik_t2 <- cpm_step$loglik_t2
      theta_prop2 <- cpm_step$theta_prop2

      theta_prop_coupled <- all(theta_prop1==theta_prop2)


      if (all(chain_state1 == chain_state2) && all(state_crn1==state_crn2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1 - 1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
      log_pdf1 <- rbind(log_pdf1, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf1)))
      log_pdf2 <- rbind(log_pdf2, matrix(NA, nrow = new_rows, ncol = ncol(log_pdf2)))
      prop_coupled <- rbind(prop_coupled, matrix(NA, nrow = new_rows, ncol = ncol(prop_coupled)))
      theta_coupled <- rbind(theta_coupled, matrix(NA, nrow = new_rows, ncol = ncol(theta_coupled)))

    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    log_pdf1[current_nsamples1+1,] <- log_pdf_state1
    log_pdf2[current_nsamples1,] <- log_pdf_state2
    prop_coupled[current_nsamples1,] <- mean(state_crn1==state_crn2)
    theta_coupled[current_nsamples1,] <- theta_prop_coupled

    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(K, tau) steps
    if (iter >= max(meetingtime, K)){
      finished <- TRUE
    }

    if(!is.null(finish_time)){
      if(Sys.time()>finish_time){
        early_stopping <- T
      }
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  log_pdf1 <- log_pdf1[1:current_nsamples1,,drop=F]
  log_pdf2 <- log_pdf2[1:(current_nsamples1-1),,drop=F]
  prop_coupled <- prop_coupled[1:(current_nsamples1-1),,drop=F]
  theta_coupled <- theta_coupled[1:(current_nsamples1-1),,drop=F]

  return(list(samples1 = samples1, samples2 = samples2, log_pdf1=log_pdf1, log_pdf2=log_pdf2,
              meetingtime = meetingtime, iteration = iter, finished = finished,
              prop_coupled = prop_coupled,theta_coupled=theta_coupled))
}









# This file consists of the general-purpose functions coupled_chains, continue_coupled_chains,
# and H_bar, which implement our debiased mcmc algorithm for general kernels and functions h(.)

# from coupled_chains -----------------------------------------------------
# Run coupled chains until max(tau, K) where tau is the meeting time and K specified by user
#'@rdname coupled_mcmc_chains
#'@title Coupled MCMC chains
#'@description sample two MCMC chains, each following 'single_kernel' marginally,
#' and 'coupled_kernel' jointly, until min(max(tau, K), max_iterations), where tau
#' is the first time the two chains meet. Or more precisely, they meet with a delay of one, i.e. X_t = Y_{t-1}.
#'@export
coupled_mcmc_chains <- function(single_kernel, coupled_kernel, rinit, ..., K = 1, max_iterations = Inf, preallocate = 10){
  chain_state1 <- rinit()
  chain_state2 <- rinit()
  p <- length(chain_state1)
  samples1 <- matrix(nrow = K+preallocate+1, ncol = p)
  nrowsamples1 <- K+preallocate+1
  samples2 <- matrix(nrow = K+preallocate, ncol = p)
  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  current_nsamples1 <- 1
  chain_state1 <- single_kernel(chain_state1, ...)$chain_state
  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      chain_state1 <- single_kernel(chain_state1, ...)$chain_state
      chain_state2 <- chain_state1
    } else {
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, iter)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    if ((current_nsamples1+1) > nrowsamples1){
      # print('increase nrow')
      new_rows <- nrowsamples1 - 1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(K, tau) steps
    if (iter >= max(meetingtime, K)){
      finished <- TRUE
    }
  }
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}



## function to continue coupled chains until step K
## c_chain should be the output of coupled_chains
## and K should be more than c_chain$iteration, otherwise returns c_chain
#'@rdname continue_coupled_chains
#'@title Continue coupled MCMC chains up to K steps
#'@description ## function to continue coupled chains until step K
#' c_chain should be the output of coupled_chains
#' and K should be more than c_chain$iteration, otherwise returns c_chain
#'@export
continue_coupled_chains <- function(c_chain, single_kernel, K = 1, ...){
  if (K <= c_chain$iteration){
    ## nothing to do
    return(c_chain)
  } else {
    niterations <- K - c_chain$iteration
    chain_state1 <- c_chain$samples1[c_chain$iteration+1,]
    p <- length(chain_state1)
    samples1 <- matrix(nrow = niterations, ncol = p)
    samples2 <- matrix(nrow = niterations, ncol = p)
    for (iteration in 1:niterations){
      chain_state1 <- single_kernel(chain_state1, iteration)$chain_state
      samples1[iteration,] <- chain_state1
      samples2[iteration,] <- chain_state1
    }
    c_chain$samples1 <- rbind(c_chain$samples1, samples1)
    c_chain$samples2 <- rbind(c_chain$samples2, samples2)
    c_chain$iteration <- K
    return(c_chain)
  }
}



# from h_bar --------------------------------------------------------------

#'@rdname H_bar
#'@title Compute unbiased estimators from coupled chains
#'@description Compute the proposed unbiased estimators, for each of the element
#'in the list 'c_chains'. The integral of interest is that of the function h,
#'which can be multivariate. The estimator uses the variance reduction technique
#'whereby the estimator is the MCMC average between times k and K, with probability
#'going to one as k increases.
#'@export
H_bar <- function(c_chains, h = function(x) x, k = 0, K = 1){
  maxiter <- c_chains$iteration
  if (k > maxiter){
    print("error: k has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  if (K > maxiter){
    print("error: K has to be less than the horizon of the coupled chains")
    return(NULL)
  }
  # test the dimension of h(X)
  p <- length(h(c_chains$samples1[1,]))
  h_of_chain <- apply(X = c_chains$samples1[(k+1):(K+1),,drop=F], MARGIN = 1, FUN = h)
  if (is.null(dim(h_of_chain))){
    h_of_chain <- matrix(h_of_chain, ncol = 1)
  } else {
    h_of_chain <- t(h_of_chain)
  }
  H_bar <- apply(X = h_of_chain, MARGIN = 2, sum)
  if (c_chains$meetingtime <= k + 1){
    # nothing else to add
  } else {
    deltas <- matrix(0, nrow = maxiter - k + 1, ncol = p)
    deltas_term <- rep(0, p)
    for (t in k:min(maxiter-1, c_chains$meetingtime-1)){ # t is as in the report, where the chains start at t=0
      coefficient <- min(t - k + 1, K - k + 1)
      delta_tp1 <- h(c_chains$samples1[t + 1 + 1,]) - h(c_chains$samples2[t+1,]) # the +1's are because R starts indexing at 1
      deltas_term <- deltas_term + coefficient * delta_tp1
    }
    H_bar <- H_bar + deltas_term
  }
  return(H_bar / (K - k + 1))
}
