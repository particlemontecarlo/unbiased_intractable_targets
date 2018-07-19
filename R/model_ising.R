#'@rdname get_ising
#'@title Square Ising model
#'@description This function returns a list with objects as
#'* d_log_lik, r_lik, log_prior, three functions necessary to perform the exchange algorithm
#'@return A list
#'@export
get_ising <- function(M,J_lim,external_field=T,H_lim=NA,test_model=F){

  stopifnot(length(J_lim)==2)
  if(external_field&is.na(H_lim)){
    stop('H_lim required if specifying external field')
  }

  update <- function(state, image, innov, J, H) {
    # Ising (parallel) update, updating two interlaced codings simultaneously.
    M <- dim(image)[1] # Assumes square image!
    line <- (1:M)+1
    work <- cbind(0, rbind(0, state, 0), 0)

    threshold <- 1 / (1 + exp(- 2 * (J * (work[line,line+1] + work[line,line-1] + work[line+1,line] + work[line-1,line]) + H * image)))
    state[1:M, 1:M] <- -1
    state[innov < threshold] <- +1

    return (state)
  }

  generate_block <- function(block_length, M) {
    # Generate block of innovations for Ising (parallel) updates.
    ilist <- list()
    for (iter in 1:block_length) {
      ilist <- c(list(matrix(data=runif(M * M), nrow=M, ncol=M)), ilist)
    }
    return (ilist)
  }


  cycle <- function(innovations, image, J=0.6, H=1.5) {
    # Cycle for Ising (parallel) CFTP.
    M <- dim(image)[1] # Assumes square image!

    upper <- matrix(data=+1, ncol=M, nrow=M)
    lower <- matrix(data=-1, ncol=M, nrow=M)

    for (innov in innovations) {
      upper <- update(upper, image, innov, J=J, H=H)
      lower <- update(lower, image, innov, J=J, H=H)
    }

    if (sum(upper - lower) != 0)
      return (NULL)
    else
      return (upper)
  }


  ising_cftp <- function (initial_time_range,image, J=0.6, H=1.5) {
    # Ising (parallel) CFTP.
    innovations <- generate_block(initial_time_range, dim(image)[1])
    result <- NULL
    while (is.null(result)) {
      innovations <- c(generate_block(length(innovations), dim(image)[1]), innovations)
      result <- cycle(innovations, image, J=J, H=H)
    }
    return (result)
  }


  sample_ising <- function(image,J,H){
    result0 <- ising_cftp(1,image,J=J, H=H)
    innov <- matrix(data=runif(M * M), nrow=M, ncol=M)
    result1 <- update(result0, image, innov, J=J, H=H)
    chess1 <- outer(1:M, 1:M, function(x, y) ((x+y)%%2))
    chess0 <- 1 - chess1
    image0 <- result0 * chess0 + result1 * chess1
    return(image0)
  }


  ising_likelihood <- function(ising_sample,true_image,J,H){
    # check dimensions are ok
    stopifnot(all(dim(ising_sample)==dim(true_image)))
    # get dimension of image
    M <- dim(ising_sample)[1]

    # estimate external field contribution
    log_external_field <- H*sum(ising_sample*true_image)

    # estimate ising model likelihood
    side_counts <- ising_sample[1:M,1:(M-1)]*ising_sample[1:M,2:M]
    right_counts <- cbind(side_counts,rep(0,M))
    left_counts <- cbind(rep(0,M),side_counts)
    updown_counts <- ising_sample[1:(M-1),1:M]*ising_sample[2:M,1:M]
    up_counts <- rbind(updown_counts,rep(0,M))
    down_counts <- rbind(rep(0,M),updown_counts)

    log_internal_field <- J*(sum(right_counts)+sum(left_counts)+sum(up_counts)+sum(down_counts))

    return(sum(log_internal_field+log_external_field))
  }


  test_ising_likelihood <- function(test_size){
    J <- 0.1
    H <- 0.1

    # more complicated testcase
    M <- test_size
    rgen1 <- matrix(runif(M^2),ncol=M)
    rgen2 <- matrix(runif(M^2),ncol=M)

    t1 <- 2*((rgen1>0.5)-0.5)
    t2 <- 2*((rgen2>0.5)-0.5)

    neighbour_mat <- matrix(NA,ncol=M,nrow=M)
    for(i in 1:M){
      for(j in 1:M){
        n1 <- i+1
        n2 <- i-1
        n3 <- j+1
        n4 <- j-1

        k_sum <- 0
        if(n1<=M & n1>0){
          k_sum <- k_sum + t1[i,j]*t1[n1,j]
        }
        if(n2<=M & n2>0){
          k_sum <- k_sum + t1[i,j]*t1[n2,j]
        }
        if(n3<=M & n3>0){
          k_sum <- k_sum + t1[i,j]*t1[i,n3]
        }
        if(n4<=M & n4>0){
          k_sum <- k_sum + t1[i,j]*t1[i,n4]
        }
        neighbour_mat[i,j] <- k_sum
      }
    }


    internal_test <- J*sum(neighbour_mat) + H*sum(t1*t2)
    like_est <- ising_likelihood(t1,t2,J,H)

    stopifnot(internal_test==like_est)
    print('test passed')
  }

  if(test_model){
    test_ising_likelihood(100)
  }


  im1 <- matrix(1,M,M)


  # specify densities
  d_log_lik <- function(y,theta){
    if(external_field){
      ll <- ising_likelihood(y,im1,theta[1],theta[2])
    }else{
      ll <- ising_likelihood(y,im1,theta,0)
    }

    return(ll)
  }

  r_lik <- function(theta){
    if(external_field){
      sample <- sample_ising(im1,theta[1],theta[2])
    }else{
      sample <- sample_ising(im1,theta,0)
    }

    return(sample)
  }

  log_prior <- function(theta){


    if(external_field){

      if((theta[1]<J_lim[1]) |(theta[1]>J_lim[2]) | (theta[2]<H_lim[1]) |(theta[2]>H_lim[2]) ){
        return(-Inf)
      }else{
        return(0)
      }

    }else{
      if((theta<J_lim[1]) |(theta>J_lim[2])){
        return(-Inf)
      }else{
        return(0)
      }
    }

  }

  return(list(d_log_lik=d_log_lik,r_lik=r_lik,log_prior=log_prior))
}
