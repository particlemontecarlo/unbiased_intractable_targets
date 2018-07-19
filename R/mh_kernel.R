#'@rdname get_mh_kernel
#'@title Get random walk Metropolis-Hastings kernels
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled MH kernels, with Normal random walks.
#' These kernels can then be used in the function \code{\link{coupled_chains}}.
#'@param logtarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_mh_kernel <- function(logtarget, Sigma_proposal, dimension){
  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))
  zeromean <- rep(0, dimension)


  # single kernel
  kernel <- function(chain_state, iteration){
    proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
    proposal_pdf <- logtarget(proposal_value)
    current_pdf <- logtarget(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))

    if (accept){
      chain_state <- proposal_value
      log_target <- proposal_pdf
    } else {
      chain_state <- chain_state
      log_target <- current_pdf
    }

    return(list(chain_state=chain_state,log_target=log_target))
  }

  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    # distance_ <- mean((chain_state1 - chain_state2)^2)
    # proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
    proposal_pdf1 <- logtarget(proposal1)
    proposal_pdf2 <- logtarget(proposal2)
    current_pdf1 <- logtarget(chain_state1)
    current_pdf2 <- logtarget(chain_state2)
    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }
    if (accept1){
      chain_state1 <- proposal1
    }
    if (accept2){
      chain_state2 <- proposal2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}


#'@rdname get_pm_kernel
#'@title Get pseudo-marginal random walk Metropolis-Hastings kernels
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled pseudo marginal kernels, with Normal random walks.
#' These kernels can then be used in the function \code{\link{coupled_chains}}.
#'@param logtarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@param state_crn_sampe function to sample the random numbers used to estimate the target
#'@param latent_state function to transform the random numbers into samples of the latent state
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_pm_kernel <- function(logtarget,state_crn_sample, latent_state, Sigma_proposal, dimension){

  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))
  zeromean <- rep(0, dimension)


  # single kernel
  kernel <- function(chain_state, log_pdf_state,iteration,N,nobs){
    proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]

    state_crn <- state_crn_sample(N,nobs)
    latent_state_values <- latent_state(proposal_value,state_crn)

    proposal_pdf <- logtarget(proposal_value,latent_state_values)
    current_pdf <- log_pdf_state

    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value,log_pdf_state = proposal_pdf))
    } else {
      return(list(chain_state = chain_state,log_pdf_state = log_pdf_state))
    }
  }

  coupled_kernel <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2, iteration,N,nobs, coupled_state=TRUE){
    # distance_ <- mean((chain_state1 - chain_state2)^2)
    # proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]

    state_rn1 <- state_crn_sample(N,nobs)
    if (coupled_state || chain_state1==chain_state2 ){
      state_rn2 <- state_rn1
    } else {
      state_rn2 <- state_crn_sample(N,nobs)
    }

    latent_state_values1 <- latent_state(proposal1,state_rn1)
    latent_state_values2 <- latent_state(proposal2,state_rn2)

    proposal_pdf1 <- logtarget(proposal1,latent_state_values1)
    proposal_pdf2 <- logtarget(proposal2,latent_state_values2)

    current_pdf1 <- log_pdf_state1
    current_pdf2 <- log_pdf_state2

    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }

    if (accept1){
      chain_state1 <- proposal1
      log_pdf_state1 <- proposal_pdf1
    }
    if (accept2){
      chain_state2 <- proposal2
      log_pdf_state2 <- proposal_pdf2
    }

    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}


#'@rdname get_pmmh_kernel
#'@title Get particle marginal Metropolis-Hastings kernel
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a function kernel
#' \code{kernel}
#'
#'@param estimate_pftarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_pmmh_kernel <- function(estimate_pftarget, Sigma_proposal, dimension){
  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))
  zeromean <- rep(0, dimension)


  # single kernel
  kernel <- function(chain_state,log_pdf_state,path_state,iteration,nparticles){

    current_pdf <- log_pdf_state

    proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]

    smc_res <- estimate_pftarget(proposal_value,nparticles)
    proposal_pdf <- smc_res$log_target
    proposal_path <- smc_res$path

    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value,log_pdf_state = proposal_pdf,path_state=proposal_path,proposal_pdf=proposal_pdf))
    } else {
      return(list(chain_state = chain_state,log_pdf_state = log_pdf_state,path_state=path_state,proposal_pdf=proposal_pdf))
    }
  }

  coupled_kernel <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,path_state1,path_state2, iteration,nparticles){
    # distance_ <- mean((chain_state1 - chain_state2)^2)
    # proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]

    smc_res1 <- estimate_pftarget(proposal1,nparticles)
    proposal_pdf1 <- smc_res1$log_target
    proposal_path1 <- smc_res1$path

    if( any(proposal1!=proposal2)){
      smc_res2 <- estimate_pftarget(proposal2,nparticles)
      proposal_pdf2 <- smc_res2$log_target
      proposal_path2 <- smc_res2$path
      coupled_prop = FALSE
    }else{
      proposal_pdf2 <- proposal_pdf1
      proposal_path2 <- proposal_path1
      coupled_prop = TRUE
    }

    current_pdf1 <- log_pdf_state1
    current_pdf2 <- log_pdf_state2

    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }

    if (accept1){
      chain_state1 <- proposal1
      log_pdf_state1 <- proposal_pdf1
      path_state1 <- proposal_path1
    }
    if (accept2){
      chain_state2 <- proposal2
      log_pdf_state2 <- proposal_pdf2
      path_state1 <- proposal_path1
    }

    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2,
                path_state1=path_state1,path_state2=path_state2,
                coupled_prop=coupled_prop))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}






#'@rdname get_pmmh_kernel_component_update
#'@title Get particle marginal Metropolis-Hastings kernel
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a function kernel
#' \code{kernel}
#'
#'@param estimate_pftarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_pmmh_kernel_component_update <- function(estimate_pftarget, Sigma_proposal, dimension){

  component_vars <- diag(Sigma_proposal)
  component_std_dev <- component_vars^0.5

  # single kernel
  kernel <- function(chain_state,log_pdf_state,path_state,iteration,nparticles){

    current_pdf <- log_pdf_state

    component_selected <- (iteration-1)%%dimension + 1
    proposal_value <- chain_state
    proposal_value[component_selected] <- proposal_value[component_selected] + rnorm(1)*component_std_dev[component_selected]

    smc_res <- estimate_pftarget(proposal_value,nparticles)
    proposal_pdf <- smc_res$log_target
    proposal_path <- smc_res$path

    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value,log_pdf_state = proposal_pdf,path_state=proposal_path,proposal_pdf=proposal_pdf))
    } else {
      return(list(chain_state = chain_state,log_pdf_state = log_pdf_state,path_state=path_state,proposal_pdf=proposal_pdf))
    }
  }

  coupled_kernel <- function(chain_state1, chain_state2,log_pdf_state1,log_pdf_state2,path_state1,path_state2, iteration,nparticles){
    # distance_ <- mean((chain_state1 - chain_state2)^2)
    # proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)

    component_selected <- (iteration-1)%%dimension + 1
    component_selected_std_dev <- component_std_dev[component_selected]
    Sigma1_chol<-matrix(component_selected_std_dev)
    Sigma2_chol<-matrix(component_selected_std_dev)
    Sigma1_chol_inv<-matrix(1/component_selected_std_dev)
    Sigma2_chol_inv<-matrix(1/component_selected_std_dev)

    chain_state1_selected <- chain_state1[component_selected]
    chain_state2_selected <- chain_state2[component_selected]

    proposal_selected <- gaussian_max_coupling_cholesky_R(chain_state1_selected, chain_state2_selected,
                                                       Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
    proposal_selected1 <- proposal_selected[,1]
    proposal_selected2 <- proposal_selected[,2]

    proposal1 <- chain_state1
    proposal2 <- chain_state2

    proposal1[component_selected] <- proposal_selected1
    proposal2[component_selected] <- proposal_selected2

    smc_res1 <- estimate_pftarget(proposal1,nparticles)
    proposal_pdf1 <- smc_res1$log_target
    proposal_path1 <- smc_res1$path

    if( any(proposal1!=proposal2)){
      smc_res2 <- estimate_pftarget(proposal2,nparticles)
      proposal_pdf2 <- smc_res2$log_target
      proposal_path2 <- smc_res2$path
    }else{
      proposal_pdf2 <- proposal_pdf1
      proposal_path2 <- proposal_path1
    }

    current_pdf1 <- log_pdf_state1
    current_pdf2 <- log_pdf_state2

    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }

    if (accept1){
      chain_state1 <- proposal1
      log_pdf_state1 <- proposal_pdf1
      path_state1 <- proposal_path1
    }
    if (accept2){
      chain_state2 <- proposal2
      log_pdf_state2 <- proposal_pdf2
      path_state1 <- proposal_path1
    }

    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2,
                path_state1=path_state1,path_state2=path_state2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}






#'@rdname get_cpm_kernel
#'@title Get correlated pseudo-marginal kernels
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a function kernel
#' \code{kernel}
#'Optionally the kernels can perform blocked gibbs updates also
#'@param logtarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_cpm_kernel <- function(logtarget, Sigma_proposal, dimension,joint_update=T,rho=NULL,dim_U=NULL,component_updates=F,rho_component=NULL,block_indices_vec=NULL,single_kernel_only=F){

  # require either joint or component updates
  if(!joint_update & !component_updates){
    stop('require either joint or component updates')
  }

  # require rho if joint_update is true
  if(joint_update & (is.null(rho) | is.null(dim_U))){
    stop('joint update requires value of rho and dimension of auxiliary random variables')
  }

  # cannot have block_indices_vec and rho_component without the other
  if(component_updates & (is.null(rho_component)|is.null(block_indices_vec))){
    stop('component_updates requires either both rho_component and block_indices_vec')
  }


  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))

  # perform error checking for block_indices_vec, this is used for the coupling strategies
  if(component_updates){
    if(!(all(block_indices_vec==rep(1:max(block_indices_vec),sum(block_indices_vec==1))))){
      stop('Equal block size implemented only')
    }
    blocksize <- sum(block_indices_vec==1)
    Sigma_chol_crn <- chol(diag(1-rho_component^2,blocksize,blocksize))
    Sigma_chol_crn_inv <- chol(diag(1/(1-rho_component^2),blocksize,blocksize))
  }

  # set chol matrices for joint update
  if(joint_update & !single_kernel_only){
    Sigma_chol_joint_crn <- chol(diag(1-rho^2,dim_U,dim_U))
    Sigma_chol_joint_crn_inv <- chol(diag(1/(1-rho^2),dim_U,dim_U))
  }

  zeromean <- rep(0, dimension)

  # single kernel
  kernel <- function(chain_state, state_crn, log_pdf_state,iteration){
    # state_crn is an nobs X nparticles matrix

    # get dimensions
    nobs <- dim(state_crn)[1]
    nparticles <- dim(state_crn)[2]

    if(component_updates){
      # get the relevant block mask
      nblocks <- max(block_indices_vec) + 1
      block_indx <- ((iteration-1)%%nblocks)+1
      block_mask <- block_indices_vec==block_indx

      if(sum(block_mask)>0){
        # propose new crn with autocorrelation parameter
        crnvec <- crn_mat2vec(state_crn)
        state_crn_prop_vec <- crnvec
        state_crn_prop_vec[block_mask] <- rho_component*state_crn[block_mask] + sqrt(1-rho_component^2)*rnorm(sum(block_mask))
        state_crn_prop <- crn_vec2mat(state_crn_prop_vec,nparticles)
        proposal_value <- chain_state
      }else{
        # propose new theta
        proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
        if(!joint_update){
          state_crn_prop <- state_crn
        }else{
          state_crn_prop <- rho*state_crn + sqrt(1-rho^2)*matrix(rnorm(nobs*nparticles),ncol=nparticles)
        }
      }
    }else{
      # is the blocking is null, the joint update must be true
      proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
      state_crn_prop <- rho*state_crn + sqrt(1-rho^2)*matrix(rnorm(nobs*nparticles),ncol=nparticles)
    }

    # map crn to latent state
    latent_state_values <- latent_state(proposal_value,state_crn_prop)

    # estimate log target
    target_res <- logtarget(proposal_value,latent_state_values)
    proposal_pdf <- target_res$ltarget
    lest_prop <- target_res$loglik
    current_pdf <- log_pdf_state

    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value,log_pdf_state = proposal_pdf,state_crn=state_crn_prop,lest_prop=lest_prop,theta_prop=proposal_value))
    } else {
      return(list(chain_state = chain_state,log_pdf_state = log_pdf_state,state_crn=state_crn,lest_prop=lest_prop,theta_prop=proposal_value))
    }
  }


  coupled_kernel <- function(chain_state1, chain_state2,
                             state_crn1,state_crn2,
                             log_pdf_state1,log_pdf_state2,
                             iteration){

    if(!all(dim(state_crn1)==dim(state_crn2))){
      stop('Incompatible CRN dimensions for coupled kernel')
    }

    # get dimensions
    nobs <- dim(state_crn1)[1]
    nparticles <- dim(state_crn1)[2]

    if(component_updates){
      # get the relevant block mask
      nblocks <- max(block_indices_vec) + 1
      block_indx <- ((iteration-1)%%nblocks)+1
      block_mask <- block_indices_vec==block_indx

      if(sum(block_mask)>0){
        # propose new crn with autocorrelation parameter
        crnvec1 <- crn_mat2vec(state_crn1)
        crnvec2 <- crn_mat2vec(state_crn2)

        state_crn_prop_vec1 <- crnvec1
        state_crn_prop_vec2 <- crnvec2

        mu1 <- rho_component*crnvec1[block_mask]
        mu2 <- rho_component*crnvec2[block_mask]

        crn_proposal_value <- gaussian_max_coupling_cholesky_R(mu1, mu2,
                                                               Sigma_chol_crn, Sigma_chol_crn, Sigma_chol_crn_inv, Sigma_chol_crn_inv)


        state_crn_prop_vec1[block_mask] <- crn_proposal_value[,1]
        state_crn_prop_vec2[block_mask] <- crn_proposal_value[,2]

        state_crn_prop1 <- crn_vec2mat(state_crn_prop_vec1,nparticles)
        state_crn_prop2 <- crn_vec2mat(state_crn_prop_vec2,nparticles)

        proposal1 <- chain_state1
        proposal2 <- chain_state2
      }else{
        # propose new theta
        proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                           Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
        proposal1 <- proposal_value[,1]
        proposal2 <- proposal_value[,2]

        if(!joint_update){
          state_crn_prop1 <- state_crn1
          state_crn_prop2 <- state_crn2
        }else{
          crnvec1 <- crn_mat2vec(state_crn1)
          crnvec2 <- crn_mat2vec(state_crn2)

          mu1 <- rho*crnvec1
          mu2 <- rho*crnvec2

          crn_proposal_value <- gaussian_max_coupling_cholesky_R(mu1, mu2,
                                                                 Sigma_chol_joint_crn, Sigma_chol_joint_crn,
                                                                 Sigma_chol_joint_crn_inv, Sigma_chol_joint_crn_inv)

          state_crn_prop_vec1 <- crn_proposal_value[,1]
          state_crn_prop_vec2 <- crn_proposal_value[,2]

          state_crn_prop1 <- crn_vec2mat(state_crn_prop_vec1,nparticles)
          state_crn_prop2 <- crn_vec2mat(state_crn_prop_vec2,nparticles)
        }
      }

    }else{
      # propose new theta
      proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                         Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
      proposal1 <- proposal_value[,1]
      proposal2 <- proposal_value[,2]

      crnvec1 <- crn_mat2vec(state_crn1)
      crnvec2 <- crn_mat2vec(state_crn2)

      mu1 <- rho*crnvec1
      mu2 <- rho*crnvec2

      crn_proposal_value <- gaussian_max_coupling_cholesky_R(mu1, mu2,
                                                             Sigma_chol_joint_crn, Sigma_chol_joint_crn,
                                                             Sigma_chol_joint_crn_inv, Sigma_chol_joint_crn_inv)

      state_crn_prop_vec1 <- crn_proposal_value[,1]
      state_crn_prop_vec2 <- crn_proposal_value[,2]

      state_crn_prop1 <- crn_vec2mat(state_crn_prop_vec1,nparticles)
      state_crn_prop2 <- crn_vec2mat(state_crn_prop_vec2,nparticles)
    }

    latent_state_values1 <- latent_state(proposal1,state_crn_prop1)
    latent_state_values2 <- latent_state(proposal2,state_crn_prop2)


    target_res1 <- logtarget(proposal1,latent_state_values1)
    proposal_pdf1 <- target_res1$ltarget
    lest_prop1 <- target_res1$loglik

    target_res2 <- logtarget(proposal2,latent_state_values2)
    proposal_pdf2 <- target_res2$ltarget
    lest_prop2 <- target_res2$loglik

    current_pdf1 <- log_pdf_state1
    current_pdf2 <- log_pdf_state2

    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(proposal_pdf1)){
      accept1 <- (logu < (proposal_pdf1 - current_pdf1))
    }
    if (is.finite(proposal_pdf2)){
      accept2 <- (logu < (proposal_pdf2 - current_pdf2))
    }

    if (accept1){
      chain_state1 <- proposal1
      state_crn1 <- state_crn_prop1
      log_pdf_state1 <- proposal_pdf1
    }
    if (accept2){
      chain_state2 <- proposal2
      state_crn2 <- state_crn_prop2
      log_pdf_state2 <- proposal_pdf2
    }

    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,state_crn1=state_crn1,state_crn2=state_crn2,log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2))
  }

  if(!single_kernel_only){
    return(list(kernel = kernel,coupled_kernel=coupled_kernel))
  }else{
    return(list(kernel = kernel))
  }

}






#'@rdname get_cpm_T_blocking
#'@title Get correlated pseudo-marginal kernels with T-blocking
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a function kernel
#' \code{kernel}
#'Optionally the kernels can perform blocked gibbs updates also
#'@param logtarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_cpm_T_blocking <- function(logtarget, Sigma_proposal, dimension,rho_component,joint_update=T,rho=NULL,single_kernel_only=T){


  # require rho if joint_update is true
  if(joint_update & (is.null(rho))){
    stop('joint update requires value of rho and dimension of auxiliary random variables')
  }


  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))


  zeromean <- rep(0, dimension)

  # single kernel
  kernel <- function(chain_state, state_crn, log_pdf_state,loglik_t,iteration){
    # state_crn is an nobs X nparticles matrix

    # get dimensions
    nobs <- dim(state_crn)[1]
    nparticles <- dim(state_crn)[2]

    # get the relevant block mask
    update_Us_locally <- iteration%%2==1
    if(update_Us_locally){
      # propose new crn with autocorrelation parameter
      state_crn_prop <- rho_component*state_crn + sqrt(1-rho_component^2)*matrix(rnorm(nobs*nparticles),ncol=nparticles)
      proposal_value <- chain_state #latent_state_values <- latent_state(chain_state,state_crn_prop)
    }else{
      # propose new theta
      proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
      if(joint_update){
        state_crn_prop <- rho*state_crn + sqrt(1-rho^2)*matrix(rnorm(nobs*nparticles),ncol=nparticles)
      }else{
        state_crn_prop <- state_crn
      }
    }


    latent_state_values <- latent_state(proposal_value,state_crn_prop)
    # estimate log target
    prop_res <- logtarget(proposal_value,latent_state_values)
    proposal_pdf <- prop_res$ltarget
    proposal_loglik_t <- prop_res$loglik_t
    lest_prop <- prop_res$loglik

    current_pdf <- log_pdf_state
    current_loglik_t <- loglik_t

    if(update_Us_locally){
      accepts <-  (log(runif(nobs)) < (proposal_loglik_t - current_loglik_t))
      state_crn_updated <- state_crn
      state_crn_updated[accepts,] <- state_crn_prop[accepts,]
      current_pdf_t_updated <- loglik_t
      current_pdf_t_updated[accepts] <- proposal_loglik_t[accepts]

      # map crn to latent state
      latent_state_values_updated <- latent_state(chain_state,state_crn_updated)
      target_res_updated <- logtarget(proposal_value,latent_state_values_updated)

      # specify return values
      ret_chain_state <- chain_state
      ret_state_crn <- state_crn_updated
      ret_log_pdf_state <- target_res_updated$ltarget
      ret_loglik_t <- target_res_updated$loglik_t
      ret_lest_prop <- target_res_updated$loglik
      ret_theta_prop <- chain_state
      ret_block_ar <- mean(accepts)
    }else{
      accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
      if (accept){
        ret_chain_state <- proposal_value
        ret_log_pdf_state <- proposal_pdf
        ret_state_crn <- state_crn_prop
        ret_loglik_t <- proposal_loglik_t
      } else {
        ret_chain_state <- chain_state
        ret_log_pdf_state <- log_pdf_state
        ret_state_crn <- state_crn
        ret_loglik_t <- loglik_t
      }
      ret_lest_prop <- lest_prop
      ret_theta_prop <- proposal_value
      ret_block_ar <- NA
    }



    return(list(chain_state = ret_chain_state,
                log_pdf_state = ret_log_pdf_state,
                state_crn = ret_state_crn,
                loglik_t = ret_loglik_t,
                lest_prop = ret_lest_prop,
                theta_prop = ret_theta_prop,
                block_ar = ret_block_ar))
  }




  coupled_kernel <- function(chain_state1, chain_state2,
                             state_crn1,state_crn2,
                             log_pdf_state1,log_pdf_state2,
                             loglik_t1,loglik_t2,
                             iteration){
    if(!all(dim(state_crn1)==dim(state_crn2))){
      stop('Incompatible CRN dimensions for coupled kernel')
    }

    # get dimensions
    nobs <- dim(state_crn1)[1]
    nparticles <- dim(state_crn1)[2]

    coupled <- all(state_crn1==state_crn2) & (chain_state1==chain_state2)

    # get the relevant block mask
    update_Us_locally <- iteration%%2==1
    if(update_Us_locally){
      # propose new crn with autocorrelation parameter
      mu1 <- rho_component*state_crn1
      mu2 <- rho_component*state_crn2
      coupling_res <- vec_gauss_max_coupling_mat(mu1,mu2,sqrt(1-rho_component^2))
      state_crn_prop1 <- coupling_res$X1
      state_crn_prop2 <- coupling_res$X2
      proposal1 <- chain_state1
      proposal2 <- chain_state2
    }else{
      # propose new theta
      # propose new theta
      proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                         Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
      proposal1 <- proposal_value[,1]
      proposal2 <- proposal_value[,2]

      if(joint_update){
        state_crn_prop1 <- rho*state_crn1 + sqrt(1-rho^2)*matrix(rnorm(nobs*nparticles),ncol=nparticles)
        if(!coupled){
          state_crn_prop2 <- rho*state_crn2 + sqrt(1-rho^2)*matrix(rnorm(nobs*nparticles),ncol=nparticles)
        }else{
          state_crn_prop2 <- state_crn_prop1
        }
      }else{
        state_crn_prop1 <- state_crn1
        state_crn_prop2 <- state_crn2
      }
    }


    # convert proposed crn to latent states
    latent_state_values1 <- latent_state(proposal1,state_crn_prop1)
    latent_state_values2 <- latent_state(proposal2,state_crn_prop2)

    # estimate log target with proposed theta and states
    prop_res1 <- logtarget(proposal1,latent_state_values1)
    proposal_pdf1 <- prop_res1$ltarget
    proposal_loglik_t1 <- prop_res1$loglik_t
    lest_prop1 <- prop_res1$loglik

    prop_res2 <- logtarget(proposal2,latent_state_values2)
    proposal_pdf2 <- prop_res2$ltarget
    proposal_loglik_t2 <- prop_res2$loglik_t
    lest_prop2 <- prop_res2$loglik

    # assign current values
    current_pdf1 <- log_pdf_state1
    current_loglik_t1 <- loglik_t1
    current_pdf2 <- log_pdf_state2
    current_loglik_t2 <- loglik_t2

    if(update_Us_locally){
      log_u <- log(runif(nobs))
      accepts1 <-  (log_u < (proposal_loglik_t1 - current_loglik_t1))
      accepts2 <-  (log_u < (proposal_loglik_t2 - current_loglik_t2))

      state_crn_updated1 <- state_crn1
      state_crn_updated1[accepts1,] <- state_crn_prop1[accepts1,]
      state_crn_updated2 <- state_crn2
      state_crn_updated2[accepts2,] <- state_crn_prop2[accepts2,]

      current_pdf_t_updated1 <- loglik_t1
      current_pdf_t_updated1[accepts1] <- proposal_loglik_t1[accepts1]
      current_pdf_t_updated2 <- loglik_t2
      current_pdf_t_updated2[accepts2] <- proposal_loglik_t2[accepts2]

      # map crn to latent state
      latent_state_values_updated1 <- latent_state(chain_state1,state_crn_updated1)
      target_res_updated1 <- logtarget(proposal1,latent_state_values_updated1)
      latent_state_values_updated2 <- latent_state(chain_state2,state_crn_updated2)
      target_res_updated2 <- logtarget(proposal2,latent_state_values_updated2)

      # specify return values
      ret_chain_state1 <- chain_state1
      ret_state_crn1 <- state_crn_updated1
      ret_log_pdf_state1 <- target_res_updated1$ltarget
      ret_loglik_t1 <- target_res_updated1$loglik_t
      ret_lest_prop1 <- target_res_updated1$loglik
      ret_theta_prop1 <- chain_state1
      ret_block_ar1 <- mean(accepts1)

      ret_chain_state2 <- chain_state2
      ret_state_crn2 <- state_crn_updated2
      ret_log_pdf_state2 <- target_res_updated2$ltarget
      ret_loglik_t2 <- target_res_updated2$loglik_t
      ret_lest_prop2 <- target_res_updated2$loglik
      ret_theta_prop2 <- chain_state2
      ret_block_ar2 <- mean(accepts2)
    }else{
      log_u <-log(runif(1))
      accept1 <- (log_u < (proposal_pdf1 - current_pdf1))
      accept2 <- (log_u < (proposal_pdf2 - current_pdf2))

      # deal with accepting first guy
      if (accept1){
        ret_chain_state1 <- proposal1
        ret_log_pdf_state1 <- proposal_pdf1
        ret_state_crn1 <- state_crn_prop1
        ret_loglik_t1 <- proposal_loglik_t1
      } else {
        ret_chain_state1 <- chain_state1
        ret_log_pdf_state1 <- log_pdf_state1
        ret_state_crn1 <- state_crn1
        ret_loglik_t1 <- loglik_t1
      }
      ret_lest_prop1 <- lest_prop1
      ret_theta_prop1 <- proposal1
      ret_block_ar1 <- NA

      # deal with accepting second guy
      if (accept2){
        ret_chain_state2 <- proposal2
        ret_log_pdf_state2 <- proposal_pdf2
        ret_state_crn2 <- state_crn_prop2
        ret_loglik_t2 <- proposal_loglik_t2
      } else {
        ret_chain_state2 <- chain_state2
        ret_log_pdf_state2 <- log_pdf_state2
        ret_state_crn2 <- state_crn2
        ret_loglik_t2 <- loglik_t2
      }
      ret_lest_prop2 <- lest_prop2
      ret_theta_prop2 <- proposal2
      ret_block_ar2 <- NA
    }

    return(list(chain_state1 = ret_chain_state1,
                chain_state2 = ret_chain_state2,
                state_crn1=ret_state_crn1,
                state_crn2=ret_state_crn2,
                log_pdf_state1=ret_log_pdf_state1,
                log_pdf_state2=ret_log_pdf_state2,
                loglik_t1 = ret_loglik_t1,
                loglik_t2 = ret_loglik_t2,
                theta_prop1 = ret_theta_prop1,
                theta_prop2 = ret_theta_prop2,
                block_ar1 = ret_block_ar1,
                block_ar2 = ret_block_ar2))
  }

  if(!single_kernel_only){
    return(list(kernel = kernel,coupled_kernel=coupled_kernel))
  }else{
    return(list(kernel = kernel))
  }

}




#'@rdname get_exchange_kernel
#'@title Get exchange kernels
#'@description This function takes descriptions of the prior, likelihood densities and samplers
#' and tuning parameters, and returns a function kernel
#' \code{kernel}
#'@param log_prior function to compute target log-prior
#'@param d_log_lik function to compute target likelihood log-density
#'@param r_log_lik function to sample target likelihood
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_exchange_kernel <- function(log_prior,d_log_lik,r_log_lik, Sigma_proposal, dimension){
  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))
  zeromean <- rep(0, dimension)
  log_prior <- log_prior



  log_accept <- function(theta_prop,chain_state){
    log_a_prior <- log_prior(theta_prop) - log_prior(chain_state)
    if(log_a_prior!=-Inf){

      w <- r_lik(theta_prop)
      log_a_likes <- d_log_lik(y,theta_prop) - d_log_lik(y,chain_state)
      log_a_simulated_likes <- d_log_lik(w,chain_state) - d_log_lik(w,theta_prop)
      log_a <- log_a_prior+log_a_likes+log_a_simulated_likes

      return(log_a)
    }else{
      return(-Inf)
    }
  }



  # single kernel
  kernel <- function(chain_state, iteration){
    proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
    log_accept_probability <- log_accept(proposal_value,chain_state)
    accept <- (log(runif(1)) < (log_accept_probability))

    if (accept){
      chain_state <- proposal_value
    } else {
      chain_state <- chain_state
    }

    return(list(chain_state=chain_state))
  }

  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    # distance_ <- mean((chain_state1 - chain_state2)^2)
    # proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]

    log_accept_probability1 <- log_accept(proposal1,chain_state1)
    if(all(proposal1==proposal2)){
      log_accept_probability2<-log_accept_probability1
    }else{
      log_accept_probability2 <- log_accept(proposal2,chain_state2)
    }



    logu <- log(runif(1))
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(log_accept_probability1)){
      accept1 <- (logu < log_accept_probability1)
    }
    if (is.finite(log_accept_probability2)){
      accept2 <- (logu < log_accept_probability2)
    }
    if (accept1){
      chain_state1 <- proposal1
    }
    if (accept2){
      chain_state2 <- proposal2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}



