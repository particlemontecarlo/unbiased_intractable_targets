#'@rdname get_mv_probit_model
#'@title Multivariate probit model
#'@description This function returns a list with objects for random effects model
#'@return A list
#'@export
get_mv_probit_model <- function(init_mu,init_cov){

  init_mu <- init_mu
  init_cov <- init_cov

  block_pm_logtarget <- function(theta,latent_values){

    logp <- logprior(theta)
    if(logp!=-Inf){
      pm_res <- block_pm_loglikelihood(theta,latent_values)
      loglik <- pm_res$ll
      loglik_t <- pm_res$ll_t
    }else{
      loglik <- -Inf
      loglik_t <- -Inf
    }

    return(list(ltarget=loglik+logp,loglik=loglik,loglik_t = loglik_t,lprior=logp))
  }


  block_pm_loglikelihood <- function(theta,latent_values){

    D1 <- dim(latent_values)[1]
    D2 <- dim(latent_values)[2]
    stopifnot(D1==(nobs))
    stopifnot(D2==(N))

    beta <- theta[1]

    Sigma_X <- SigmaVec2Mat(theta[3:length(theta)])

    beta_X <- rowSums(t(t(X)*beta))

    mean_mat <- matrix(NA,ncol=T_length, nrow=nobs)
    for(i in 1:T_length){
      mean_mat[,i] <- beta_X[seq(i,T_length*nobs,T_length)] + theta[2]
    }

    stopifnot(all(dim(mean_mat)==dim(y_obs)))

    log_W_est <- ghk_obs_vectorise(N,y_obs,mean_mat,latent_values,Sigma_X)

    return(list(ll=sum(log_W_est),ll_t=log_W_est))
  }


  # Test single chain first
  rinit <- function(){
    theta_prop <- init_mu + 1*sqrt(diag(init_cov))*rnorm(length(init_mu))
#
#     theta_prop_u <- init_mu + 2*sqrt(diag(init_cov))
#     theta_prop_l <- init_mu - 2*sqrt(diag(init_cov))

    while(logprior(theta_prop)==-Inf){
      theta_prop <- init_mu + 1*sqrt(diag(init_cov))*rnorm(length(init_mu))
    }
    return(theta_prop)
  }

  block_coupled_init <- function( nparticles,nobs,coupled_state=TRUE){
    loglik_t1 <- -Inf
    loglik_t2 <- -Inf
    while(any(loglik_t1==-Inf)| any(loglik_t2==-Inf)){

      chain_state1 <- rinit()
      chain_state2 <- rinit()

      state_crn1 <- state_crn_sample( nparticles,nobs)
      if(coupled_state){
        state_crn2 <- state_crn1
      }else{
        state_crn2 <- state_crn_sample( nparticles,nobs)
      }

      latent_state_values1 <- latent_state(chain_state1,state_crn1)
      latent_state_values2 <- latent_state(chain_state2,state_crn2)

      target_res1 <- block_pm_logtarget(chain_state1,latent_state_values1)
      target_res2 <- block_pm_logtarget(chain_state2,latent_state_values2)

      log_pdf_state1 <- target_res1$ltarget
      log_pdf_state2 <- target_res2$ltarget

      loglik1 <- target_res1$loglik
      loglik2 <- target_res2$loglik

      loglik_t1 <- target_res1$loglik_t
      loglik_t2 <- target_res2$loglik_t

    }

    return(list(chain_state1=chain_state1,
                chain_state2=chain_state2,
                log_pdf_state1=log_pdf_state1,
                log_pdf_state2=log_pdf_state2,
                state_crn1=state_crn1,
                state_crn2=state_crn2,
                loglik1=loglik1,
                loglik2=loglik2,
                loglik_t1=loglik_t1,
                loglik_t2=loglik_t2))
  }



  state_crn_sample <- function( nparticles,nobs){
    return(matrix(rnorm(nobs*nparticles),nrow=nobs))
  }

  latent_state <- function(theta,state_crn){
    return(pnorm(state_crn))
  }

  logprior <- function(theta){
    Sigma_X <- SigmaVec2Mat(theta[3:length(theta)])

    if(any(eigen(Sigma_X)$values<0)){
      log_prior_res <- -Inf
    }else{
      log_prior_res <- -0.5*sum(theta[1:2]^2)/10^2

    }

    return(log_prior_res)

  }

  # functions for standard pseudo marginal chains
  pm_loglikelihood <- function(theta,latent_values){

    pm_res <- block_pm_loglikelihood(theta,latent_values)
    loglik <- pm_res$ll

    return(loglik)
  }

  # posterior density function up to normalising constant
  pm_logtarget <- function(theta,latent_values){
    logp <- logprior(theta)

    if(logp!=-Inf){
      pm_ll <- pm_loglikelihood(theta,latent_values)
      pm_logt_val <- pm_ll + logp
    }else{
      pm_logt_val <- -Inf
    }
    return(pm_logt_val)
  }


  pm_coupled_init <- function(N,nobs,coupled_state=TRUE){
    chain_state1 <- rinit()
    chain_state2 <- rinit()

    state_rn1 <- state_crn_sample(N,nobs)
    if(coupled_state){
      state_rn2 <-state_rn1
    }else{
      state_rn2 <- state_crn_sample(N,nobs)
    }

    latent_state_values1 <- latent_state(chain_state1,state_rn1)
    latent_state_values2 <- latent_state(chain_state2,state_rn2)

    log_pdf_state1 <- pm_logtarget(chain_state1,latent_state_values1)
    log_pdf_state2 <- pm_logtarget(chain_state2,latent_state_values2)

    return(list(chain_state1=chain_state1,chain_state2=chain_state2,log_pdf_state1=log_pdf_state1,log_pdf_state2=log_pdf_state2))
  }

  SigmaVec2Mat <- function(sigma_vec){
    m1 <- diag(1, T_length, T_length)
    m1[lower.tri(m1, diag=F)] <- sigma_vec
    m2 <- t(m1)
    m2[lower.tri(m2, diag=F)] <- sigma_vec
    return(m2)
  }
  SigmaMat2Vec <- function(sigma_mat){
    sigma_vec <- sigma_mat[lower.tri(sigma_mat,diag=F)]
    return(sigma_vec)
  }

  cpm_model <- list(block_pm_loglikelihood = block_pm_loglikelihood,
                    block_pm_logtarget = block_pm_logtarget,
                    state_crn_sample = state_crn_sample,
                    latent_state = latent_state,
                    logprior = logprior,
                    rinit = rinit,
                    block_coupled_init = block_coupled_init,
                    pm_logtarget=pm_logtarget,
                    pm_coupled_init=pm_coupled_init)
  return(cpm_model)
}
