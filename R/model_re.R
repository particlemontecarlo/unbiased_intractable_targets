#'@rdname get_re_model
#'@title Random effects model of a stochastic frontier
#'@description This function returns a list with objects for random effects model
#'@return A list
#'@export
get_re_model <- function(init_mu,init_cov){

  init_mu <- init_mu
  init_cov <- init_cov
  beta_variables <- c('W1','W2','W3','W4','Q1','Q2','Q3','Q4','Q5','T')
  # posterior density function up to normalising constant


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
    nobs <- dim(latent_values)[1]
    nparticles <- dim(latent_values)[2]
    like_ests = rep(0,nobs)

    ll_i_mat <- matrix(NA,nobs,nparticles)

    alpha_param <- theta[1]
    beta_param <- theta[2:5]
    gamma_param <- theta[6:10]
    delta_param <- theta[11]
    lambda_param <- theta[12]
    sigma_param <- theta[13]
    sigma_w_param <- theta[14]

    beta_X_elements <- t(t(df[beta_variables])*c(beta_param,gamma_param,delta_param))
    beta_X_it <- rowSums(beta_X_elements)

    # crn_vals <- matrix(rnorm(nobs*N),ncols=N)

    for(i in 1:nparticles){

      w_i <- sigma_w_param*latent_values[,i][df$BANK]
      eps_it <- y - (alpha_param+w_i) - beta_X_it

      l1 <- log(2)
      l2 <- dnorm(eps_it,mean=0,sd=sigma_param,log=T)
      l3 <- pnorm(eps_it,mean=0,sd=sigma_param/lambda_param,log=T)

      ll <- l1 + l2  + l3
      ll_i_mat[,i] <- rowSums(matrix(ll,nrow=nobs,byrow=T))
    }
    ll_max <- apply(ll_i_mat,1,max)
    ll_normalised <- ll_i_mat - ll_max
    like_ests <- log(rowMeans(exp(ll_normalised))) + ll_max

    return(list(ll=sum(like_ests),ll_t=like_ests))
  }


  # Test single chain first
  rinit <- function(){
    theta_prop <- init_mu + 1*sqrt(diag(init_cov))*rnorm(length(theta_true))

    theta_prop_u <- init_mu + 2*sqrt(diag(init_cov))
    theta_prop_l <- init_mu - 2*sqrt(diag(init_cov))

    while(any(theta_prop[12:14]<0) | any(theta_prop>theta_prop_u) | any(theta_prop<theta_prop_l)){
      theta_prop <- init_mu + 1*sqrt(diag(init_cov))*rnorm(length(theta_true))
    }
    return(theta_prop)
  }

  block_coupled_init <- function( nparticles,nobs,coupled_state=TRUE){
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
    return(matrix(rnorm(nobs* nparticles),nrow=nobs,ncol= nparticles))
  }

  latent_state <- function(theta,state_crn){
    return(state_crn)
  }


  # prior function
  logprior <- function(theta){
    if(any(theta[12:14]<0)){
      log_prior_res <- -Inf
    }else{
      log_prior_res <- sum(dnorm(theta[1:11],mean=0,sd=10,log=T)) + sum(dexp(theta[12:13],rate=1,log=T))
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


