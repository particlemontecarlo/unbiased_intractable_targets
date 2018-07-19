#'@rdname get_iid_cpm
#'@title IID gaussian LVM for CPM
#'@description This function returns a list with objects for cpm iid
#'@return A list
#'@export
get_iid_cpm <- function(dimension){
  # posterior density function up to normalising constant
  pm_logtarget <- function(theta,latent_values){
    pm_res <- pm_loglikelihood(latent_values)
    loglik <- pm_res$ll
    loglik_t <- pm_res$ll_t
    logp <- logprior(theta)
    return(list(ltarget=loglik+logp,loglik=loglik,loglik_t = loglik_t,lprior=logp))
  }


  pm_loglikelihood <- function(latent_values){
    like_ests = rep(0,length(y))
    delta = y-latent_values

    log_Gs = -0.5*(log(2*pi*sigma2)+delta^2/sigma2)
    max_logGs = apply(log_Gs,1,FUN=max)
    rescaled_Gs = log_Gs-max_logGs

    like_ests = log(rowMeans(exp(rescaled_Gs)))+max_logGs

    return(list(ll=sum(like_ests),ll_t=like_ests))
  }


  # Test single chain first
  rinit <- function() rnorm(1, 0, sd = 2)

  coupled_init <- function( nparticles,nobs,coupled_state=TRUE){
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

    target_res1 <- pm_logtarget(chain_state1,latent_state_values1)
    target_res2 <- pm_logtarget(chain_state2,latent_state_values2)

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



  # posterior density function up to normalising constant
  mh_logtarget <- function(theta) mh_loglikelihood(theta) + logprior(theta)


  mh_loglikelihood <- function(theta){
    return(sum(dnorm(y, theta, sqrt(1 + sigma2), log = TRUE)))
  }

  state_crn_sample <- function( nparticles,nobs){
    return(matrix(rnorm(nobs* nparticles),nrow=nobs,ncol= nparticles))
  }

  latent_state <- function(theta,state_crn){
    return(state_crn+theta)
  }


  # prior function
  logprior <- function(theta){
    return(dnorm(theta, mean = 0, sd = var_prior^0.5, log = TRUE))
  }



  cpm_model <- list(pm_loglikelihood =pm_loglikelihood,
                   pm_logtarget = pm_logtarget,
                   state_crn_sample = state_crn_sample,
                   latent_state = latent_state,
                   logprior = logprior,
                   rinit = rinit,
                   coupled_init = coupled_init,
                   mh_logtarget = mh_logtarget,
                   mh_loglikelihood = mh_loglikelihood)
  return(cpm_model)
}
