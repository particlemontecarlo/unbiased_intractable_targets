# Vectorised implementation of GHK algorithm
#'@rdname ghk_obs_vectorise
#'@title GHK algorithm (vectorised)
#'@description runs the GHK algorithm to obtain importance sampling estimates for the observations
#'
#'@param nsamp number of particles x dimension of each observation
#'@param y_obs matrix of observations, nobs x dimension of each observation
#'@param mean_mat matrix of means for each observation in y_obs
#'@param uniform_rvs nobs x nsamp uniform random numbers
#'@param Sigma covariance matrix used in logistic regression
#'@export
# code to perform vectorised GHK estimator
ghk_obs_vectorise <- function(nsamp,y_obs,mean_mat,uniform_rvs,Sigma,K_units=NULL){


  # takes a matrix of nobs x D_y observations
  # a matrix of means for each observation and a covariance matrix across D_y
  # the idea is to vectorise samples of GHK over the observations
  nobs <- dim(y_obs)[1]
  D_y <- dim(y_obs)[2]
  C_ghk<- t(chol(Sigma))
  stopifnot(dim(C_ghk)[1]==D_y)
  stopifnot(all(dim(uniform_rvs)==c(nobs,nsamp)))
  stopifnot((nsamp%%D_y)==0)
  nsamp_per_N <- nsamp/D_y

  large_matrix <- (nsamp*nobs)>=1e6
  if(is.null(K_units)){
    K_units <- 10
  }

  # infer likelihoods for first observation
  lower_truncation <- rep(-Inf,nobs)
  lower_truncation[(y_obs[,1]==1)] <- 0
  upper_truncation <- rep(Inf,nobs)
  upper_truncation[!(y_obs[,1]==1)] <- 0



  # want to loop over D_y
  uniform_rvs_i <- uniform_rvs[,1:nsamp_per_N]

  l_eta <- replicate(nsamp_per_N,(lower_truncation - mean_mat[,1])/C_ghk[1,1])
  u_eta <- replicate(nsamp_per_N,(upper_truncation - mean_mat[,1])/C_ghk[1,1])

  eta_arr <- array(NA,dim=c(nobs,nsamp_per_N,D_y))
  ll_res  <- array(NA,dim=c(nobs,nsamp_per_N,D_y))
  eta_arr[,,1] <- matrix(norminvp(c(uniform_rvs_i),c(l_eta),c(u_eta)),ncol=nsamp_per_N)


  ll_res[,,1] <- log(pnorm(u_eta)-pnorm(l_eta))

  for(d in 2:D_y){
    uniform_rvs_i <- uniform_rvs[,(1+nsamp_per_N*(d-1)):(d*nsamp_per_N)]


    eta_arr_tmp <- eta_arr[,,1]*C_ghk[d,1]
    if(d>2){
      for(d_i in 2:(d-1)){
        eta_arr_tmp <- eta_arr_tmp + eta_arr[,,d_i]*C_ghk[d,d_i]
      }
    }
    mean_d <- eta_arr_tmp + mean_mat[,d]

    lower_truncation <- rep(-Inf,nobs)
    lower_truncation[(y_obs[,d]==1)] <- 0
    upper_truncation <- rep(Inf,nobs)
    upper_truncation[!(y_obs[,d]==1)] <- 0


    l_eta <- (lower_truncation - mean_d)/C_ghk[d,d]
    u_eta <- (upper_truncation - mean_d)/C_ghk[d,d]

    uniform_rvs_vec <-c(uniform_rvs_i)
    l_eta_vec <- c(l_eta)
    u_eta_vec <- c(u_eta)
    if(!large_matrix){
      inverse_res_vec <- norminvp(uniform_rvs_vec,l_eta_vec,u_eta_vec)
    }else{
      vec_length <- length(uniform_rvs_vec)
      stopifnot((vec_length%%K_units)==0)

      unit_length <- vec_length/K_units
      inverse_res_vec <- rep(NA,vec_length)
      for(kk in 1:K_units){
        unit_mask <- (1+(kk-1)*unit_length):(kk*unit_length)
        inverse_res_vec[unit_mask] <- norminvp(uniform_rvs_vec[unit_mask],l_eta_vec[unit_mask],u_eta_vec[unit_mask])
      }
    }

    eta_arr_mat <- matrix(inverse_res_vec,ncol=nsamp_per_N)

    eta_arr[,,d] <- eta_arr_mat

    ll_res[,,d] <- log(pnorm(u_eta)-pnorm(l_eta))
  }

  ll_est_ti <- rowSums(ll_res,dim=2)
  ll_rownormalisers <- apply(ll_est_ti,1,max)

  ll_est <- log(rowMeans(exp(ll_est_ti-ll_rownormalisers))) + ll_rownormalisers
  ll_est[ll_rownormalisers==-Inf] <- -Inf

  return(ll_est)
}








