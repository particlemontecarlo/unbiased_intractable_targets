#'@rdname rnorm_max_coupling
#'@title Maximal coupling of two univariate Normal distributions
#'@description Sample from maximal coupling of two univariate Normal distributions,
#'specified through their means and standard deviations.
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param sigma1 standard deviation of first distribution
#'@param sigma2 standard deviation of second distribution
#'
#'@export
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  x <- rnorm(1, mu1, sigma1)
  if (dnorm(x, mu1, sigma1, log = TRUE) + log(runif(1)) < dnorm(x, mu2, sigma2, log = TRUE)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rnorm(1, mu2, sigma2)
      reject <- (dnorm(y, mu2, sigma2, log = TRUE) + log(runif(1)) < dnorm(y, mu1, sigma1, log = TRUE))
    }
    return(c(x,y))
  }
}



# from gaussian_max_coupling ----------------------------------------------

#'@rdname gaussian_max_coupling_cholesky_R
#'@title Maximal coupling of two multivariate Normal distributions
#'@description Sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means, the Cholesky factors of their covariance matrices,
#'and the Cholesky factors of the inverse covariance matrices (i.e. the precision matrices).
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param Cholesky1 Cholesky factor of variance of first distribution
#'@param Cholesky2 Cholesky factor of variance of second distribution
#'@param Cholesky_inverse1 Cholesky factor of precision of first distribution
#'@param Cholesky_inverse2 Cholesky factor of precision of second distribution
#'
#'@export
gaussian_max_coupling_cholesky_R <- function(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2){
  # we need Cholesky <- chol(Sigma), not the transpose
  valid_Cholesky1 <- all(Cholesky1[lower.tri(Cholesky1)]==0)
  valid_Cholesky2 <- all(Cholesky2[lower.tri(Cholesky2)]==0)
  stopifnot(valid_Cholesky1, valid_Cholesky2)

  return(gaussian_max_coupling_cholesky(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2))
}

#'@rdname gaussian_max_coupling
#'@title Maximal coupling of two multivariate Normal distributions
#'@description Sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means and covariance matrices.
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param Sigma1 Variance of first distribution
#'@param Sigma2 Variance of second distribution
#'@export
gaussian_max_coupling <- function(mu1, mu2, Sigma1, Sigma2){
  return(gaussian_max_couplingC(mu1, mu2, Sigma1, Sigma2))
}






#'@rdname vec_gaussian_max_coupling
#'@title Maximal coupling of two normals (vectorised code)
#'@description Sample maximal couplings for gaussians with different means but same std. dev.
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param s std dev of couplings
#'@export
vec_gauss_max_coupling <- function(mu1,mu2,s){

  stopifnot(length(mu1)==length(mu2))
  D <- length(mu1)

  output <- rep(NA,D)

  X <- mu1 + s*rnorm(D)
  pX <- dnorm(X,mean=mu1,sd=s)
  W <- runif(D)*pX

  qX <- dnorm(X,mean=mu2,sd=s)

  output_mask <- W<qX

  output[output_mask] <- X[output_mask]

  not_output <- !output_mask
  t <- 1
  while(any(not_output)){
    n_not_output <- sum(not_output)
    Y_star <- mu2[not_output] + s*rnorm(n_not_output)

    qY_star <- dnorm(Y_star,mean=mu2[not_output],sd=s)
    pY_star <- dnorm(Y_star,mean=mu1[not_output],sd=s)

    W_star <- runif(n_not_output)*qY_star

    new_output_mask <- W_star>pY_star

    newly_output_mask <- rep(F,D)
    newly_output_mask[not_output][new_output_mask] <- T
    output[newly_output_mask] <- Y_star[new_output_mask]

    output_mask <- output_mask | newly_output_mask
    not_output <- !output_mask
    t <- t+1
  }

  return(list(X1=X ,X2=output,t=t))
}






#'@rdname vec_gauss_max_coupling_mat
#'@title Maximal coupling of two normals (vectorised code) where each row is maximally coupled with a diagonal covariance
#'@description Sample maximal couplings for gaussians with different means but same std. dev.
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param s std dev of couplings
#'@export
vec_gauss_max_coupling_mat <- function (mu1, mu2, s)
{
  mat_dim1 <- dim(mu1)
  mat_dim2 <- dim(mu2)
  stopifnot(all(mat_dim1==mat_dim2)&is.matrix(mu1))

  n <- mat_dim1[1]
  D <- mat_dim1[2]
  output <- matrix(NA,n, D)
  X <- mu1 + s * matrix(rnorm(n*D),ncol=D)
  log_pX <- rowSums(dnorm(X, mean = mu1, sd = s,log=T))
  log_W <- log(runif(n)) +log_pX
  log_qX <- rowSums(dnorm(X, mean = mu2, sd = s,log=T))
  output_mask <- log_W < log_qX
  output[output_mask,] <- X[output_mask,]
  not_output <- !output_mask
  t <- 1
  while (any(not_output)) {
    n_not_output <- sum(not_output)
    Y_star <- mu2[not_output,] + s * matrix(rnorm(n_not_output*D),ncol=D)
    log_qY_star <- rowSums(dnorm(Y_star, mean = mu2[not_output,], sd = s,log=T))
    log_pY_star <- rowSums(dnorm(Y_star, mean = mu1[not_output,], sd = s,log=T))
    log_W_star <- log(runif(n_not_output)) + log_qY_star
    new_output_mask <- log_W_star > log_pY_star
    newly_output_mask <- rep(F, n)
    newly_output_mask[not_output][new_output_mask] <- T
    output[newly_output_mask,] <- Y_star[new_output_mask,]
    output_mask <- output_mask | newly_output_mask
    not_output <- !output_mask
    t <- t + 1
  }
  return(list(X1 = X, X2 = output, t = t))
}

