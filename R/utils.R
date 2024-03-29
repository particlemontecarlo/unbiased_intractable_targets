
# This file consists of short utility functions for working with the MVN distribution,
# propensity scores, as well as the ggplot2-related  setmytheme.

# MCMC helpers ------------------------------------------------------------
#'@rdname acc_rate
#'@title acc_rate
#'@description mcmc acceptance rate
#'@export
acc_rate <- function(chain){
  if(is.vector(chain)){
    return (sum(diff(chain)!=0)/length(chain))
  }else{
    return (sum(apply(diff(chain)!=0,1,any))/nrow(chain))
  }
}

# from util_dmvnorm -------------------------------------------------------
#'@rdname fast_dmvnorm
#'@title fast_dmvnorm
#'@description fast_dmvnorm
#'@export
fast_dmvnorm <- function(x, mean, covariance){
  return(dmvnorm(x, mean, covariance))
}

#'@rdname fast_dmvnorm_chol_inverse
#'@title fast_dmvnorm_chol_inverse
#'@description fast_dmvnorm_chol_inverse
#'@export
fast_dmvnorm_chol_inverse <- function(x, mean, chol_inverse){
  return(dmvnorm_cholesky_inverse(x, mean, chol_inverse))
}


# from util_rmvnorm -------------------------------------------------------
#'@rdname fast_rmvnorm
#'@title fast_rmvnorm
#'@description fast_rmvnorm
#'@export
fast_rmvnorm <- function(nparticles, mean, covariance){
  return(rmvnorm(nparticles, mean, covariance))
}

#'@rdname fast_rmvnorm_chol
#'@title fast_rmvnorm_chol
#'@description fast_rmvnorm_chol
#'@export
fast_rmvnorm_chol <- function(nparticles, mean, chol){
  return(rmvnorm_cholesky(nparticles, mean, chol))
}


# from util_expit ---------------------------------------------------------
#'@rdname expit
#'@title expit
#'@description expit function
#'@export
expit <- function(z) 1 / (1 + exp(-z))

# from util_setmytheme ----------------------------------------------------

#'@rdname setmytheme
#'@title Customize graphical settings
#'@description This function customizes the theme used by ggplot.
#'@export
setmytheme <- function(){
  theme_set(theme_bw())
  theme_update(axis.text.x = element_text(size = 20),
               axis.text.y = element_text(size = 20),
               axis.title.x = element_text(size = 25, margin=margin(20,0,0,0)),
               axis.title.y = element_text(size = 25, angle = 90, margin = margin(0,20,0,0)),
               legend.text = element_text(size = 20),
               legend.title = element_text(size = 20),
               title = element_text(size = 30),
               strip.text = element_text(size = 25),
               strip.background = element_rect(fill="white"),
               legend.position = "bottom")
}
