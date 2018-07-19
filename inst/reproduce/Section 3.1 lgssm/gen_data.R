# installs the rcpp - useful if R CMD INSTALL requires more priveleges
if (!require(debiasedpmcmc)){
  library(devtools)
  devtools::document()
}
library(debiasedpmcmc)
library(ggplot2)
rm(list = ls())
set.seed(17)
setmytheme()


source("inst/exp_settings.R")



##
## Try a possible coupling with correlated latent variables
##
## Settings
settings <- lgssm_model_2params_100obs()

nobservations<- settings$nobservations
dimension<- settings$dimension
mu_0<- settings$mu_0
Sigma_0<- settings$Sigma_0
theta<- settings$theta
data_file <- settings$data_file
D_theta <- settings$D_theta
sigma_y <- settings$sigma_y


nobs_gendata <- nobservations*5

#
x <- matrix(0, nrow = nobs_gendata+1, ncol = dimension)
y <- matrix(0, nrow = nobs_gendata, ncol = dimension)
x[1,] <- fast_rmvnorm(1, mu_0, Sigma_0)
for (t in 1:nobs_gendata){
  x[t+1,] <- x[t,,drop=F] %*% diag(theta[1], dimension, dimension)+ fast_rmvnorm_chol(1, rep(0, dimension), diag(theta[2], dimension, dimension))
  y[t,]  <- x[t+1,] + fast_rmvnorm_chol(1, rep(0, dimension), diag(sigma_y, dimension, dimension))
}

save(x,y, file = data_file)


