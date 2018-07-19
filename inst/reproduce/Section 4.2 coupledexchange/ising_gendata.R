if (!require(debiasedpmcmc)){
  library(devtools)
  devtools::document()
}
library(debiasedpmcmc)
library(coda)
library(doMC)
library(MASS)
library(latex2exp)
library(ggplot2)
rm(list = ls())
set.seed(17)
setmytheme()

source("inst/reproduce/Section 4.2 coupledexchange/exp_settings.R")

##
##
## Settings
# settings <- lgssm_model()
#settings <- ising_model()
settings <- ising_model_no_external()

M <- settings$M
J_true <- settings$J_true
J_lim <- settings$J_lim
external_field <- settings$external_field
if(external_field){
  H_true <- settings$H_true
  ising_model <- get_ising(M,J_lim,external_field=external_field,H_lim = H_lim,test_model=T)

  data <- ising_model$r_lik(c(J_true,H_true))
}else{
  ising_model <- get_ising(M,J_lim,external_field=external_field,test_model=T)

  data <- ising_model$r_lik(J_true)
}
datafile <- settings$datafile



image(data)


save(data,file=datafile)
print(sprintf('written to %s ',datafile))





