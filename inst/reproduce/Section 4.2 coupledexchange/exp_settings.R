
ising_model_no_external <- function(){

  M <- 80
  K <- 1100
  max_iterations <- 1.1*K
  mh_prop <- 1e-4*diag(1,1)
  nmcmc <- 200000
  nrep <- 1000
  nrep_initial <- 1000
  nrep_serial <- 10
  external_field <- F

  datafile <- 'inst/reproduce/Section 4.2 coupledexchange/data/ising_data_noext.RData'
  serial_datafile <- 'inst/reproduce/Section 4.2 coupledexchange/data/ising_serial_noext.RData'
  coupled_datafile <- 'inst/reproduce/Section 4.2 coupledexchange/data/ising_coupled_noext.RData'

  image_res_folder <- 'inst/reproduce/Section 4.2 coupledexchange/data/'

  beta_critical <- 0.5*log(1+sqrt(2))

  J_true <- beta_critical*(1/2)
  theta_true <- J_true

  rinit <- function(){
    return(beta_critical*runif(1))
  }

  J_lim <- c(0,beta_critical)

  return(list(
    M=M,
    K=K,
    max_iterations=max_iterations,
    mh_prop=mh_prop,
    nmcmc=nmcmc,
    nrep=nrep,
    nrep_initial=nrep_initial,
    nrep_serial=nrep_serial,
    external_field=external_field,
    J_true=J_true,
    J_lim=J_lim,
    datafile=datafile,
    theta_true=theta_true,
    serial_datafile=serial_datafile,
    coupled_datafile=coupled_datafile,
    image_res_folder=image_res_folder,
    rinit=rinit

  ))
}
