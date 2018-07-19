
election_model_settings <- function(){


  # cov_file <- 'inst/election_panel_data/cov_data.RData'
  election_datafile <- 'inst/reproduce/Section 4.1 blockpm/elect_data.RData'
  save_mt_fname <- 'inst/reproduce/Section 4.1 blockpm/data/bpm_pm_mt.RData'
  save_Hkm_fname <- 'inst/reproduce/Section 4.1 blockpm/data/bpm_pm_Hkm.RData'
  save_serial_fname <- 'inst/reproduce/Section 4.1 blockpm/data/bpm_serial.RData'
  save_serial_PMMH_fname <- 'inst/reproduce/Section 4.1 blockpm/data/pm_PMMH_serial.RData'

  nobs_selection <- 2000
  N_cov_est <- 10*3
  nmcmc_serial_cov_est <- 500000

  # meeting time settings
  N_arr_block_coupled <- c(5,10,20)*3
  nrep_max <- 10000
  max_iterations_block <- 50000

  N_arr_pm_coupled <- c(600,700,800)*3#c(550,700,850)*3
  max_iterations_pm <-50000

  N_block_Hkm <- 10*3
  N_pm_Hkm <- 700*3

  nrep_Hkm <- 2*cores_requested
  k <- 500
  K <- 5000
  K_large <- 20000

  nmcmc_serial <- 250000
  nmcmc_serial_PMMH <- 5000
  nrep_serial <- 20
  nrep_serial_PMMH <- 40

  mean_est <- c(0.01464665, -0.08091029,  0.99063885,  0.95055937,  0.96655334 )
  cov_est <- matrix(c(  6.964512e-05, -2.323131e-04,  1.292314e-06, -5.004271e-07, 3.558656e-07,
                       -2.323131e-04,  1.432540e-03, -5.199375e-06,  4.784542e-06, 1.553754e-06,
                        1.292314e-06, -5.199375e-06,  3.948509e-06,  4.837124e-06, 1.348633e-06,
                       -5.004271e-07,  4.784542e-06,  4.837124e-06,  4.056879e-05, 2.256164e-05,
                        3.558656e-07,  1.553754e-06,  1.348633e-06,  2.256164e-05, 2.352386e-05),nrow=length(mean_est))
  return(list(#cov_file=cov_file,
              election_datafile=election_datafile,
              save_serial_fname=save_serial_fname,
              save_Hkm_fname=save_Hkm_fname,
              nobs_selection=nobs_selection,
              N_cov_est=N_cov_est,
              nmcmc_serial_cov_est=nmcmc_serial_cov_est,
              save_mt_fname=save_mt_fname,
              N_arr_block_coupled=N_arr_block_coupled,
              nrep_max=nrep_max,
              max_iterations_block=max_iterations_block,
              N_arr_pm_coupled=N_arr_pm_coupled,
              max_iterations_pm=max_iterations_pm,
              nrep_Hkm=nrep_Hkm,
              nmcmc_serial_PMMH=nmcmc_serial_PMMH,
              k=k,
              K=K,
              nmcmc_serial=nmcmc_serial,
              nrep_serial=nrep_serial,
              nrep_serial_PMMH=nrep_serial_PMMH,
              mean_est=mean_est,
              cov_est=cov_est,
              N_block_Hkm=N_block_Hkm,
              N_pm_Hkm=N_pm_Hkm))


}
