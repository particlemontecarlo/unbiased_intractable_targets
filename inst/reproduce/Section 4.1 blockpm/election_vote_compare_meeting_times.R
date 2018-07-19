if (!require(debiasedpmcmc)){
  library(devtools)
  devtools::document()
}
library(debiasedpmcmc)
library(coda)
library(ggplot2)
library(doMC)
library(MASS)
library(foreign)
library(TruncatedNormal)
rm(list = ls())
set.seed(1)
setmytheme()




cores_requested <- min(100,detectCores()-1)
registerDoMC(cores = cores_requested)

source("inst/reproduce/Section 4.1 blockpm/exp_settings.R")
settings <- election_model_settings()
# settings
# cov_file <- settings$cov_file
election_datafile <- settings$election_datafile
save_mt_fname <- settings$save_mt_fname
nobs_selection <- settings$nobs_selection

N_arr_block_coupled <- settings$N_arr_block_coupled
nrep_max <- settings$nrep_max
max_iterations_block <- settings$max_iterations_block
N_arr_pm_coupled <- settings$N_arr_pm_coupled
max_iterations_pm <- settings$max_iterations_pm

mean_est <- settings$mean_est
cov_est <- settings$cov_est

load(file=election_datafile)

X_save <- X_save[1:nobs_selection,1:3]
y_save <- y_save[1:nobs_selection,1:3]
X <- matrix(c(t(X_save)))

y <- c(t(y_save))
y <- 1*(y==1)


T_length <- 3
nobs <- length(y)/T_length

# parameter ests
sigma_vec <- c(0.5,0.5,0.5)


theta <- c(-0.01266113,  0,0.95,0.95,0.95)
D_theta <- length(theta)



y_obs <- matrix(NA,ncol=T_length, nrow=nobs)
for(i in 1:T_length){
  y_obs[,i] <- y[seq(i,T_length*nobs,T_length)]
}



# load(file=cov_file)
# cov_est <- cov_est
# mean_est <- mean_est
# cov_est <- diag(5e-5,D_theta)
# mean_est <- theta


init_cov <- 0.01^2*diag(1,D_theta)

re_model <- get_mv_probit_model(mean_est,init_cov)
block_pm_loglikelihood <- re_model$block_pm_loglikelihood
block_pm_logtarget <- re_model$block_pm_logtarget
state_crn_sample <- re_model$state_crn_sample
latent_state <- re_model$latent_state
logprior <- re_model$logprior
rinit <- re_model$rinit
block_coupled_init <- re_model$block_coupled_init
pm_logtarget <- re_model$pm_logtarget
pm_coupled_init <- re_model$pm_coupled_init



mh_prop <- (2.38^2/D_theta)*cov_est
D_theta <- length(theta)



# estimate variance of the log-likelihood estimate
R <- 10000

var_logz_arr_block <- rep(NA,length(N_arr_block_coupled))
var_logz_arr_like_t <- matrix(NA,length(N_arr_block_coupled),nobs)

for(ss in 1:length(N_arr_block_coupled)){
  N <- N_arr_block_coupled[ss]
  # N <- 100*3
  print(sprintf('getting var logz block pm, %i particles',N))

  logz_est <- foreach(pp = 1:R) %dopar% {

    state_crn <- state_crn_sample(N,nobs)
    latent_values <- latent_state(mean_est,state_crn)
    logz <- block_pm_logtarget(mean_est,latent_values)

  }

  logz_t <-sapply(logz_est,function(x) x$loglik_t)
  logz_var_t <- apply(logz_t,1,var)
  var_logz_arr_block[ss] <- var(unlist(logz_est))
  var_logz_arr_like_t[ss,] <- logz_var_t
}


var_logz_arr <- rep(NA,length(N_arr_pm_coupled))

for(ss in 1:length(N_arr_pm_coupled)){
  N <- N_arr_pm_coupled[ss]
  # N <- 100*3
  print(sprintf('getting var logz pm, %i particles',N))

  logz_est <- foreach(pp = 1:R) %dopar% {

    state_crn <- state_crn_sample(N,nobs)
    latent_values <- latent_state(mean_est,state_crn)
    logz <- pm_logtarget(mean_est,latent_values)

  }
  var_logz_arr[ss] <- var(unlist(logz_est))
}


print('variance ll:')
print(var_logz_arr)



# run block pm
t1 <- Sys.time()

coupled_state <- F

# compare the distribution of meeting times
couple_block_res <- list()
for(ss in 1:length(N_arr_block_coupled)){

  N <- N_arr_block_coupled[ss]
  print(sprintf('getting meeting times block pm, %i particles',N))

  # marginal kernel with joint updates only
  cpm_kernels <- get_cpm_T_blocking(logtarget = block_pm_logtarget, Sigma_proposal=mh_prop, dimension=D_theta,
                                    rho_component = 0,joint_update=F,rho = 0,
                                    single_kernel_only = F)
  single_cpm_kernel_component <- cpm_kernels$kernel
  coupled_cpm_kernel_component <- cpm_kernels$coupled_kernel


  block_start_time <- Sys.time()
  block_end_time <- block_start_time + 0.5*60*60

  couple_block <- foreach(pp = 1:nrep_max) %dopar% {
    ran_experiment <- F
    t1 <- Sys.time()
    if(t1<block_end_time){
      couple_res <- coupled_cpm_chains_Tblocking(single_cpm_kernel_component, coupled_cpm_kernel_component, block_coupled_init,N,nobs,iters_per_cycle=1,K = 1, max_iterations = max_iterations_block,verbose=F)
      ran_experiment <- T
    }else{
      couple_res <- NA
    }
    couple_res_time <- Sys.time() - t1

    return(list(couple_res=couple_res,couple_res_time=couple_res_time,ran_experiment=ran_experiment))
  }

  couple_block_res[[ss]] <- couple_block
  save.image(file=save_mt_fname)
}

# hist(mt_lst[[1]])

# now do the same again but for the pseudo marginal algorithm



# run coupled chain for pm
pmmh_kernels <- get_pm_kernel(logtarget = pm_logtarget,state_crn_sample, latent_state, Sigma_proposal = mh_prop, dimension = D_theta)
single_pmmh_kernel <- pmmh_kernels$kernel
coupled_pmmh_kernel <- pmmh_kernels$coupled_kernel

pm_couple_block_res <- list()
for(ss in 1:length(N_arr_pm_coupled)){
  N <- N_arr_pm_coupled[ss]
  print(sprintf('getting meeting times pm, %i particles',N))

  block_start_time <- Sys.time()
  block_end_time <- block_start_time + 0.5*60*60
  # finish_time <- block_start_time + (1+1/6)*60*60


  couple_block <- foreach(pp = 1:nrep_max) %dopar% {
    ran_experiment <- F
    t1 <- Sys.time()
    if(t1<block_end_time){

      couple_res <- tryCatch({
        res <- coupled_pm_chains(single_pmmh_kernel, coupled_pmmh_kernel, pm_coupled_init,N,nobs,K=1,max_iterations=max_iterations_pm,coupled_state = coupled_state)#,finish_time = finish_time)
      }, error = function(err) {
        return(err)
      })

      ran_experiment <- T
    }else{
      couple_res <- NA
    }
    couple_res_time <- Sys.time() - t1

    return(list(couple_res=couple_res,couple_res_time=couple_res_time,ran_experiment=ran_experiment))
  }


  pm_couple_block_res[[ss]] <- couple_block
  save.image(file=save_mt_fname)

}






save.image(file=save_mt_fname)


load(save_mt_fname)


# get important information

ntau_arr <- rep(NA,length(N_arr_block_coupled))
mt_lst <- list()
timing_per_iter_lst <- list()
timing_arr_lst <- list()
block_ar <- rep(NA,length(N_arr_block_coupled))

for(ss in 1:length(N_arr_block_coupled)){
  block_res <- couple_block_res[[ss]]
  block_res_chains <- lapply(block_res,function(x) x$couple_res)


  finished_mask <- sapply(block_res,function(x) x$ran_experiment)
  error_mask <- sapply(block_res,function(x) any(class(x$couple_res) == "error"))

  ntau_arr[[ss]] <- sum(finished_mask)



  mt_arr <- sapply(block_res_chains[finished_mask & !error_mask],function(x) x$meetingtime)
  timing_arr <- sapply(block_res[finished_mask & !error_mask],function(x) as.numeric(x$couple_res_time,units='secs'))

  mt_lst[[ss]] <- mt_arr
  timing_arr_lst[[ss]] <- timing_arr
  timing_per_iter_lst[[ss]] <- mean(timing_arr/mt_arr)
  block_ar[ss] <- mean(sapply(block_res_chains[finished_mask & !error_mask],function(x) acc_rate(x$samples1)))


}




pm_ntau_arr <- rep(NA,length(N_arr_pm_coupled))
pm_err_arr <- rep(NA,length(N_arr_pm_coupled))
pm_mt_lst <- list()
pm_timing_arr_lst <- list()
pm_timing_per_iter_lst <- list()
pm_ar <- rep(NA,length(N_arr_pm_coupled))
for(ss in 1:length(N_arr_pm_coupled)){
  block_res <- pm_couple_block_res[[ss]]
  block_res_chains <- lapply(block_res,function(x) x$couple_res)
  finished_mask <- sapply(block_res,function(x) x$ran_experiment)
  pm_ntau_arr[ss] <- sum(finished_mask)


  error_mask <- sapply(block_res,function(x) any(class(x$couple_res) == "error"))
  pm_err_arr[ss] <- sum(error_mask)


  mt_arr <- sapply(block_res_chains[finished_mask & !error_mask],function(x) x$meetingtime)
  timing_arr <- sapply(block_res[finished_mask & !error_mask],function(x) as.numeric(x$couple_res_time,units='secs'))

  pm_mt_lst[[ss]] <- mt_arr
  pm_timing_arr_lst[[ss]] <- timing_arr
  pm_timing_per_iter_lst[[ss]] <- mean(timing_arr/mt_arr)
  pm_ar[ss] <- mean(sapply(block_res_chains[finished_mask & !error_mask],function(x) acc_rate(x$samples1)))
}





# plot some interesting results
mt_cost_weighted <- list()
for(i in 1:length(N_arr_block_coupled)){
  mt_cost_weighted[[i]] <- mt_lst[[i]]*N_arr_block_coupled[i]
}
boxplot(mt_cost_weighted)


# plot some interesting results
pm_mt_cost_weighted <- list()
for(i in 1:length(N_arr_pm_coupled)){
  pm_mt_cost_weighted[[i]] <- pm_mt_lst[[i]]*N_arr_pm_coupled[i]
}
boxplot(pm_mt_cost_weighted)

both_cost_weighted <- c(mt_cost_weighted,pm_mt_cost_weighted)

# plot meeting times
boxplot(c(mt_lst,pm_mt_lst),main='meeting times both')


# plot cost weighted meeting times
boxplot(both_cost_weighted,log='y',main='cost weighted mt')


# plot acceptance rates for both guys
barplot(rbind(block_ar*2,pm_ar),beside=T,legend=c('block','pm'),main='acc rate')





# do ggplots
# cost weighted
mt_df_lst <- list()
for(i in 1:length(N_arr_block_coupled)){
  N_i <- N_arr_block_coupled[i]/T_length
  mt_df_i <- data.frame(mt_lst[[i]]/2,mt_lst[[i]]*N_i)
  mt_df_i$N <-N_i
  mt_df_lst[[i]] <- mt_df_i
}
mt_df <- do.call(rbind,mt_df_lst)
names(mt_df) <- c('mt','mtcw','N')
mt_df$N <- as.factor(mt_df$N)
mt_df$alg <- 'block'

pm_mt_df_lst <- list()
for(i in 1:length(N_arr_block_coupled)){
  N_i <- N_arr_pm_coupled[i]/T_length
  mt_df_i <- data.frame(pm_mt_lst[[i]],pm_mt_lst[[i]]*N_i)
  mt_df_i$N <- N_i
  pm_mt_df_lst[[i]] <- mt_df_i
}
pm_mt_df <- do.call(rbind,pm_mt_df_lst)
names(pm_mt_df) <- c('mt','mtcw','N')
pm_mt_df$N <- as.factor(pm_mt_df$N)
pm_mt_df$alg <- 'pm'


joint_df <- rbind(mt_df,pm_mt_df)



library(latex2exp)

g <- ggplot(data=joint_df) +
  geom_boxplot(aes(x=N,y=mt,fill=alg),notch=T)+
  ylab(TeX('meeting time$') ) +
  scale_fill_manual(values = c("grey", "white"))+
  theme(legend.title=element_blank())
g
ggsave('mt.pdf',width=7,height=5)


g <- ggplot(data=joint_df) +
  geom_boxplot(aes(x=N,y=mtcw,fill=alg),notch=T)+
  scale_y_continuous(trans='log10') +
  ylab(TeX('cost $\\times$ \n meeting time') ) +
  scale_fill_manual(values = c("grey", "white"))+
  theme(legend.title=element_blank())
g

ggsave('mt_cw.pdf',width=7,height=5)







