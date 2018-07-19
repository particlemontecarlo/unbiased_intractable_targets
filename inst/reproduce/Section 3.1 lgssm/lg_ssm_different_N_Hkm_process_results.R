library(BH)
if (!require(debiasedpmcmc)){
  library(devtools)
  devtools::document()
}
library(debiasedpmcmc)
library(MASS)
library(latex2exp)
library(doMC)
library(coda)
library(ggplot2)
library(reshape2)
rm(list = ls())
set.seed(1)

# experiment settings
source("inst/reproduce/Section 3.1 lgssm/exp_settings.R")



settings <- lgssm_model_2params_100obs()



nobservations<- settings$nobservations
dimension<- settings$dimension
mu_0<- settings$mu_0
Sigma_0<- settings$Sigma_0
theta<- settings$theta
N_arr_Hkm<- settings$N_arr_Hkm
D_theta <- settings$D_theta
exp_name <- settings$exp_name
exp_name_serial <- settings$exp_name_serial
logprior <- settings$logprior
sigma_y <- settings$sigma_y
serial_data_file <- settings$serial_data_file
N_arr_serial <- settings$N_arr_serial
N_arr_mt<- settings$N_arr_mt
mt_data_file <- settings$mt_data_file
nrep_mt <- settings$nrep_mt
setmytheme()

cores_requested <- min(100,detectCores()-5)
registerDoMC(cores = cores_requested)

n_N_arr_Hkm <- length(N_arr_Hkm)

data_file <- settings$data_file
lowTcalibration_file <- settings$lowTcalibration_file
nrep <- 1000




load(data_file)
x <- as.matrix(x[1:(nobservations+1),],ncol=1)
y <- as.matrix(y[1:nobservations,],ncol=1)
print(sprintf('number of obs : %i',nobservations))


N_linecolors <- rainbow(n_N_arr_Hkm)



# coupling settings
n_rep_Hkm <- 10000
K <- 20000
max_iterations <- 1.1*K
m <- K

k_plot <- 200
K_plot <- 4000
n_bs <- 3



# data folders
datafolder <- ''          # folder containing the results of the unbiased estimators, having run lg_ssm_different_N_Hkm.R
mt_res_datafile <- ''     # folder containing the results of the meeting times, having run lg_ssm_different_N.R
serial_res_datafile <-''  # folder containing the results of the serial algorithm
fig_folder <- ''          # folder used to specify where the resulting figures are directed to



fig_folder <- 'inst/lg_ssm/2params/unif_prior/figs/'


# get meeting times
load(mt_res_datafile)

mh_mts <- sapply(mh_mt_batch,function(x) x$meetingtime)

n_mtN <- length(N_arr_mt)
pmmh_mts <- matrix(NA,n_mtN,nrep_mt)

for(i in 1:n_mtN){
  pmmh_mts[i,] <- sapply(pmmh_mt_batches[[i]],function(x) x$meetingtime)
}

tail_start <- 150
tail_mts <- list()
for( i in 1:n_mtN){
  mts_i <- pmmh_mts[i,]
  tail_mts[[i]] <- mts_i[mts_i>tail_start]
}
boxplot(tail_mts,log='y')
matplot(rep(k_plot,8),add=T,type='l',color=rgb(0,0,1))


pm_mt_quantiles1 <- apply(pmmh_mts,1,function(x) quantile(x,0.99))
pm_mt_quantiles2 <- apply(pmmh_mts,1,function(x) quantile(x,0.999))
pm_mt_quantiles3 <- apply(pmmh_mts,1,function(x) quantile(x,0.9999))
pm_mt_max <- apply(pmmh_mts,1,function(x) max(x))


# estimate survival probabilities
max_val <- 1000
t_vals <- 10^seq(2,3,length.out=20)
survival_probs <- matrix(NA,length(t_vals),n_mtN)

for(i in 1:length(t_vals)){
  survival_probs[i,] <- rowSums(pmmh_mts>t_vals[i])/nrep_mt
}

matplot(log(t_vals),log(survival_probs),type='l')
matplot(log(t_vals),log(1300*t_vals^-2),add=T,pch=1)
matplot(log(t_vals),log(1e10*t_vals^-6),add=T,pch=2)
grid()



name <- paste(fig_folder,'lgssm_mt.pdf',sep="")
ggplot_df = data.frame(N=N_arr_mt,q0=pm_mt_quantiles1,q1=pm_mt_quantiles2,q2=pm_mt_quantiles3)

ggplot_melt <- reshape2::melt((ggplot_df),id='N')
ggplot_melt$quantile <- as.factor(ggplot_melt$variable)
levels(ggplot_melt$quantile) <- c(0.99,0.999,0.9999)

breaks <- 10^(-10:10)
minor_breaks <- rep(1:9, 21)*(10^rep(-10:10, each=9))

g <- ggplot(data=ggplot_melt) +
  geom_point(aes(x=N,y=value,shape=quantile),cex=2) +
  geom_line(aes(x=N,y=value,linetype=quantile)) +
  scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks) +
  ylab('meeting time')
g
#  scale_y_log10()
ggsave(name,plot=g,width=7,height=5)
# get hkm estimators









avg_varests <- matrix(NA,length(N_arr_Hkm),1)
theta1_varests <- matrix(NA,length(N_arr_Hkm),1)
theta2_varests <- matrix(NA,length(N_arr_Hkm),1)
theta1_var_bs <- matrix(NA,length(N_arr_Hkm),n_bs)
theta2_var_bs <- matrix(NA,length(N_arr_Hkm),n_bs)
theta1_estimators <- matrix(NA,length(N_arr_Hkm),n_rep_Hkm)
theta2_estimators <- matrix(NA,length(N_arr_Hkm),n_rep_Hkm)

theta_h_estimators <- matrix(NA,length(N_arr_Hkm),n_rep_Hkm)
theta_h_varests <- matrix(NA,length(N_arr_Hkm),1)
h <- function(x){
  return(sum(x)+sum(x^2))
}
mt_mat <- matrix(NA,length(N_arr_Hkm),n_rep_Hkm)

for(n_i in 1:length(N_arr_Hkm)){
  nparticles <- N_arr_Hkm[n_i]
  print(sprintf('Analysing %i particles',nparticles))

  load(sprintf('%slgssm_Hkm_%i.RData',datafolder,n_i))
  # k_plot <- round(pm_mt_quantiles2[n_i])
  # K_plot <- k_plot*5


  print(sprintf('k = %i',k_plot))

  error_mask <- sapply(pmmh_batch,function(x)class(x)=='try-error')

  nrep_Hkm_i <- length(pmmh_batch)
  for(i in 1:nrep_Hkm_i){ if(error_mask[i]){print(sprintf('warning error occurred in iteration: %i',i))}  }

  meeting_times_initial <- sapply(pmmh_batch[!error_mask],function(x) x$meetingtime)
  mt_mat[n_i,] <- meeting_times_initial

  coupled_proposals <- sapply(pmmh_batch[!error_mask],function(x) x$coupled_prop)

  estimators_lst <- foreach(k_i = 1:nrep_Hkm_i) %dopar% {
    H_bar(pmmh_batch[[k_i]],k=k_plot,K=K_plot)
  }
  estimators <- matrix(unlist(estimators_lst),ncol=D_theta,byrow=T)
  var_estimators <- diag(cov(estimators))


  estimators_lst_h <- foreach(k_i = 1:nrep_Hkm_i) %dopar% {
    H_bar(pmmh_batch[[k_i]],h=h,k=k_plot,K=K_plot)
  }
  estimators_h <- unlist(estimators_lst_h)
  var_estimators_h <- var(estimators_h)


  theta1_varests[n_i,] <- var_estimators[1]
  theta2_varests[n_i,] <- var_estimators[2]
  theta_h_varests[n_i,] <- var_estimators_h
  theta1_estimators[n_i,] <- estimators[,1]
  theta2_estimators[n_i,] <- estimators[,2]
  theta_h_estimators[n_i,] <- estimators_h
}



# get the Hkm estimators for the mh algorithm

load(sprintf('%slgssm_Hkm_mh_batch.RData',datafolder))
# k_plot <- round(pm_mt_quantiles2[n_i])
# K_plot <- k_plot*5


print(sprintf('processing mh batch k = %i',k_plot))
error_mask <- sapply(mh_batch,function(x)class(x)=='try-error')

nrep_Hkm_i <- length(mh_batch)
for(i in 1:nrep_Hkm_i){ if(error_mask[i]){print(sprintf('warning error occurred in iteration: %i',i))}  }

meeting_times_initial <- sapply(mh_batch[!error_mask],function(x) x$meetingtime)
mt_mh <- meeting_times_initial

estimators_lst_mh <- foreach(i = 1:nrep_Hkm_i) %dopar% {
  H_bar(mh_batch[[i]],k=k_plot,K=K_plot)
}
estimators_mh <- matrix(unlist(estimators_lst_mh),ncol=D_theta,byrow=T)
var_estimators_mh <- diag(cov(estimators_mh))

estimators_lst_h_mh <- foreach(i = 1:nrep_Hkm_i) %dopar% {
  H_bar(mh_batch[[i]],h=h,k=k_plot,K=K_plot)
}
estimators_h_mh <- unlist(estimators_lst_h_mh)
var_estimators_h_mh <- var(estimators_h_mh)




# get the inefficiencies for the estimators for a few test functions
var_ests <- theta_h_varests
std_err <- var_ests*N_arr_Hkm*sqrt(2/(2000-1))
plot((var_ests*N_arr_Hkm),log='y')
matplot((var_ests*N_arr_Hkm)+2*std_err[1:8],type='l',add=T)
matplot(var_ests*N_arr_Hkm-2*std_err,type='l',add=T)


var_ests <- diag(cov(t(theta2_estimators)))
std_err <- var_ests*sqrt(2/(n_rep_Hkm-1))
v2_plt <- var_ests*N_arr_Hkm*K_plot
v2_u_plt <- (var_ests+2*std_err)*N_arr_Hkm*K_plot
v2_l_plt <- (var_ests-2*std_err)*N_arr_Hkm*K_plot
plot(v2_plt,log='y')
matplot(v2_u_plt,type='l',add=T)
matplot(v2_l_plt,type='l',add=T)



var_ests <- diag(cov(t(theta1_estimators)))
std_err <- var_ests*sqrt(2/(n_rep_Hkm-1))
v1_plt <- var_ests*N_arr_Hkm*K_plot
v1_u_plt <- (var_ests+2*std_err)*N_arr_Hkm*K_plot
v1_l_plt <- (var_ests-2*std_err)*N_arr_Hkm*K_plot
plot(v1_plt)
matplot(v1_u_plt,type='l',add=T)
matplot(v1_l_plt,type='l',add=T)


var_ests <- theta_h_varests
std_err <- var_ests*sqrt(2/(n_rep_Hkm-1))
vh_plt <- var_ests*N_arr_Hkm*K_plot
vh_u_plt <- (var_ests+2*std_err)*N_arr_Hkm*K_plot
vh_l_plt <- (var_ests-2*std_err)*N_arr_Hkm*K_plot
plot(v1_plt)
matplot(v1_u_plt,type='l',add=T)
matplot(v1_l_plt,type='l',add=T)




# plot certain figures
name1 <- 'var_est_combinedN.pdf'
name2 <- 'var_est_combinedN_costweighted.pdf'
name3 <- 'var_est_theta1.pdf'
name4 <- 'var_est_theta1_costweighted.pdf'
name_mt <- 'lgssm_mt.pdf'
name_serial <- 'ineff_serial.pdf'
name_bs <- 'ineff_bs.pdf'
name_bs_simple <- 'ineff_bs_simple.pdf'
name_ineff_compare <- 'ineff_compare.pdf'
name_ineff_raw <- 'ineff_raw.pdf'
name_ineff_N <- 'ineff_N.pdf'
name_ineff_N_zoomed <- 'ineff_N_zoom.pdf'


# plot of meeting times
tail_start <- 150
tail_mts_Hkm <- list()
for( i in 1:n_mtN){
  mts_i <- mt_mat[i,]
  tail_mts_Hkm[[i]] <- mts_i[mts_i>tail_start]
}
n_tails <- sapply(tail_mts_Hkm,length)
boxplot(tail_mts_Hkm,log='y')
matplot(rep(200,8),add=T,type='l')
matplot(pm_mt_quantiles3,add=T,type='l')



# ggplot for serial inefficiency
# load serial results and compare
load(serial_res_datafile)
n_serial <- length(N_arr_serial)

library(coda)
acfs <- matrix(NA,n_serial,D_theta)
acfs_h <- matrix(NA,n_serial,1)
for(i in 1:n_serial) {
  print(i)
  pmmh_run <- pmmh_runs[[i]]$pmcmc_chain
  nmcmc <- dim(pmmh_run)[1]
  acfs[i,] <- (spectrum0.ar(pmmh_run[(0.1*nmcmc):nmcmc,])$spec)
  acfs_h[i,] <- (spectrum0.ar(apply(pmmh_run,1,h)[(0.1*nmcmc):nmcmc])$spec)
}

acf1 <- acfs[,1]
acf2 <- acfs[,2]
acf_h <- acfs_h[,1]

N_arr_unique <- unique(N_arr_serial)
eff_serial_mat1 <- matrix(NA,length(N_arr_unique),n_serial/length(N_arr_unique))
eff_serial_mat2 <- matrix(NA,length(N_arr_unique),n_serial/length(N_arr_unique))
eff_serial_mat_h <- matrix(NA,length(N_arr_unique),n_serial/length(N_arr_unique))
for(i in 1:length(N_arr_unique)){
  eff_serial_mat1[i,] <- acf1[N_arr_serial==N_arr_unique[i]]*N_arr_unique[i]
  eff_serial_mat2[i,] <- acf2[N_arr_serial==N_arr_unique[i]]*N_arr_unique[i]
  eff_serial_mat_h[i,] <- acf_h[N_arr_serial==N_arr_unique[i]]*N_arr_unique[i]
}
eff_unique1 <- rowMeans(eff_serial_mat1)
eff_unique2 <- rowMeans(eff_serial_mat2)
eff_unique_h <- rowMeans(eff_serial_mat_h)

plot(N_arr_unique,eff_unique_h)

v2_ser <- eff_unique2
v2_u_ser <- eff_unique2+2*sqrt(diag(var(t(eff_serial_mat2))/10))
v2_l_ser <- eff_unique2-2*sqrt(diag(var(t(eff_serial_mat2))/10))

v1_ser <- eff_unique1
v1_u_ser <- eff_unique1+2*sqrt(diag(var(t(eff_serial_mat1))/10))
v1_l_ser <- eff_unique1-2*sqrt(diag(var(t(eff_serial_mat1))/10))

v_h_ser <- eff_unique_h
v_h_u_ser <-eff_unique_h+2*sqrt(diag(var(t(eff_serial_mat_h))/(n_serial/length(N_arr_unique) )))
v_h_l_ser <- eff_unique_h-2*sqrt(diag(var(t(eff_serial_mat_h))/(n_serial/length(N_arr_unique) )))


plot(v_h_ser)
matplot(v_h_u_ser ,pch=2,add=T)
matplot(v_h_l_ser ,pch=2,add=T)



ggplot_df <- data.frame(t(eff_serial_mat_h))
names(ggplot_df) <- N_arr_unique
ggplot_df_melt <- melt(ggplot_df)
ggplot_df_melt$N <- ggplot_df_melt$variable
ggplot_df_melt$N <- as.factor(ggplot_df_melt$N)

name <- sprintf('%s%s',fig_folder,name_serial)
g <- ggplot(data=ggplot_df_melt)+
  geom_boxplot(aes(x=N,y=value)) +
  xlab('N') +
  ylab('N x inefficiency')
g
ggsave(g,file=name,width=7,height=5)







# produce inefficiency plots for both the serial and parallel algorithm

# begin with the inefficiency of the parallel algorithm
name <- sprintf('%s%s',fig_folder,name_ineff_raw)

ggplot_df1 <- data.frame(N=N_arr_unique,ineff=v_h_ser/N_arr_unique,ineff_u=v_h_u_ser/N_arr_unique,ineff_l=v_h_l_ser /N_arr_unique)
ggplot_df1$type='serial'
ggplot_df2 <- data.frame(ineff=vh_plt/N_arr_Hkm,ineff_u=vh_u_plt/N_arr_Hkm,ineff_l=vh_l_plt/N_arr_Hkm,N=N_arr_Hkm)
ggplot_df2$type='coupled PMMH'

ggplot_df_coupled_all <- rbind(ggplot_df2)

g<- ggplot(ggplot_df_coupled_all, aes(x=N, y=ineff,pch=type)) +
  geom_errorbar(aes(ymin=ineff_l, ymax=ineff_u,linetype=type), width=5)+
  geom_line(aes(linetype=type))+
  geom_point(cex=2) +
  ylab('inefficiency')+
  geom_hline(aes(yintercept=var_estimators_h_mh*K_plot,linetype = "coupled MH") )+
  scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks)+
  guides(pch = "none")+
  theme(legend.title=element_blank())
g

ggsave(plot=g,file=name,width=7,height=5)


# plot the inefficiency compared to the serial algorithm
name <- sprintf('%s%s',fig_folder,name_ineff_N)

ggplot_df1 <- data.frame(N=N_arr_unique,ineff=v_h_ser,ineff_u=v_h_u_ser,ineff_l=v_h_l_ser )
ggplot_df1$type='serial'
ggplot_df2 <- data.frame(ineff=vh_plt,ineff_u=vh_u_plt,ineff_l=vh_l_plt,N=N_arr_Hkm)
ggplot_df2$type='coupled PMMH'
ggplot_df_both_all_ineff <- rbind(ggplot_df1,ggplot_df2)


g<- ggplot(ggplot_df_both_all_ineff, aes(x=N, y=ineff,pch=type)) +
  geom_errorbar(aes(ymin=ineff_l, ymax=ineff_u,linetype=type), width=5)+
  geom_line(aes(linetype=type))+
  geom_point(cex=2) +
  ylab('N x inefficiency')+
  scale_y_log10(breaks =breaks ,minor_breaks=minor_breaks)+
  theme(legend.title=element_blank())
g
ggsave(plot=g,file=name,width=7,height=5)




# plot the inefficiency compared to the serial algorithm though zoomed in
name <- sprintf('%s%s',fig_folder,name_ineff_N_zoomed)

ggplot_df <- rbind(ggplot_df1,ggplot_df2[ggplot_df2[,1]<5000,])

g<- ggplot(ggplot_df, aes(x=N, y=ineff,pch=type)) +
  geom_errorbar(aes(ymin=ineff_l, ymax=ineff_u,linetype=type), width=5)+
  geom_line(aes(linetype=type))+
  geom_point(cex=2) +
  coord_cartesian(ylim=c(550,3000))+
  ylab('N x inefficiency')+
  theme(legend.title=element_blank())
g

ggsave(plot=g,file=name,width=7,height=5)



