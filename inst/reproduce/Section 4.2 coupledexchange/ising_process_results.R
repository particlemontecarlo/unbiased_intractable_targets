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


## Settings
source("inst/reproduce/Section 4.2 coupledexchange/exp_settings.R")

settings <- ising_model_no_external()

M <- settings$M
K <- settings$K
max_iterations <- settings$max_iterations
mh_prop <- settings$mh_prop
nmcmc <- settings$nmcmc
nrep <- settings$nrep
nrep_initial <- settings$nrep_initial
J_true <- settings$J_true
J_lim <- settings$J_lim
external_field <- settings$external_field
datafile <- settings$datafile
coupled_datafile <- settings$coupled_datafile
serial_datafile <- settings$serial_datafile
theta_true <- settings$theta_true
D_theta <- length(theta_true)
rinit <- settings$rinit
image_res_folder <-settings$image_res_folder


ising_model <- get_ising(M,J_lim,external_field=external_field,test_model=T)
log_prior <- ising_model$log_prior
d_log_lik <- ising_model$d_log_lik
r_lik <- ising_model$r_lik

k <- 100

run_initial_kernels <- F

# load data
load(datafile)
y <- data

# parallel settings
cores_requested <- min(100,detectCores()-1)
registerDoMC(cores = cores_requested)


Sigma_proposal <- mh_prop
kernels <- get_exchange_kernel(log_prior,d_log_lik,r_log_lik, Sigma_proposal, dimension=D_theta)
kernel <- kernels$kernel
coupled_kernel <- kernels$coupled_kernel

# load results of coupling experiment
load(coupled_datafile)
nrep <- length(mh_batch)

print(sprintf('iterations %i',nrep))

mt <- sapply(mh_batch_initial,function(x) x$meetingtime)
mt_initial <- sapply(mh_batch,function(x) x$meetingtime)
batch_times_secs <- sapply(batch_times,function(x) as.numeric(x,units='secs'))

k <- 100
K <- 1000
# run coupling algorithm to look at meeting times
ests <- foreach(i = 1:nrep) %dopar% {
  return(H_bar(mh_batch[[i]],k=0.1*K,K=K))
}

# run coupling algorithm to look at meeting times
ests_x2 <- foreach(i = 1:nrep) %dopar% {
  return(H_bar(mh_batch[[i]],function(x) x^2,k=k,K=K))
}


ests_mat <- do.call(rbind,ests)
ests_x2_mat <- do.call(rbind,ests_x2)



name <-  'ising_J.pdf'
pdf(name,width=10,height=10)
hist(ests_mat[,1])
dev.off()
print(sprintf('written to %s',name))



histogram_c_chains <- function(c_chains, component, k, K, breaks = NULL, nclass = 30){
  nsamples <- length(c_chains)
  if (is.null(breaks)){
    breaks <- find_breaks(c_chains, component, nclass, k)
  }
  mids <- create_mids(breaks)
  width <- diff(breaks)[1]
  # ## compute histogram
  res__ <- foreach (ibreak = 2:length(breaks), .combine = rbind) %dopar% {
    estimators <- rep(0, nsamples)
    for (irep in 1:nsamples){
      estimators[irep] <- estimator_bin_R(c_chains[[irep]], component, breaks[ibreak-1], breaks[ibreak], k, K)
    }
    prop <- mean(estimators)
    sd_prop <- sd(estimators) / sqrt(nsamples)
    c(prop, sd_prop)
  }
  prop <- res__[,1]
  sd_prop <- res__[,2]
  return(list(mids = mids, breaks = breaks, proportions = prop, sd = sd_prop, width = width))
}


name <- sprintf('%sunbiased_hist_ising.pdf',image_res_folder)
h <- histogram_c_chains(mh_batch,1,k=0.1*K,K=K)
g <- plot_histogram(h, with_bar = T)
g <- g + labs(x = TeX('$\\beta$'),y='density')
ggsave(filename =name,width=7,height=5)




name <-  sprintf('%sising_mt.pdf',image_res_folder)

mt_df <- data.frame(mt_initial)
hist(mt_df$mt, nclass = 100)
dryhist <- hist(mt_df$mt, plot = FALSE, breaks=50)
breaks <- dryhist$breaks
mids <- dryhist$mids
breaks <- c(breaks)
mids <- c(mids)

meetingtimes <- mt_df$mt
counts <- hist(meetingtimes, breaks=50,plot = FALSE)$counts
counts <- counts / sum(counts)

dim(counts)
freqmeetingtimes <-counts
g <- qplot(x = mids, y = freqmeetingtimes, geom = "col") + xlab("meeting time") + ylab("estimated frequency")
g
# g <- g + geom_segment(aes(y = freqmeetingtimes, yend = freqmeetingtimes, x = mids-50, xend = mids))
ggsave(filename = name,width=7,height=5)





name <-  sprintf('%sising_Hkm_times.pdf',image_res_folder)

batch_times_vals <- batch_times_secs
mt_df <- data.frame(batch_times_vals)
hist(mt_df$batch_times_vals, nclass = 50)
dryhist <- hist(mt_df$batch_times_vals, plot = FALSE, breaks=50)
breaks <- dryhist$breaks
mids <- dryhist$mids
breaks <- c(breaks)
mids <- c(mids)

meetingtimes <- mt_df$batch_times_vals
counts <- hist(meetingtimes, breaks=50,plot = FALSE)$counts
counts <- counts / sum(counts)

dim(counts)
freqmeetingtimes <-counts
g <- qplot(x = mids, y = freqmeetingtimes, geom = "col") +
  xlab(TeX("clock time for $H_{k:m}$ (secs)")) +
  ylab("estimated frequency") +
  theme(legend.title=element_blank())


g
# g <- g + geom_segment(aes(y = freqmeetingtimes, yend = freqmeetingtimes, x = mids-50, xend = mids))
ggsave(filename = name,width=7,height=5)





# plot results
load(serial_datafile)
library(coda)

serial_lengths <- sapply(mh_serial_batches,function(x) dim(x$chain))
serial_length <- serial_lengths[1,1]
acfs <- sapply(mh_serial_batches,function(x) spectrum0.ar(x$chain[(0.1*serial_length):serial_length,])$spec)
serial_ineff <- mean(acfs)



# plot some results
chain <- mh_serial_batches[[1]]$chain
hist(chain[,1],probability=T,breaks=200)
abline(v=J_true,col='blue')

acf(chain[,1],300)

res_times_serial <- sapply(mh_serial_batches,function(x) as.numeric(x$res_time,units='secs'))

serial_ineff_timed <- mean(res_times_serial)*serial_ineff/(0.9*serial_length)
coupled_ineff_timed <-  var(ests_mat)*mean(batch_times_secs)



# compare inefficiencies
ineff1 <- coupled_ineff_timed/serial_ineff_timed

print(sprintf('Increase in inefficiency : beta=%.4f',ineff1[1]))


