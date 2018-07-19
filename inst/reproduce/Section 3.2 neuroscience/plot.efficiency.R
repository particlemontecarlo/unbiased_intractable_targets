# processes files from runs on the cluster,
# and does the calculation of asymptotic inefficiencies
# so that numbers can be reported as in Section 3.2.4
# also creates Figure 8

library(debiasedpmcmc)
rm(list = ls())
set.seed(1)
library(ggplot2)
library(dplyr)
library(doParallel)
library(doRNG)
registerDoParallel(cores = 4)
setmytheme()
setwd(dir = "outputfromcluster/")
load_df <- function(datafiles){
  df <- data.frame()
  for (ifile in 1:length(datafiles)){
    load(file = datafiles[ifile])
    if (length(cchains_) == 0){
      df <- rbind(df, data.frame(iproc = ifile, isample = 1,
                                 meetingtime = NA,
                                 starttime = NA,
                                 endtime = NA,
                                 duration = NA))

    } else {
      df_ <- data.frame()
      for (isample in 1:nsamples_){
        if (isample == 1){
          starttime <- 0
        } else {
          starttime <- starttime + durations_[isample-1]
        }
        meetingtime <- cchains_[[isample]]$meetingtime
        df_ <- rbind(df_, data.frame(iproc = ifile, isample = isample,
                                     meetingtime = meetingtime,
                                     starttime = starttime,
                                     endtime = starttime + durations_[isample],
                                     duration = durations_[isample]))
      }
      df <- rbind(df, df_)
    }
  }
  return(df)
}

## load meeting times associated with cSMC
datafiles <- list.files(pattern = "cchains.m10000.*")
datafiles
df.csmc <- load_df(datafiles)

ns <- (df.csmc %>% group_by(iproc) %>% summarise(ns = n()))$ns
table(ns)
df.csmc$starttime <- df.csmc$starttime/3600
df.csmc$duration <- df.csmc$duration/3600
df.csmc$endtime <- df.csmc$endtime/3600
g <- ggplot(df.csmc, aes(y = iproc, yend = iproc, x = starttime+0.01, xend = endtime-0.01, colour = factor(isample))) + geom_segment(lineend = "round") + geom_point()
g <- g + theme(legend.position = "none") + xlab("time (hours)") + ylab("processor") + geom_vline(xintercept = 23, linetype = 2)
g
# ggsave(filename = "neuro.csmc.chronology.png", plot = g, width = 7, height = 5)

## Histogram of parameters
br1 <- seq(from = 0.992, to = 1, length.out = 50)
prop1 <- rep(0, length(br1)-1)
mid1 <- rep(0, length(br1)-1)
br2 <- seq(from = 0.05, to = 0.2, length.out = 50)
prop2 <- rep(0, length(br2)-1)
mid2 <- rep(0, length(br2)-1)

for (i in 1:length(datafiles)){
  load(datafiles[i])
  his <- histogram_c_chains(cchains_, 1, 1e3, 1e4, breaks = br1)
  prop1 <- prop1 + his$proportions
  mid1 <- his$mids
  his <- histogram_c_chains(cchains_, 2, 1e3, 1e4, breaks = br2)
  prop2 <- prop2 + his$proportions
  mid2 <- his$mids
}

load("pmmh.csmc.T250000.N128.niter3.RData")
burnin <- 1e4
pmmh.csmc.chain <- chain[burnin:nmcmc,]
mcmc.counts1 <- hist(pmmh.csmc.chain[,1], breaks = br1, prob = TRUE)$counts
mcmc.proportions1 <- mcmc.counts1 / sum(mcmc.counts1)
mcmc.counts2 <- hist(pmmh.csmc.chain[,2], breaks = br2, prob = TRUE)$counts
mcmc.proportions2 <- mcmc.counts2 / sum(mcmc.counts2)

g1 <- qplot(mid1, prop1/100, geom = "col") + ylim(0, 0.1) + xlab(expression(a)) + ylab("frequency")
g1 <- g1 + geom_line(aes(y = mcmc.proportions1), colour = 'red', linetype = 1)
g1

g2 <- qplot(mid2, prop2/100, geom = "col") + ylim(0, 0.2) + xlab(expression(sigma[X]^2)) + ylab("frequency")
g2 <- g2 + geom_line(aes(y = mcmc.proportions2), colour = 'red', linetype = 1)
g2

library(gridExtra)
g12 <- grid.arrange(g1, g2, nrow = 2)
# ggsave(filename = "neuro.csmc.estimparameters.png", plot = g12, width = 7, height = 5)


### inefficiency
load("pmmh.csmc.T250000.N128.niter3.RData")
burnin <- 1e4
elapsed / 3600
elapsed / nmcmc # time per mcmc step
testchain <- apply(chain[burnin:nmcmc,], MARGIN = 1, function(x) sum(x) + sum(x^2))
mean(testchain)
library(coda)
var(testchain)
# inefficiency parametrized by time
elapsed * spectrum0(testchain)$spec / (nmcmc - burnin)
# inefficiency parametrized by number of iterations
nmcmc * spectrum0(testchain)$spec / (nmcmc - burnin)

estimators <- rep(0, length(datafiles))
durations <- rep(0, length(datafiles))
costs <- rep(0, length(datafiles))
for (ifile in 1:length(datafiles)){
  load(datafiles[ifile])
  costs[ifile] <- sum(sapply(1:length(cchains_), function(i) 2*cchains_[[i]]$meetingtime + max(1, cchains_[[i]]$iteration-cchains_[[i]]$meetingtime+1)))
  durations[ifile] <- sum(durations_)
  estimators[ifile] <- mean(sapply(1:length(cchains_), function(i) H_bar(cchains_[[i]], h = function(x) sum(x) + sum(x^2), k = 1e3, K = 1e4)))
}

mean(costs) * var(estimators)
23*3600 * var(estimators)
mean(estimators)
# inefficiency
# parametrized by time
23 * 3600 * var(estimators)
# parametrized by number of iterations
mean(costs) * var(estimators)

