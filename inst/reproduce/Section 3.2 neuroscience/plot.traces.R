# plots traces coming from the longest runs on the cluster;
# this creates plots for Figure 6 in the manuscript
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

source("binomialmodel.R")

# load 2d grid of likelihood estimates
load("csmc2dgrid.RData")
thetas.df$lls <- lls
thetas.df$dpr <- apply(thetas.df, 1, function(v) dprior(as.numeric(v[1:2])))
thetas.df$dpost <- thetas.df$lls + thetas.df$dpr

g <- ggplot(thetas.df, aes(x = Var1, y = Var2, z = dpost, fill = dpost))  # + geom_raster()
g <- g + theme(legend.position = "none") + geom_contour(data = thetas.df, bins = 100, alpha = 0.3, colour = "black")
g <- g + xlab(expression(a)) + ylab(expression(sigma[X]^2))
g


# load long traces obtained with small ("original") standard deviation of random walk proposals
pathtoodysseyresults <- "outputfromcluster"
datafiles <- list.files(path = pathtoodysseyresults, pattern = "*sd[.]m[0-9]*")

df <- data.frame()
for (ifile in 1:length(datafiles)){
  load(file = paste0(pathtoodysseyresults, datafiles[ifile]))
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

head(df)
summary(df$meetingtime)
df[which.max(df$meetingtime),]
load(paste0(pathtoodysseyresults, datafiles[270]))
isample <- 1
trace1.df <- data.frame(x = cchains_[[isample]]$samples1[1:nrow(cchains_[[isample]]$samples1),1], y = cchains_[[isample]]$samples1[1:nrow(cchains_[[isample]]$samples1),2])
trace2.df <- data.frame(x = cchains_[[isample]]$samples2[1:nrow(cchains_[[isample]]$samples2),1], y = cchains_[[isample]]$samples2[1:nrow(cchains_[[isample]]$samples2),2])
gtrace <- g + geom_path(data=trace1.df, aes(x = x, y = y, z = NULL, fill = NULL))
gtrace <- gtrace + geom_path(data=trace2.df, aes(x = x, y = y, z = NULL, fill = NULL))
gtrace
# ggsave(filename = "neuro.longesttrace.originalsd.png", plot = gtrace, width = 7, height = 5)

head(trace1.df)
g1 <-  ggplot(trace1.df, aes(x = 1:nrow(trace1.df), y = x)) + geom_line() + xlab("iteration") + ylab(expression(a))
g2 <- ggplot(trace1.df, aes(x = 1:nrow(trace1.df), y = y)) + geom_line()+ xlab("iteration") + ylab(expression(sigma[X]^2))
library(gridExtra)
g12 <- grid.arrange(g1, g2, nrow = 2)
g12
# ggsave(filename = "neuro.longtraces.pdf", plot = g12, width = 7, height = 5)


mean(abs(diff(trace1.df$x)) > 1e-5)
mean(abs(diff(trace2.df$x)) > 1e-5)

datafiles <- list.files(path = pathtoodysseyresults, pattern = "*cchains[.]m[0-9]*")
df <- data.frame()
for (ifile in 1:length(datafiles)){
  load(file = paste0(pathtoodysseyresults, datafiles[ifile]))
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

head(df)
summary(df$meetingtime)
df[which(df$meetingtime >  900),]
load(paste0(pathtoodysseyresults, datafiles[37]))
isample <- 33
trace1.df <- data.frame(x = cchains_[[isample]]$samples1[1:nrow(cchains_[[isample]]$samples1),1], y = cchains_[[isample]]$samples1[1:nrow(cchains_[[isample]]$samples1),2])
trace2.df <- data.frame(x = cchains_[[isample]]$samples2[1:nrow(cchains_[[isample]]$samples2),1], y = cchains_[[isample]]$samples2[1:nrow(cchains_[[isample]]$samples2),2])
gtrace <- g + geom_path(data=trace1.df, aes(x = x, y = y, z = NULL, fill = NULL))
gtrace <- gtrace + geom_path(data=trace2.df, aes(x = x, y = y, z = NULL, fill = NULL))
gtrace
# ggsave(filename = "neuro.longesttrace.largersd.png", plot = gtrace, width = 7, height = 5)
# mean(abs(diff(trace1.df$x)) > 1e-5)
# mean(abs(diff(trace2.df$x)) > 1e-5)

