# Miscellaneous manipulations of the files produced on the cluster

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
# datafiles <- list.files(pattern = "*originalsd[.]m100*")
# datafiles <- list.files(pattern = "*originalsd[.]m0[.]*")
# datafiles <- list.files(pattern = "cchains.m0.*")
datafiles <- list.files(pattern = "*ns[.]m1[0-9]*")
# datafiles <- list.files(pattern = "neuro.bpf[.]N*")
datafiles

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

head(df)

ns <- (df %>% group_by(iproc) %>% summarise(ns = n()))$ns
table(ns)

summary(df$duration)
hist(df$duration, prob = TRUE, nclass = 50)
#
g <- ggplot(df, aes(y = iproc, yend = iproc, x = starttime+0.01, xend = endtime-0.01, colour = factor(isample))) + geom_segment(lineend = "round") + geom_point()
g <- g + theme(legend.position = "none") + xlab("time") + ylab("processor")
g

library(dplyr)
summary((df %>% group_by(iproc) %>% summarise(meanmeet = mean(meetingtime)))$meanmeet)
g <- ggplot(df, aes(x = meetingtime, y = duration / 3600, colour = factor(iproc))) + geom_point() + theme(legend.position = "none")
g

# hist(df$meetingtime, nclass = 100)
dryhist <- hist(df$meetingtime, plot = FALSE, nclass = 100)
breaks <- dryhist$breaks
mids <- dryhist$mids
freqmeetingtimes <- sapply(1:length(datafiles), function(i){
  meetingtimes <- (df %>% filter(iproc == i))$meetingtime
  counts <- hist(meetingtimes, breaks = breaks, plot = FALSE)$counts
  counts <- counts / sum(counts)
  return(counts)
} )

freqmeetingtimes <- rowMeans(freqmeetingtimes)
g <- qplot(x = mids, y = freqmeetingtimes, geom = "point") + xlab("meeting time") + ylab("estimated frequency")
# g <- g + geom_segment(aes(y = freqmeetingtimes, yend = freqmeetingtimes, x = mids-50, xend = mids))
g <- g + geom_line()
g

g <- ggplot(df, aes(x = meetingtime, y = duration, label = iproc)) + geom_text() + theme(legend.position = "none")
g

# matplot(cchains_[[3]]$samples1[,2], type = "l")
br1 <- seq(from = 0.99, to = 1, length.out = 100)
prop1 <- rep(0, 99)
mid1 <- rep(0, 99)
br2 <- seq(from = 0.05, to = 0.2, length.out = 100)
prop2 <- rep(0, 99)
mid2 <- rep(0, 99)

for (i in 1:length(datafiles)){
  load(datafiles[i])
  his <- histogram_c_chains(cchains_, 1, 8e3, 1e4, breaks = br1)
  prop1 <- prop1 + his$proportions
  mid1 <- his$mids
  his <- histogram_c_chains(cchains_, 2, 8e3, 1e4, breaks = br2)
  prop2 <- prop2 + his$proportions
  mid2 <- his$mids
}


plot(mid1, prop1/100, type = "l")
plot(mid2, prop2/100, type = "l")

df[which(df$meetingtime <  5000),]
df %>% filter(iproc == 104)
load(datafiles[72])
head(cchains_[[1]]$samples2[,1])
matplot(cchains_[[1]]$samples1[4000:nrow(cchains_[[1]]$samples1),1], type = "l", ylim = c(0.99,1))
matplot(cchains_[[1]]$samples2[3999:nrow(cchains_[[1]]$samples2),1], add = TRUE, col = "red", type = "l")

matplot(cchains_[[1]]$samples1[1:nrow(cchains_[[1]]$samples2),1], col = "red", type = "l")
matplot(cchains_[[1]]$samples1[1:nrow(cchains_[[1]]$samples2),2], col = "red", type = "l")

plot(cchains_[[1]]$samples1[1:nrow(cchains_[[1]]$samples1),1], cchains_[[1]]$samples1[1:nrow(cchains_[[1]]$samples1),2], type = "l")


estim <- c()
for (i in 1:length(datafiles)){
  load(datafiles[i])
  # for (is in 1:length(cchains_)){
    estim <- c(estim, H_bar(cchains_[[1]], h = function(x) x[1], k = 1, K = 1e4))
  # }
}

summary(estim)
hist(estim, nclass = 100)
mean(estim)
median(estim)
