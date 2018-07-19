# plots meeting times obtained by running coupled chains (on a large cluster)
# leads to figures 5 and 7
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

## load meeting times associated with BPF
datafiles <- list.files(pattern = "bpf[.]")
df.bpf <- load_df(datafiles)
datafiles2 <- list.files(pattern = "bpf2[.]")
df.bpf2 <- load_df(datafiles2)
df.bpf2$iproc <- df.bpf2$iproc + 200
df.bpf <- rbind(df.bpf, df.bpf2)

# number of samples per processors
ns.bpf <- (df.bpf %>% group_by(iproc) %>% summarise(ns = n()))$ns
table(ns.bpf)
head(df.bpf)
g <- ggplot(df.bpf, aes(x = meetingtime, y = duration/3600, colour = factor(iproc))) + geom_point() + theme(legend.position = "none")
g <- g + xlab("meeting time") + ylab("duration (hours)")
g
# ggsave(filename = "neuro.bpf.durationvsmeeting.png", plot = g, width = 7, height = 5)
#
summary(df.bpf$meetingtime)

hist(df.bpf$meetingtime, nclass = 100)
dryhist <- hist(df.bpf$meetingtime, plot = FALSE, nclass = 100)
breaks <- dryhist$breaks
mids <- dryhist$mids
breaks <- c(breaks, 9500)
mids <- c(mids, 9450)
freqmeetingtimes <- sapply(1:length(datafiles), function(i){
  meetingtimes <- (df.bpf %>% filter(iproc == i))$meetingtime
  counts <- hist(meetingtimes, breaks = breaks, plot = FALSE)$counts
  counts <- counts / sum(counts)
  return(counts)
} )

dim(freqmeetingtimes)
freqmeetingtimes <- rowMeans(freqmeetingtimes)
g <- qplot(x = mids, y = freqmeetingtimes, geom = "col") + xlab("meeting time") + ylab("frequency")
g
# ggsave(filename = "neuro.bpf.meetingfrequencies.png", plot = g, width = 7, height = 5)

## load meeting times associated with cSMC
datafiles <- list.files(pattern = "*cchains[.]m[0-9]*")
df.csmc <- load_df(datafiles)
# number of samples per processors
ns.csmc <- (df.csmc %>% group_by(iproc) %>% summarise(ns = n()))$ns
table(ns.csmc)
df.csmc$iproc %>% unique
#
head(df.csmc)
g <- ggplot(df.csmc, aes(x = meetingtime, y = duration/3600, colour = factor(iproc))) + geom_point() + theme(legend.position = "none")
g <- g + xlab("meeting time") + ylab("duration (hours)")
g
# ggsave(filename = "neuro.csmc.durationvsmeeting.png", plot = g, width = 7, height = 5)
#

hist(df.csmc$meetingtime, nclass = 100)
dryhist <- hist(df.csmc$meetingtime, plot = FALSE, nclass = 50)
breaks <- dryhist$breaks
mids <- dryhist$mids
freqmeetingtimes <- sapply(1:length(datafiles), function(i){
  meetingtimes <- (df.csmc %>% filter(iproc == i))$meetingtime
  counts <- hist(meetingtimes, breaks = breaks, plot = FALSE)$counts
  counts <- counts / sum(counts)
  return(counts)
} )

dim(freqmeetingtimes)
freqmeetingtimes <- rowMeans(freqmeetingtimes)
g <- qplot(x = mids, y = freqmeetingtimes, geom = "col") + xlab("meeting time") + ylab("frequency")
g
# ggsave(filename = "neuro.csmc.meetingfrequencies.png", plot = g, width = 7, height = 5)


## load meeting times associated with cSMC and small sd
datafiles <- list.files(pattern = "*sd[.]m[0-9]*")
df.csmc <- load_df(datafiles)
# number of samples per processors
ns.csmc <- (df.csmc %>% group_by(iproc) %>% summarise(ns = n()))$ns
table(ns.csmc)
df.csmc$iproc %>% unique
#
head(df.csmc)
g <- ggplot(df.csmc, aes(x = meetingtime, y = duration/3600, colour = factor(iproc))) + geom_point() + theme(legend.position = "none")
g <- g + xlab("meeting time") + ylab("duration (hours)")
g
# ggsave(filename = "neuro.csmc.durationvsmeeting.png", plot = g, width = 7, height = 5)
#
summary(df.csmc$meetingtime)
hist(df.csmc$meetingtime, nclass = 100)
dryhist <- hist(df.csmc$meetingtime, plot = FALSE, nclass = 100)
breaks <- dryhist$breaks
mids <- dryhist$mids
freqmeetingtimes <- sapply(1:length(datafiles), function(i){
  meetingtimes <- (df.csmc %>% filter(iproc == i))$meetingtime
  counts <- hist(meetingtimes, breaks = breaks, plot = FALSE)$counts
  counts <- counts / sum(counts)
  return(counts)
} )

dim(freqmeetingtimes)
freqmeetingtimes <- rowMeans(freqmeetingtimes)

g <- qplot(x = mids, y = freqmeetingtimes, geom = "col") + xlab("meeting time") + ylab("frequency")
g
# ggsave(filename = "neuro.csmc.originalsd.meetingfrequencies.png", plot = g, width = 7, height = 5)
