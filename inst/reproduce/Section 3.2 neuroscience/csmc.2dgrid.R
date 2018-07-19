# runs controlled SMC on a grid of parameter values, in order to create
# contour plots of the likelihood as in Figure 4 of the manuscript

# remove all objects from R environment
rm(list = ls())
# load package
library(debiasedpmcmc)
library(doRNG)
library(doParallel)
registerDoParallel(cores = 10)
### fix the random seed
set.seed(17)
# set path
# setwd("")
source("binomialmodel.R")
safell <- function(theta){
  ll <- try(csmc(nparticles, as.numeric(theta), observations, niter)$ll)
  if (inherits(ll, "try-error")){
    return(-Inf)
  } else {
    return(ll)
  }
}


nthetas1 <- 500
theta1s <- seq(from = 0.01, to = 0.999, length.out = nthetas1)
nthetas2 <- 500
theta2s <- seq(from = 0.01, to = 4, length.out = nthetas2)
thetas.df <- expand.grid(theta1s, theta2s)
## RUN
# lls <- foreach(itheta = 1:(dim(thetas.df)[1]), .combine = c) %dorng% {
#   safell(thetas.df[itheta,])
# }
## SAVE
# save(theta1s, theta2s, thetas.df, lls, file = "csmc2dgrid.RData")
## LOAD
load("csmc2dgrid.RData")
# any(!is.finite(lls))
# summary(as.numeric(lls))
# as.numeric(lls)

thetas.df$lls <- lls
thetas.df$dpr <- apply(thetas.df, 1, function(v) dprior(as.numeric(v[1:2])))
thetas.df$dpost <- thetas.df$lls + thetas.df$dpr

library(ggplot2)
# plot contours
g <- ggplot(thetas.df, aes(x = Var1, y = Var2, z = dpost, fill = dpost))  + geom_raster()
g <- g + scale_fill_gradient2(low = "navy", mid = "white", high = "red", midpoint = -10000,
                              limits = c(-20000,-5000))
g + theme(legend.position = "none") + geom_contour(data = thetas.df %>% filter(dpost > -20000), bins = 20, alpha = 0.6, colour = "black")

g <- ggplot(thetas.df, aes(x = Var1, y = Var2, z = dpr)) + geom_contour(bins = 200, alpha = 0.6, colour = "black")
g <- g + geom_point(data=NULL, aes(x = thetas.df[which.max(lls),1], y = thetas.df[which.max(lls),2]), colour = "black")
g
