# creates Figure 4 of the manuscript

library(debiasedpmcmc)
rm(list = ls())
set.seed(1)
library(ggplot2)
library(dplyr)
library(doParallel)
library(doRNG)
registerDoParallel(cores = 4)
setmytheme()

# set path
# setwd("")
# load model and relevant functions
source("binomialmodel.R")

## Plot data
g <- qplot(x = 1:3000, y = observations[,1], geom = "line") + xlab("time (ms)") + ylab("observations")
g
# ggsave(filename = "~/Dropbox/UPMCMC/experiments/neuro.data.png", plot = g, width = 7, height = 5)


## Plot likelihood
## (needs to have previously run "csmc.2dgrid.R")
load("csmc2dgrid.RData")
# any(!is.finite(lls))
# summary(as.numeric(lls))
# as.numeric(lls)
thetas.df$lls <- lls
thetas.df$dpr <- apply(thetas.df, 1, function(v) dprior(as.numeric(v[1:2])))
thetas.df$dpost <- thetas.df$lls + thetas.df$dpr

g <- ggplot(thetas.df, aes(x = Var1, y = Var2, z = dpost, fill = dpost))  + geom_raster()
g <- g + scale_fill_gradientn(colours = c("navy", "white", "red"),
                              values = scales::rescale(c(-30000, -7000, -5000)))
g <- g + theme(legend.position = "none") + geom_contour(data = thetas.df, bins = 100, alpha = 0.6, colour = "black")
g <- g + xlab(expression(a)) + ylab(expression(sigma[X]^2))
g <- g + geom_point(data=NULL, aes(x = thetas.df[which.max(lls),1], y = thetas.df[which.max(lls),2]), colour = "black")
g
# ggsave(filename = "~/Dropbox/UPMCMC/experiments/neuro.likelihood.png", plot = g, width = 7, height = 5)


## Extra...
## We can also look at it in 3d using rgl's function persp3d
# z <- matrix(thetas.df$dpost, nrow = length(theta1s), byrow = FALSE)
# library(rgl)
# persp3d(theta1s, theta2s, z, xlab = "alpha", ylab = "sigma2", col = "red")
