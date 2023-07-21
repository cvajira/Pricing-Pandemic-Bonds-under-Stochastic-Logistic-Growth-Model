## This code is used to calculate probability of pandemic for numerical example 1, 2 as to create Figure 4
## both weibull and exponential implemented


library(Rfast2)
library(ggplot2)
library(fitdistrplus)
library(EnvStats)

graphics.off()  # clear all graphs
#rm(list = ls()) # remove all files from your workspace

# change directory to where the script located

my_d <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(my_d)
      

set.seed(123)
data<-c(3, 49, 7, 42, 3, 3, 53, 10, 19, 39, 11, 9, 26,3, 5,3)
censored<-c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,0)

fit<-censweibull.mle(data,censored)

#tau=shape
#theta=scale

qqPlot(data, distribution = "weibull", param.list = list(shape = 1/fit$param[[2]], scale = fit$param[[1]]),add.line=TRUE,pch=19, line.lwd=2)

ks.test(data, "pweibull", shape = 1/fit$param[[2]], scale = fit$param[[1]])

P_H<-(pweibull(6,shape = 1/fit$param[[2]], scale = fit$param[[1]])-pweibull(3,shape = 1/fit$param[[2]], scale = fit$param[[1]]))/(1-pweibull(3,shape = 1/fit$param[[2]], scale = fit$param[[1]]))
P_H_c<-1-P_H
P_H
P_H_c

###############For second example

set.seed(123)
data<-c(3, 49, 7, 42, 3, 3, 53, 10, 19, 39, 11, 9, 26,3, 5,5,4)
censored<-c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,1,0)

fit<-censweibull.mle(data,censored)

#tau=shape
#theta=scale

qqPlot(data, distribution = "weibull", param.list = list(shape = 1/fit$param[[2]], scale = fit$param[[1]]),add.line=TRUE,pch=19, line.lwd=2)

ks.test(data, "pweibull", shape = 1/fit$param[[2]], scale = fit$param[[1]])

P_H<-(pweibull(6,shape = 1/fit$param[[2]], scale = fit$param[[1]])-pweibull(3,shape = 1/fit$param[[2]], scale = fit$param[[1]]))/(1-pweibull(3,shape = 1/fit$param[[2]], scale = fit$param[[1]]))
P_H_c<-1-P_H
P_H
P_H_c

############################Fit an exponential distirbution using mle. Read page 261 in loss models(4th ed)

set.seed(123)
data<-c(3, 49, 7, 42, 3, 3, 53, 10, 19, 39, 11, 9, 26,3, 5,3)
censored<-c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,1,0)
theta_mle<-sum(data)/sum(censored)

qqPlot(data, distribution = "exp", param.list = list(rate=1/theta_mle),add.line=TRUE,pch=19, line.lwd=2)
ks.test(data, pexp, rate=1/theta_mle)

P_H<-(pexp(6,rate=1/theta_mle)-pexp(3,rate=1/theta_mle))/(1-pexp(3,rate=1/theta_mle))
P_H_c<-1-P_H
P_H
P_H_c

