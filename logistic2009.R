##Numerical Example 1
## generate Figure 5 and paths : 5000 simulated paths of number of infected & death under 2009 H1N1 pandemic


library(COVID19)
library(demodelr)
library(matrixStats)
library(ggplot2)
library(forecast)
library(tseries)
library(RQuantLib)
library(quantmod)
library(dplyr)
library(tidyverse)
library(zoo)
library(Sim.DiffProc)
library(rstudioapi)

# Getting the path of your current open file
current_path = rstudioapi::getActiveDocumentContext()$path 
setwd(dirname(current_path ))
print( getwd() )



#####################Generate plots of number of infected and death under 2009 pandemic scenarios##########################################################
set.seed(123)
Total_infected<-833*10^6
Total_death<-0.4*10^6
g_I<-0.04

Bond_Price<-c()
P_C<-c()


y1_max<- Total_infected     #maximum infected
y2_max<-Total_death# maximum deaths
y1_min<-1 # initial number of infected
y2_min<-1 #initial number of death
numsim<-5000 # number of simulations
n<-1104 # number of days between July 07, 2017 -July 15, 2020 (including start date)
theta_I <-5000
theta_D<-2500

s0_I<-0.1  #volatility for infected
r0_I<-g_I #growth rate for infected
K0_I<-y1_max  #carrying capacity 
x0_I<-y1_min  #N0_I

s0_D<-0.1 #volatility for death
r0_D<-g_I #growth rate for death, we set this same to g_D
K0_D<-y2_max
x0_D<-y2_min     #N0_D

dt<-1   #time step
m<-n  #number of steps
M<-numsim #number of paths(simulations)
rho<-0.5 #correlation coefficient between death and infected
#######################################################################

fun_st_logistic<- function(dt,m,rho,M){
fx<-expression(r0_I*x*(1-x/(K0_I )),r0_D*y*(1-y/(K0_D)))
gx<-expression(s0_I*x*(1-x/(K0_I )),s0_D*y*(1-y/(K0_D)))
x0<-c(x0_I,x0_D)
t0<-0
Sigma<-matrix(c(1,rho,rho,1),nrow=2)
set.seed(123)
res<-snssde2d(N=m,M=M,x0=x0,t0=t0,Dt=dt,drift=fx,diffusion=gx,corr=Sigma)

 return(res )
}

res<-fun_st_logistic(dt,m,rho,M)


simulated_I<-matrix(unlist(res$X),ncol=numsim,byrow=F)
simulated_D<-matrix(unlist(res$Y),ncol=numsim,byrow=F)


#plot simulated paths#################################################

gbm<-simulated_I


gbm_dfI <- as.data.frame(gbm) %>%
  mutate(ix = 1:nrow(gbm)) %>%
  pivot_longer(-ix, names_to = 'sim', values_to = 'price')

#par(mfrow=c(1,2))

p3<-ggplot(data=gbm_dfI,aes(x=ix, y=price, color=sim)) +
  geom_line(color="black") +
  #geom_line(data = filter(gbm_df, sim == "y1"), size = 1, color="black")+
  theme(legend.position = 'none')+      
  labs(x="Days", y="Number of infected")
#geom_hline(yintercept = theta_I, linetype="dashed", color = "red", size=1)

x11();p3




gbm<-simulated_D

gbm_dfD <- as.data.frame(gbm) %>%
  mutate(ix = 1:nrow(gbm)) %>%
  pivot_longer(-ix, names_to = 'sim', values_to = 'price')

p4<-ggplot(data=gbm_dfD, aes(x=ix, y=price, color=sim)) +
  geom_line(color="black") +
  #geom_line(data = filter(gbm_df, sim == "y2"), size = 1, color="black")+
  theme(legend.position = 'none')+      
  labs(x="Days", y="Number of death")
#geom_hline(yintercept = theta_D, linetype="dashed", color = "red", size=1)

x11();p4






