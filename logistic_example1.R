#Numerical example 1
# need numerical example 1 coupons
#Calculate Table 5 values: World Bank issued pandemic bond, priced under each historic pandemic if required yield
#rate is 8.6734%


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

graphics.off()  # clear all graphs
rm(list = ls()) # remove all files from your workspace

# change directory to where the script located

my_d <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(my_d)


###############################################################################
set.seed(123)
Total_infected<-c(2000*10^6,1000*10^6,1000*10^6,833*10^6)
Total_death<-c(50*10^6,4*10^6,4*10^6,0.4*10^6)
g_I<-c(0.049,0.125,0.075,0.04)

Bond_Price<-c()
P_C<-c()

for(j in 1:4){
  print("Simulating scenario",j)

y1_max<- Total_infected[j]     #maximum infected
y2_max<-Total_death[j]# maximum deaths
y1_min<-1 # initial number of infected
y2_min<-1 #initial number of death
numsim<-5000 # number of simulations
n<-1104 # number of days between July 07, 2017 -July 15, 2020 (including start date)
theta_I <-5000
theta_D<-2500

s0_I<-0.1  #volatility for infected
r0_I<-g_I[j] #growth rate for infected
K0_I<-y1_max  #carrying capacity 
x0_I<-y1_min  #N0_I

s0_D<-0.1 #volatility for death
r0_D<-g_I[j] #growth rate for death
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



###########################We calculate seven day average for the data set ################
###First find the real increment from previous day to today ##############

Increment_I<-colDiffs(simulated_I, lag=1)
Increment_I<-pmax(Increment_I,0)
nn<-nrow(Increment_I)


index<-1:nn

sim_avg7_I<- matrix(ncol = numsim, nrow = length(7:nn))
sim_GR_I<-matrix(ncol = numsim, nrow = length(7:nn))

for(i in 1:numsim){
  obj1<-zoo(Increment_I[,i],index)
  seven_avg_I<-rollmean(obj1,7,align='right',fill=0)
  seven_sd_I<-rollapply(data = obj1,width=7,FUN=sd,align='right',fill=0)
  
  y11<-as.data.frame(seven_avg_I[7:nn])
  y11<-y11[["seven_avg_I[7:nn]"]]
  sim_avg7_I[,i]<-y11
  
  y33<-as.data.frame(seven_sd_I[7:nn])
  y33<-y33[["seven_sd_I[7:nn]"]]
  sim_GR_I[,i]<-(y11-(1.533*y33))
  
}

Increment_D<-colDiffs(simulated_D, lag=1)
Increment_D<-pmax(Increment_D,0)
sim_avg7_D<- matrix(ncol = numsim, nrow = length(7:nn))


for(i in 1:numsim){
  obj1<-zoo(Increment_D[,i],index)
  seven_avg_D<-rollmean(obj1,7,align='right',fill=0)
  
  
  y22<-as.data.frame(seven_avg_D[7:nn])
  y22<-y22[["seven_avg_D[7:nn]"]]
  sim_avg7_D[,i]<-y22
  
  
}



##################################################################


trigger<-0

set_I<-0
set_D<-0
set_GR<-0
common<-0
common_index<-c()
common_time<-c()

for(i in 1:numsim){
  set_I<-which(sim_avg7_I[,i]>theta_I)
  set_D<-which(sim_avg7_D[,i]>theta_D)
  set_GR<-which(sim_GR_I[,i]>0)
  common<-Reduce(intersect, list(set_I,set_D,set_GR)) 
  if(length(common)>0){
    trigger<-trigger+1
    common_index<-append(common_index,i) #in which simulation this happen
    common_time<-append(common_time,common[1]) #at what time this happen, then add first time all criteria met
  }
  common<-NULL
}


#read predicted coupons. last elements is coupon+redemption amount
coupon_data<-readRDS(file = "coupon.rds") 
#A function to calculate present value
pv <- function(y, cf, t){ sum(cf/(1+y)^t)}  

#y1 is the yield rate calculated early assuming no pandemic. 
y1<-0.08673402

BondValue_partial<-c() 

for(i in 1:length(common_index)){
  trigger_time<-common_time[i]+6 #we add 6, since these are obtained from 7 day average, which set 7->1. So add 6 to bring back to original
  coupondf<-subset(coupon_data, coupon_data$Time<=trigger_time)
  cf<-as.numeric(coupondf$CF)
  
  t<-coupondf$Time/360
  presentV<-pv(y1,cf,t)
  BondValue_partial<-append(BondValue_partial,presentV)
}

#E[BondPrice|Trigger]
price1<-mean(BondValue_partial)
#E[BondPrice|No Trigger]
price2<-225 #this is the price advertised by World Bank/paid by investors

#P(H) probability of pandemic 
prob1<-0.1393
#P(C|H) probability that trigger activated given pandemic
print(trigger)
prob2<-(trigger/numsim)
# P(C) probability of trigger activated
prob3<-prob1*prob2
#P(C|H') probability that trigger activated when no pandemic is zero.
#P(C') probability that trigger not activated
prob4<-1-prob3

# The price of the bond
#E[Bond Price|Trigger]P (Trigger)+E[Bond Price| No Trigger]P (No Trigger)

price3<-(price1*prob3)+(price2*prob4)
Bond_Price<-append(Bond_Price,price3)
P_C<-append(P_C,prob3)
}


round(Bond_Price,4)
round(P_C,4)

mean(Bond_Price) # this is the price investor should pay
