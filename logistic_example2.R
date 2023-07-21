#Numerical example 2
# Require coupon for example 2
# Figure 8, table 8 and table 9 was created.


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

USdata<-covid19( country = "US", level = 1, start = "2020-03-01", end = "2022-12-31")
USdata[is.na(USdata)]<-0

#extract only cumulative "confirmed" and "death".
covid_data<-USdata[,c("date","confirmed","deaths" )]
#write.csv(covid_data,"covid_data.csv")


y1<-USdata$confirmed
y2<-USdata$deaths
#y1<-(y1-min(y1))/max(y1)



n<-length(y1)
days<-0:(n-1)
x1<-days
x2<-days

#First we calculate daily growth rate for infected(confirmed) ####################
#daily growth (y_new - y_old/y_old)

daily_infected<-diff(y1,lag=1)
daily_growth_I<-daily_infected/head(y1,-1)
dat1<-as.data.frame(cbind(1:length(daily_growth_I),daily_growth_I))
names(dat1)<-c("Date","Growth.rate")

p1 <- ggplot()+
  geom_line(data=dat1, aes(x =Date , y = Growth.rate),color = 'blue', size = 0.75)+
  labs(title = '',
       x = 'Day',
       y = 'Daily Growth rate: Infections',
       caption = '') + 
  theme_bw(base_family = "TT Times New Roman") 


x11(); p1;ggsave(filename ='growth_rate_I.png', plot = p1)



######same for death
daily_death<-diff(y2,lag=1)
daily_growth_D<-daily_death/head(y1,-1)
dat2<-as.data.frame(cbind(1:length(daily_growth_D),daily_growth_D))
names(dat2)<-c("Date","Growth.rate")

p2 <- ggplot()+
  geom_line(data=dat2, aes(x =Date , y = Growth.rate),color = 'blue', size = 0.75)+
  labs(title = '',
       x = 'Day',
       y = 'Daily Growth rate: Death',
       caption = '') + 
  theme_bw(base_family = "TT Times New Roman") 


x11(); p2;ggsave(filename ='growth_rate_D.png', plot = p2)

p3<-ggarrange(p1, p2,ncol = 2, nrow = 1)

x11(); p3;ggsave(filename ='growth_rate_combine.png', plot = p3)


###################################################Calculating g_I,g_D#####################
gI_1<-round(max(daily_growth_I),4)
gD_1<-round(max(daily_growth_D),4)

gI_2<-round(mean(daily_growth_I[daily_growth_I>0]),4)
gD_2<-round(mean(daily_growth_D[daily_growth_D>0]),4)


gI_3<-round(mean(daily_growth_I[1:50]),4)
gD_3<-round(mean(daily_growth_D[1:50]),4)

###############################Calculating sigma_I and sigmaD###############
sI_1<-round(sd(daily_growth_I),4)
sD_1<-round(sd(daily_growth_D),4)

sI_2<-round(sd(daily_growth_I[daily_growth_I>0]),4)
sD_2<-round(sd(daily_growth_D[daily_growth_D>0]),4)


sI_3<-round(sd(daily_growth_I[1:50]),4)
sD_3<-round(sd(daily_growth_D[1:50]),4)

###########################################################correlation coefficient #########
#cor(y1,y2) #not the correct way. this would imply there is a linear relationship and anyone get sick going to die, since rate is high

###############################################################################
set.seed(123)
y1<-USdata$confirmed
y2<-USdata$deaths
Total_infected<-max(y1)
Total_death<-max(y2)
g_I<-c(gI_1,gI_2,gI_3)
g_D<-c(gD_1,gD_2,gD_3)
s_I<-c(sI_1,sI_2,sI_3)
s_D<-c(sD_1,sD_2,sD_3)

Bond_Price<-c()
P_C<-c()

for(j in 1:3){
  print(paste0("Simulating scenario ",j))
  
  y1_max<- Total_infected     #maximum infected
  y2_max<-Total_death# maximum deaths
  y1_min<-1 # initial number of infected
  y2_min<-1 #initial number of death
  numsim<-5000 # number of simulations
  n<-1104 # number of days between January 06, 2023 -January 15, 2026 (including start date)
  theta_I <-5000
  theta_D<-2500
  
  s0_I<-s_I[j]  #volatility for infected
  r0_I<-g_I[j] #growth rate for infected
  K0_I<-y1_max  #carrying capacity 
  x0_I<-y1_min  #N0_I
  
  s0_D<-s_D[j] #volatility for death
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
  coupon_data<-readRDS(file = "coupon_example2.rds") 
  #A function to calculate present value
  # i discount rate, cf cash flows, t time
  
  pv <- function(i, cf, t){ sum(cf*(i)^t)}
  

  
  BondValue_partial<-c() 
  
  for(i in 1:length(common_index)){
    trigger_time<-common_time[i]+6 #we add 6, since these are obtained from 7 day average, which set 7->1. So add 6 to bring back to original
    coupondf<-subset(coupon_data, coupon_data$Time<=trigger_time)
    cf<-as.numeric(coupondf$CF)
    y_rate<-as.numeric(coupondf$Disc.rate)
    t<-coupondf$Time/360
    presentV<-pv(y_rate,cf,t)
    BondValue_partial<-append(BondValue_partial,presentV)
  }
  
  #E[BondPrice|Trigger]
  price1<-mean(BondValue_partial)
  #E[BondPrice|No Trigger]
  price2<-225 #this is the price advertised by World Bank/paid by investors
  
  #P(H) probability of pandemic with covid-19
  prob1<-0.1468
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

#############################################################
#############################################################
#############################################################
##############################################################
set.seed(123)
y1<-USdata$confirmed
y2<-USdata$deaths
Total_infected<-max(y1)
Total_death<-max(y2)
g_I<-seq(0.1,0.4,0.1)
g_D<-seq(0.1,0.4,0.1)
s_I<-seq(0.02,0.2,0.02)
s_D<-seq(0.02,0.2,0.02)
rho<-seq(0.1,0.9,0.1)
d1<-expand.grid(g_I,g_D,s_I,s_D,rho)
len_d1<-nrow(d1)
Bond_Price<-c()
P_C<-c()

for(j in 1:len_d1){
  print(paste0("Simulating scenario ",j," out of ", len_d1))
  
  y1_max<- Total_infected     #maximum infected
  y2_max<-Total_death# maximum deaths
  y1_min<-1 # initial number of infected
  y2_min<-1 #initial number of death
  numsim<-100 # number of simulations
  n<-1104 # number of days between January 06, 2023 -January 15, 2026 (including start date)
  theta_I <-5000
  theta_D<-2500
  
  s0_I<-d1[j,3]  #volatility for infected
  r0_I<-d1[j,1] #growth rate for infected
  K0_I<-y1_max  #carrying capacity 
  x0_I<-y1_min  #N0_I
  
  s0_D<-d1[j,4] #volatility for death
  r0_D<-d1[j,2] #growth rate for death
  K0_D<-y2_max
  x0_D<-y2_min     #N0_D
  
  dt<-1   #time step
  m<-n  #number of steps
  M<-numsim #number of paths(simulations)
  rho<-d1[j,5] #correlation coefficient between death and infected
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
  coupon_data<-readRDS(file = "coupon_example2.rds") 
  #A function to calculate present value
  # i discount rate, cf cash flows, t time
  
  pv <- function(i, cf, t){ sum(cf*(i)^t)}
  
  
  
  BondValue_partial<-c() 
  
  for(i in 1:length(common_index)){
    trigger_time<-common_time[i]+6 #we add 6, since these are obtained from 7 day average, which set 7->1. So add 6 to bring back to original
    coupondf<-subset(coupon_data, coupon_data$Time<=trigger_time)
    cf<-as.numeric(coupondf$CF)
    y_rate<-as.numeric(coupondf$Disc.rate)
    t<-coupondf$Time/360
    presentV<-pv(y_rate,cf,t)
    BondValue_partial<-append(BondValue_partial,presentV)
  }
  
  #E[BondPrice|Trigger]
  price1<-mean(BondValue_partial)
  #E[BondPrice|No Trigger]
  price2<-225 #this is the price advertised by World Bank/paid by investors
  
  #P(H) probability of pandemic with covid-19
  prob1<-0.1468
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

df_expand<-data.frame(cbind(Bond_Price,P_C))
names(df_expand)<-c("Price","Trig.Prob")

#saveRDS(df_expand,"df_expand.rds")

round(Bond_Price,4)
round(P_C,4)

mean(Bond_Price) # this is the price investor should pay


#####################################################################

