# Numerical Example 1
# Hull White Model
# Figure 1, Figure 2, Figure 3, fixed yield rate, table 3
# Hull white model paramters calibrated.



# This is not my R code except the last part.
#All credits should go to Sang-Heon Lee 
#=========================================================================#
# Financial Econometrics & Derivatives, ML/DL using R, Python, Tensorflow 
# by Sang-Heon Lee 
# https://www.r-bloggers.com/2021/06/hull-white-1-factor-model-using-r-code/
# https://kiandlee.blogspot.com 
#————————————————————————-#
# Numerical Simulation for Hull-White 1 factor model
#=========================================================================#

library(Rfast)  # colCumProds
library(rstudioapi)

graphics.off()  # clear all graphs
rm(list = ls()) # remove all files from your workspace

# change directory to where the script located

my_d <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(my_d)


# Functions for numerical Integration

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  I(t) = Int_0^t sigma(s)^2 A exp(Bs) ds
#————————————————————-#
#           t
#  I(t) = ∫  σ(u)^2 A exp(Bu) du  
#          0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fI<-function(t, A, B, lt.HW) {
  M <- 0
  value <- 0
  
  tVol <- lt.HW$tsig   # volatility tenor
  Vol  <- lt.HW$sigma  # volatility vector
  nVol <- lt.HW$nsig   # # of volatility
  
  # find Maximum M from j which is t_j < t
  M <- ifelse(length(which(tVol<=t))==0,1,max(which(tVol<=t))+1)
  
  # summation part
  if (B==0) {
    if (M==1) value <- value + Vol[1]^2*A*t
    else {
      for (i in 1:(M-1)) {
        add <- Vol[i]^2*A*(tVol[i] - ifelse(i==1,0,tVol[i-1]))
        value <- value + add
      }
      add <- Vol[ifelse(M==(nVol+1),M-1,M)]^2*A*(t-tVol[M-1])
      value <- value + add
    }
  }
  else {
    if (M==1) { value <- value + Vol[1]^2*A/B*(exp(B*t)-1)}
    else {
      for (i in 1:(M-1)) {
        add <- Vol[i]^2*A/B*
          (exp(B*tVol[i])-ifelse(i==1,1,exp(B*tVol[i-1])))
        value <- value + add
      }
      add <- Vol[ifelse(M==(nVol+1),M-1,M)]^2*A/B*
        (exp(B*t)-exp(B*tVol[M-1]))
      value <- value + add
    }
  }
  return(value)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  A(s,t)=e^(-Int_s^t a(v) dv)
#————————————————————-#
#                   s
#  A(s,t) = exp( -∫ a(v)dv )  
#                  t
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fA<-function(s, t, lt.HW) {
  tau <- lt.HW$tkap       # tau
  K1  <- lt.HW$kappa[1]   # short-term kappa
  K2  <- lt.HW$kappa[2]   # long-term kappa
  
  if      (tau <= s) f <- exp(-K2*(t-s))
  else if (t < tau ) f <- exp(-K1*(t-s))
  else               f <- exp(-K1*(tau-s)-K2*(t-tau))
  
  return(f)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  B(s,t)=Int_s^t e^(-Int_t^u a(v) dv) du
#————————————————————-#
#             t       u
#  B(s,t) = ∫ exp( -∫ a(v)dv ) du  
#            s       t
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fB1<-function(s, t, kappa) {return((1 - exp(-kappa*(t-s)))/ kappa)}

fB<-function(s, t, lt.HW) {
  tau <- lt.HW$tkap       # tau
  K1  <- lt.HW$kappa[1]   # short-term kappa
  K2  <- lt.HW$kappa[2]   # long-term kappa
  
  if      (tau <= s) f <- fB1(s, t, K2)
  else if (t < tau ) f <- fB1(s, t, K1)
  else  f <- fB1(s,tau,K1)+exp(-K1*(tau-s))*fB1(tau,t,K2)
  
  return(f)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Zeta(t) = Int_0^t σ(u)^2 e^(-2 Int_u^t a(v) dv) du
#————————————————————-#
#             t               t
#  Zeta(t) = ∫ σ(u)^2 exp( -2∫ a(v)dv ) du  
#            0               u
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fZeta<-function(t, lt.HW) {
  tau <- lt.HW$tkap       # tau
  K1  <- lt.HW$kappa[1]   # short-term kappa
  K2  <- lt.HW$kappa[2]   # long-term kappa
  
  if (t < tau) f = exp(-2*K1*t)*fI(t,1,2*K1,lt.HW) 
  else  f = exp(-2*K2*(t-tau)-2*K1*tau)*fI(tau,1,2*K1,lt.HW)+ 
    exp(-2*K2*t)*(fI(t,1,2*K2,lt.HW)-fI(tau,1,2*K2,lt.HW))
  
  return(f)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Z(t) = Int_0^t σ(u)^2 e^(-Int_u^t a(v) dv) B(u,t) du 
#————————————————————-#
#          t              t
#  Z(t) = ∫ σ(u)^2 exp( -∫ a(v)dv ) B(u,t) du  
#         0              u
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fZ1<-function(t, kappa, lt.HW) {
  I1 = exp(  -kappa*t)*fI(t,1,  kappa, lt.HW) / kappa
  I2 = exp(-2*kappa*t)*fI(t,1,2*kappa, lt.HW) / kappa
  return(I1 - I2)
}

fZ<-function(t, lt.HW) {
  tau <- lt.HW$tkap       # tau
  K1  <- lt.HW$kappa[1]   # short-term kappa
  K2  <- lt.HW$kappa[2]   # long-term kappa
  
  if (t < tau) 
    f = fZ1(t, K1, lt.HW)
  else {
    I1 = exp(-K2*(t-tau))*fZ1(tau, K1, lt.HW)
    I2 = exp(-K2*(t-tau))*fB(tau,t,lt.HW)*
      exp(-2*K1*tau)*fI(tau,1,2*K1,lt.HW)
    I3 = exp(-K2*t) * fI(tau, 1, K2, lt.HW) / K2
    I4 = exp(-2*K2*t) * fI(tau, 1, 2*K2, lt.HW) / K2
    f =  I1 + I2 + fZ1(t, K2, lt.HW) - I3 + I4
  }
  return(f)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#  Omega(t,T) = Int_0^t sigma(s)^2 [B(s,t)^2 - B(s,T)^2] ds
#————————————————————-#
#                t                     
#  Omega(t,T) = ∫ σ(s)^2 [B(s,t)^2 - B(s,T)^2] ds
#               0                    
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
fOmega<-function(t, T, lt.HW) {
  return(-fB(t,T,lt.HW) * (2.0*fZ(t,lt.HW) + 
                             fB(t,T,lt.HW)*fZeta(t,lt.HW)))
}

#=========================================================================#
#             Main : Hull-White 1 Factor Model Simulation
#=========================================================================#

#—————————————————————-#
# Information List for the Hull-White model
#—————————————————————-#
# - tkap  : threshold year which divide mean-reversion speed
# - kappa : mean-reversion speed parameters
# - tsig  : maturity vector for volatility parameters
# - sigma : volatility parameter vector
# - tDF   : maturity vector for spot rates
# - rc    :spot rates curve
#—————————————————————-#

# list object which contain Hull-White model related information
#rates are set as of July 07, 2017 to match market LIBOR and ICE SWAP rates
#https://www.theice.com/iba/ice-swap-rate
#https://www.theice.com/iba/libor

lt.HW <- list(
  tkap  = 10, 
  #this value obtained from Kladívko, Kamil, and Tomáš Rusý. "Maximum likelihood estimation of the Hull–White model." Journal of Empirical Finance 70 (2023): 227-247.
  #see Table 5 alpha values for 2014, 2015, 2016. So set it at 0.0001 since setting 0 is not possible. Also no differentiation between two regime
  kappa = c(0.0001, 0.0001), 
  tsig  = c(1.0, 2.0, 3.0, 4.0, 5.0, 6.0,7.0, 8.0, 9.0, 10.0,15),
  #these values calculated using historical ICE Swap rates,. Not he best way. But nothing to no since no access to swaption implied volatltities
  sigma = c(0.003447103,0.002813146,	0.003613892,	0.00443339,	0.004694461,	0.004693241,	0.004644499,	0.004538977,	0.004480984,	0.004420839,	0.004263847),
  tDF   = c(1/365,1/52,1/12,2/12,3/12,6/12,1.0, 2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, 10.0,15.0),
  #These are spot rates from USD LIBOR+ ICE SWAP rates. The spot rates taken on July -7, 2017
  rc    = c(0.0117389,0.0119,0.0122633,0.01255,0.0130522,0.0146544,0.01471,0.01653,0.018,0.0192,0.02022,0.02109,0.02184,0.02247,0.02302,0.02351,0.02508)
)

# Add other information to list 
lt.HW$nDF  <- length(lt.HW$tDF)   # # of spot
lt.HW$nsig <- length(lt.HW$sigma) # # of vol 
lt.HW$nkap <- length(lt.HW$kappa) # # of kappa 

# Check for Numerical Integration Functions for HW1F
m.temp <- matrix(NA,15,5)
colnames(m.temp) <- c("I", "B", "Zeta", "Z", "Omega")
for(i in 1:15) {
  m.temp[i,1] <- fI    (i, 2, 3, lt.HW)
  m.temp[i,2] <- fB    (0.5,  i, lt.HW)
  m.temp[i,3] <- fZeta (i,       lt.HW)
  m.temp[i,4] <- fZ    (i,       lt.HW)
  m.temp[i,5] <- fOmega(0.5,  i, lt.HW)
}
print("Check for Numerical Integration Functions for HW1F")
print(m.temp)

# Discount Factor 
lt.HW$DF   <- exp(-lt.HW$tDF*lt.HW$rc)

#—————————————————————-#
# Preprocessing for simulation
#—————————————————————-#

# Simulation information
denom.1y     <- 252    # # of dt in 1-year # we use 252 since that is average number of trading dates per year

# t : valuation date, T : maturity
lt.HW.sim    <- list(t=0, T=4, dt=1/denom.1y, nscenario =5000)

lt.HW.sim$nt <- round(lt.HW.sim$t*denom.1y,0)
lt.HW.sim$nT <- round(lt.HW.sim$T*denom.1y,0)

# spit the time axis by dt
v.Ti <- seq(lt.HW.sim$dt, lt.HW.sim$T, length = lt.HW.sim$nT)  

#—————————————————————-#
# Linear Interpolation of spot rate curve
#—————————————————————-#
# rule=2 : For outside the interval [min(x), max(x)], 
#          the value at the closest data extremeis used.
#—————————————————————-#
frci <-approxfun(x=lt.HW$tDF, y=lt.HW$rc, rule=2)

v.rci <- frci(v.Ti)             # interpolated spot rates
v.DFi <- exp(-v.Ti*v.rci) # interpolated DF

#—————————————————————-#
# temporary use for blog width adjustment
#—————————————————————-#
sim <- lt.HW.sim
par <- lt.HW
dt  <- lt.HW.sim$dt

# standard normal random error
set.seed(123456)

# predetermined vector
v.A <- v.Zeta <- v.dZeta.sqrt <- v.B <- v.Omega <- rep(0, sim$nT)

for (n in 1:sim$nT) {
  v.A[n]     <- fA    (v.Ti[n]-dt, v.Ti[n], par)
  v.Zeta[n]  <- fZeta (v.Ti[n],             par)
  v.B[n]     <- fB    (v.Ti[n]-dt, v.Ti[n], par)
  v.Omega[n] <- fOmega(v.Ti[n]-dt, v.Ti[n], par)
}

v.dZeta.sqrt <- c(sqrt(v.Zeta[1]),
                  sqrt(v.Zeta[-1]-v.A[-1]^2*v.Zeta[-sim$nT])) 

# selecting some indices because plotting is time-consuming 
v.idx.sample <- sample(1:sim$nscenario, 500)

#—————————————————————-#
# Simulation Part
#—————————————————————-#

# interpolated discount factor from initial yield curve
v.P0 <- v.DFi 
# ratio of bond price P(0,t+dt)/P(0,t)
v.P0T_P0T1 <- c(v.P0[1]/1,v.P0[-1]/v.P0[-sim$nT])

m.P.ts   <- matrix(0, sim$nT, sim$nscenario ) # P(t,t+dt)
m.Rsc.ts <- matrix(0, sim$nT, sim$nscenario ) # short rate

# Simulate from now on.

# for n=1
m.P.ts  [1,] <- v.P0T_P0T1[1]
m.Rsc.ts[1,] <- -log(m.P.ts[1,])/dt
xt <- rnorm(sim$nscenario, 0, 1)*v.dZeta.sqrt[1]

for(n in 2:sim$nT) {
  print(n)
  m.P.ts[n,] <- v.P0T_P0T1[n]*exp(-xt*v.B[n]+0.5*v.Omega[n])  
  xt <- xt*v.A[n] + rnorm(sim$nscenario, 0, 1)*v.dZeta.sqrt[n]
}

m.Rsc.ts <- -log(m.P.ts)/dt      # short rates(current spot rates)
m.DF.ts  <- colCumProds(m.P.ts)  # Discount Factors
m.R0T.ts <- -log(m.DF.ts)/v.Ti   # future spot rates

## plot paths
t <- seq(dt, lt.HW.sim$T, dt)

x11(width=6, height=5);
matplot(m.P.ts[,v.idx.sample], type="l", lty=1,
        xlab="Mat",ylab="P(t,t+dt)",main="Simulated ZCB")
x11(width=6, height=5);
matplot(m.Rsc.ts[,v.idx.sample], type="l", lty=1,
        xlab="Days",ylab="R(t,t+dt)",main="Simulated Short Rate")
x11(width=6, height=5);
matplot(m.DF.ts[,v.idx.sample], type="l", lty=1,
        xlab="Mat",ylab="DF(0,T)"  ,main="Simulated Discount Factor")
x11(width=6, height=5);
matplot(m.R0T.ts[,v.idx.sample], type="l", lty=1,
        xlab="Mat",ylab="R(0,T)"   ,main="Simulated Spot Rate")

######################################################################################################
#######################################################################################################
###My coding######################################################################################
##################################################################################################
# We want 6M USD LIBOR based on short rates
# For each short rate path, starting from a given day, we calculate average of short rate predicted for next 6-month. This is 6M USD LIBOR rate for that date.If average rate is negative, we still use it.
# Since we have simulated 5000 paths, we get 5000 possible 6M USD LIBOR rates for that given date
# Finally we take average of those rates.
#we assume year is 252 days (trading days) and six month consists of 126 days
library(zoo)
dim(m.Rsc.ts)

# Define a new function that computes the rolling mean and handles NAs
matrix_short_rate<-m.Rsc.ts
#matrix_short_rate<-matrix(c(1:16), nrow = 4, byrow = TRUE) # for test purposes

rollmean_na <- function(x, k) {
  if (all(is.na(x))) {
    return(rep(NA, length(x)))
  } else {
    return(rollmean(x, k = k, fill = NA, align = "left"))
  }
}

# Apply the function to each column
matrix_avg_rate <- apply(matrix_short_rate, 2, rollmean_na, k = 126)
#you want to omit NA values at the end due to roll mean, as well as convert it to just a numeric vector for further processing.
PRED_6M_LIBOR<-as.numeric(na.omit(apply(matrix_avg_rate, 1, mean, na.omit=T)))


# Load the data (6-month LIBOR rates)
# Note: In this example, we assume that the data is available in a CSV file
# Please adjust this to fit your actual data source
data <- read.csv("libor_rates_6M_ICE.csv")
names(data)<-c("Date","LIBOR.6M")

# The date format in the CSV file should be "YYYY-MM-DD"
data$Date <- as.Date(data$Date, format = "%m/%d/%Y")
df<-data[(data$Date<=as.Date("07/15/2020",format = "%m/%d/%Y")),]
df<-na.omit(df)

df2<-data[(data$Date<=as.Date("07/15/2020",format = "%m/%d/%Y")) & (data$Date>=as.Date("07/08/2017",format = "%m/%d/%Y")),]
df2<-na.omit(df2)
#Check how many predictions you need
len<-length(df2$LIBOR.6M) #You need 764 predictions
df2$PRED.LIBOR.6M<-(PRED_6M_LIBOR[1:len])*100 # convert rates to percentages


Coupon_dates<-as.Date(c("08/15/2017","09/15/2017","10/16/2017","11/15/2017","12/15/2017",
                        "01/15/2018","02/15/2018","03/15/2018","04/16/2018","05/15/2018","06/15/2018","07/16/2018","08/15/2018","09/17/2018","10/15/2018","11/15/2018","12/17/2018",
                        "01/15/2019","02/15/2019","03/15/2019","04/15/2019","05/15/2019","06/17/2019","07/15/2019","08/15/2019","09/16/2019","10/15/2019","11/15/2019","12/16/2019",
                        "01/15/2020","02/17/2020","03/16/2020","04/15/2020","05/15/2020","06/15/2020","07/15/2020"),format = "%m/%d/%Y")
Coupon_rates<-df2[df2$Date%in%Coupon_dates,]
Coupon_rates

#find the number of dsys between successive coupon payment date. We add bond issue date at this point
dataset1<-append(as.Date("07/07/2017",format = "%m/%d/%Y"),Coupon_dates[-length(Coupon_dates)])
dataset2<-Coupon_dates
diff_time<-abs(difftime(dataset1, dataset2, units = "days"))

#Face value of the bond
FaceValue<-225 #in millions
spread<-6.5 # how much extra over LIBOR Rate

rates<-Coupon_rates$PRED.LIBOR.6M

# Remember on the issuance date July 07, 2017 we knew true value of libor rate. So we need to add  that. Also we do no need rates on July 15, 2017 since 
# that day coupon based on June 15, 2015 6M LIBOR rates

r1<-df[df$Date%in%c(as.Date("07/07/2017",format="%m/%d/%Y")),"LIBOR.6M"]
rates2<-append(r1,rates[-length(rates)])
rates3<-(rates2+spread)/100

# Calculate the predicted coupon amount based on FaceValue *rate3*(actual days/360)

Coupon_Amount<-as.numeric(FaceValue*rates3*(diff_time/360))
Coupon_Amount #in millions

###################
#calculate real coupon amount based on real LIBOR rates
rates_re<-Coupon_rates$LIBOR.6M

# Remember on the issuance date July 07, 2017 we knew tru value of libor rate. So we need to add  that. Also we do no need rates on July 15, 2017 since 
# that day coupon based on June 15, 2015 6M LIBOR rates

r1_re<-df[df$Date%in%c(as.Date("07/07/2017",format="%m/%d/%Y")),"LIBOR.6M"]
rates2_re<-append(r1_re,rates_re[-length(rates_re)])
rates3_re<-(rates2_re+spread)/100

# Calculate the real coupon amount based on FaceValue *rate3*(actual days/360)

Coupon_Amount_re<-as.numeric(FaceValue*rates3_re*(diff_time/360))
Coupon_Amount_re #in millions

library(dplyr)
df_coupon<-as.data.frame(bind_cols(Coupon_dates,Coupon_rates$LIBOR.6M,Coupon_Amount_re,Coupon_rates$PRED.LIBOR.6M,Coupon_Amount))
names(df_coupon)<-c("Date","Real.6M.LIBOR","Real.Coupon","Pred.6M.LIBOR","Pred.Coupon")
#add a row on July 07, 2017 data
df_coupon<-rbind(data.frame(Date=as.Date("07/07/2017",format="%m/%d/%Y"),Real.6M.LIBOR=1.4,Real.Coupon=NA,Pred.6M.LIBOR=NA, Pred.Coupon=NA),df_coupon)
#Remove LIBOR rates for July 15, since is is not needed
df_coupon[length(df_coupon$Date),2]<-NA
df_coupon[length(df_coupon$Date),4]<-NA
###########################

#Now calculate the yield rate investor will receive if no pandemic occur. The investor calculating this based on the price he/she paid for the bond
# World bank issue price is 100%. That mean investors paid full face price. usually one investor may not have leverage to affect on price the bond the agency willing to sell.
# The price based on bids and number of buyers. 100% issuance mean there may have been more buyers who are willing to buy, so no choice but buy at face value
# Once the person buy the bond, the real uncertainty is coupons and whether it default. Coupons are uncertain because floating rate.

#Based on predicted coupon amounts lets calculate yield rate if no pandemic occur.
# Now we are discounting to July 07, 2017

diff_time2<-abs(difftime(as.Date("07/07/2017",format = "%m/%d/%Y"), dataset2, units = "days"))
year_fraction_360<-diff_time2/360



#find yield to maturity/internal rate of return
cash_flows<-Coupon_Amount
#Add Amount paid to the begining as negative value
cash_flows<-append(-225,cash_flows)
#add zero to year fraction to denote when price was paid
year_fraction_360_2<-append(0,year_fraction_360)

# Add 225M (redemption amount) to last coupon, 
cash_flows[length(cash_flows)]<-cash_flows[length(cash_flows)]+225

#find irr https://stackoverflow.com/questions/11660187/any-r-package-available-to-calculate-irr-from-uneven-payments-on-specific-dates 

cf<-cash_flows
npv <- function(i, cf=cf, t=year_fraction_360_2){ sum(cf/(1+i)^t)} 
irr <- function(cf) { uniroot(npv, c(0,1), cf=cf)$root } 
irr(cf)
round(irr(cf)*100, 4)

#lets create a NPV profile

rate_i<-seq(0, 0.3, by=0.01)
NPV_set<-sapply(rate_i,npv,cf=cf, t=year_fraction_360_2 )

NPV_data<-as.data.frame(cbind(rate_i,NPV_set))
names(NPV_data)<-c("rate","NPV")

#############plot NPV #########################################
p2<- ggplot()+
  geom_line(data=NPV_data, aes(x =rate , y = NPV),color = 'black', size = 0.75)+
  #geom_hline(yintercept=0,  color='black', size=0.75)
  labs(title = 'NPV of World Bank pandemic bond',
       x = 'rate',
       y = 'NPV(USD millions)',
       caption = '') + 
  theme_bw(base_family = "TT Times New Roman") 

x11();p2;ggsave(filename ='npv.png', plot = p2)



########################################################################
#to price under pandemic scenarios we create following data frame

df3<-cbind(diff_time2,Coupon_Amount)
df3<-as.data.frame(df3)
names(df3)<-c("Time","CF")
#add redemption amount to last coupon
df3$CF[length(df3$CF)]<-df3$CF[length(df3$CF)]+225

saveRDS(df3, file = "coupon.rds")

##################Plot#################################################
library(ggplot2)
p1 <- ggplot()+
geom_line(data=df2, aes(x =Date , y = LIBOR.6M),color = 'black', size = 0.75)+
  geom_line(data=df2, aes(x =df2$Date , y = df2$PRED.LIBOR.6M),color = 'blue', size = 0.75)+
  labs(title = 'Observed vs Predicted 6Month USD LIBOR Rates',
       x = '',
       y = 'Interest Rate',
       caption = '') + 
  theme_bw(base_family = "TT Times New Roman") 


x11(); p1;ggsave(filename ='obs_pred_libor.png', plot = p1)


 