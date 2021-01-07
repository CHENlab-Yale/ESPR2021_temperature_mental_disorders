################################################################################
#### Effects of extreme ambient temperature on mental disorders            #####
#### Using daily maximum temperature and ER visits for mental disorders    #####
#### in Erie county during 2009-2015 as an example                         #####
#### Enki Yoo* & Kai Chen**, January 2021                                  #####
#### *Department of Geography, University at Buffalo, SUNY                 #####
#### **Yale School of Public Health                                        #####
################################################################################
##load the packages
library(dlnm) ##Main package for DLNM
library(splines); library(mgcv); library(lubridate)
library(ggplot2); library(knitr);library(kableExtra)


### 1.preparation
load('RData/ED_Weather_PM25.RData')
data      <- ed.df
data$year <- year(data$date)
data$dow  <- as.factor(weekdays(data$date))
data$date.num <- as.numeric(data$date)

###2.DLNM
##model parameter setting, note that these parameters need change in sensitivity analyses
df.temp <- 5   ##note using alternative df.temp for temperature
dfseas  <- 7   ##note using alternativedf for long-term time trend and seasonality
maxlag  <- 14  ##note using alternative maxlag days for sensitivity analyses

data$temp <- data$TMAX

##formula
formula <- 'n ~ cb+dow+ns(date,df=dfseas*length(unique(year)))+ns(PRCP,df=3)' 

#define the crossbasis for the temperature term
argvar <- list(fun="ns", df=df.temp)
lagnk <- 3  ### 3 internal knots equally spaced in the log scale of the lag days
arglag=list(knots=logknots(maxlag,lagnk)) 
cb <- crossbasis(data$temp,lag=maxlag,argvar=argvar,
                 arglag=arglag)
summary(cb)

#run the model and obtain predictions
model <- glm(formula,data,family=quasipoisson,na.action="na.exclude")

#reduction to overall cumulative
red <- crossreduce(cb,model)
coef <- coef(red)
vcov <- vcov(red)

##minimum mortality temperature in this city
predvar <- quantile(data$temp,1:99/100,na.rm=T)
argvar <- list(x=predvar, fun="ns", df=df.temp, 
               Bound=range(data$temp,na.rm=T))
bvar <- do.call(onebasis,argvar)

##get the minimum emergency room visit temperature (MERT)
MMT <- predvar[which.min((bvar%*%coef))]

###get the 2.5% and 97.5% percentile of temperature of the city
tempexcold <- quantile(data$temp,2.5/100,na.rm=T)
tempexhot  <- quantile(data$temp,97.5/100,na.rm=T)
##prediction
pred <- crosspred(bvar,coef=coef,vcov=vcov, model.link="log",by=0.1,cen=MMT,
                  at=c(min(data$temp):max(data$temp),tempexcold,tempexhot))

#########################################################################################
####Figure 3. Maximum daily temperature for mental disorder exposure-response curve ####
#########################################################################################
#plot the exposure-response curve
plot(pred,type="n", ylim=c(0.5,2),yaxt="n",lab=c(6,5,7),
     xlab= expression(paste("Temperature (",degree,"C)")),ylab="RR", 
     main="")
ind1 <- pred$predvar<=MMT
ind2 <- pred$predvar>=MMT
lines(pred$predvar[ind1],pred$allRRfit[ind1],col=4,lwd=1.5)
lines(pred$predvar[ind2],pred$allRRfit[ind2],col=2,lwd=1.5)
axis(2,at=1:10*0.5)

##add temperature distribution (histogram)
breaks <- c(min(data$temp,na.rm=T)-1,seq(pred$predvar[1],
                                             pred$predvar[length(pred$predvar)],length=30),max(data$temp,na.rm=T)+1)
abline(v=MMT,lty=3)
abline(v=c(quantile(data$temp,c(0.025,0.975),na.rm=T)),lty=2)

#######################################################
##Get the RR for heat (97.5th percentile) vs MMT, and 
##cold (2.5th percentile) vs MMT
#######################################################
#Heat
Heat.RR <- pred$allRRfit[as.character(tempexhot)]
Heat.RRlow <- pred$allRRlow[as.character(tempexhot)]
Heat.RRhigh <- pred$allRRhigh[as.character(tempexhot)]
#Cold
Cold.RR <- pred$allRRfit[as.character(tempexcold)]
Cold.RRlow <- pred$allRRlow[as.character(tempexcold)]
Cold.RRhigh <- pred$allRRhigh[as.character(tempexcold)]

##Show the results
#heat
paste0(format(round(Heat.RR,digits = 2), nsmall=2), " (", 
       format(round(Heat.RRlow,digits = 2), nsmall=2),  ", ",
       format(round(Heat.RRhigh,digits = 2), nsmall=2),  ")")
#cold
paste0(format(round(Cold.RR,digits = 2), nsmall=2),  " (", 
       format(round(Cold.RRlow,digits = 2), nsmall=2), ", ",
       format(round(Cold.RRhigh,digits = 2), nsmall=2), ")")