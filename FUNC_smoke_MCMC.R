
library(deSolve)
rm(list = ls())
setwd("/Users/Natasha/Dropbox/Smoking Results/Smokecode-R")
#1965, 1970, 1975, 1980, 1985, 1990
period<-c(1,2,3,4,5)

#betacc=0.3, betajj=0.278, betayy=0.7, betaja=0.7, betaya=0.7
unknown.pars<-c(0.278, 0.7, 0.7, 0.7)

all.pars<<-function(beta.start)
  {
  c(betacc=0.1,
    betajj=beta.start[2],
    betayy=beta.start[3],
    betaja=beta.start[4],
    betaya=beta.start[5],
    lc=0.125, 
    lj=0.143, 
    ly=0.25, 
    sigma=0.5, 
    q=(0.074), 
    uj=0.000883, 
    uy=0.00125, 
    ua=0.00801, 
    uc=0.00042,
    phica=runif(1,0.1,1), 
    phicj=runif(1,0.1,0.5))
  }


######################################################################################
#Import DATA sources
####################################################################################

#from https://www.census.gov/popest/data/national/asrh/pre-1980/PE-11.html
pop.data<-read.csv(file="/Users/Natasha/Dropbox/Smoking Results/pop1965.csv",header=FALSE,sep=",",stringsAsFactors=FALSE)
pop.data<-na.omit(as.data.frame(as.numeric(gsub(",","",pop.data[-(1:7),2]))))


initial.pop <<- c(B=sum(pop.data[1:7,]), #ages 0-6
                 Sc=sum(pop.data[8:15,])*0.5, #age 7-14 
                 Ec=sum(pop.data[8:15,])*0.5,
                 Sj=sum(pop.data[16:23,])*0.5,
                 Ej=sum(pop.data[16:23,])*0.5,
                 Sy=sum(pop.data[24:27,])*0.5,
                 Ey=sum(pop.data[24:27,])*0.5,
                 Non=sum(pop.data[28:dim(pop.data)[1],])*0.6,
                 Ha=sum(pop.data[28:dim(pop.data)[1],])*0.1,
                 Na=sum(pop.data[28:dim(pop.data)[1],])*0.2,
                 Qa=sum(pop.data[28:dim(pop.data)[1],])*0.1)

###ODE plots
par(mfrow=c(1,1))
matplot(ode.output$time,ode.output[ ,2:length(ode.output)], type='l')
legend('topleft',legend=2:length(ode.output),col=1:length(ode.output), pch=2)


###################################################
##Likelihood calculation
##################################################

#load likelihood formular
llformula <- dget(file="FUNC_smoke_ll.R", keep.source=TRUE) 
llikelihood<-llformula()


####################################################
#MCMC
#####################################################
library(fluEvidenceSynthesis)
require(adaptMCMC)
burnin<-100; #potatoes
out<-200; #meat
saveiteration<-1; #keeper

ptm <- proc.time() 

mcmc.result<-Metro_Hastings(li_func=llikelihood, pars=c(0.3, 0.278, 0.7, 0.7, 0.7), prop_sigma = new.hess, 
               par_names = NULL, iterations = 50, burn_in = 0, 
               adapt_par = c(100, 20, 0.5, 0.75), quiet = FALSE, init=initial.pop)
proc.time() - ptm

new.hess<-mcmc.result$prop_sigma

new.hess==mcmc.result$prop_sigma

plotMH(mcmc.result)



ptm <- proc.time() 

mcmc.result<-MCMC(p=llikelihood, n=length(test.pars), prop_sigma = NULL, 
                            par_names = NULL, iterations = 2, burn_in = 0, 
                            adapt_par = c(100, 20, 0.5, 0.75), quiet = FALSE, init=initial.pop)
proc.time() - ptm


ptm <- proc.time() 

mcmc.result<-MCMCmetrop1R(fun=llformula, theta.init=c(0.3, 0.278, 0.7, 0.7, 0.7), mcmc = 2, burnin = 0, 
                         tune=1, logfun=TRUE, init=initial.pop)
proc.time() - ptm

#Hessian matrix
[,1]          [,2]          [,3]          [,4]          [,5]
[1,] -1.573330e-09  3.475703e-10  8.925389e-11  9.794644e-10 -1.200455e-09
[2,]  3.475703e-10  3.417890e-10  5.981638e-10 -4.784591e-10  1.026024e-10
[3,]  8.925389e-11  5.981638e-10 -4.654548e-10  7.994820e-10  2.296867e-11
[4,]  9.794644e-10 -4.784591e-10  7.994820e-10 -2.839998e-10  1.553476e-10
[5,] -1.200455e-09  1.026024e-10  2.296867e-11  1.553476e-10 -1.362636e-09
library(maxLik)

[.1]  736531047   125875176 -932214275 -698028398 -734683086
[2,]  125875176 -1176604457 -654329026  489766702 -154681688
[3,] -932214275  -654329026  462803975 -408797623  733188773
[4,] -698028398   489766702 -408797623 -543259962  583001573
[5,] -734683086  -154681688  733188773  583001573 1448289034
maxLik
