library(foreach)
library(doParallel)
library(parallel)

numCores <- detectCores()
cl <- makeCluster(numCores)
registerDoParallel(cl)


results <- foreach(i=inputs) %dopar% {
  processInput(i)
}

#Smoking code re-written from mathematica code. Written by nwenzel August 6, 2016.
#1965, 1970, 1975, 1980, 1985, 1990
library(MHadaptive)
library(mcmc)
library(deSolve)
rm(list = ls())
setwd("/Users/Natasha/Dropbox/Smoking Results/Smokecode-R")

beta.start<-c(0.299,0.2769,0.699,0.689,0.689) #betacc betajj betayy betaja=0, betaya

set.years<-c(1965,1970,1975,1980,1985,1990)
#set.years=1965
################################################################################
#START FOR LOOP
##############################################################################
mcmc.smoke<-function(i)
#for(i in 1:length(set.years))
 {
  year<-set.years[i]

pars4year<-dget(file="/Users/Natasha/Dropbox/Smoking Results/Smokecode-R/FUNC_smoke_pars4year.R",keep.source = TRUE)
#load betas which are dependent on different years
all.pars<<-function(beta,year,...)
  {
  betas<-
  c(betacc=beta[1],
    betajj=beta[2],
    betayy=beta[3],
    betaja=beta[4],
    betaya=beta[5])
  
  parameters<-pars4year(year)
  
  return(c(betas,parameters))
  }


##################################################
#Initial population values
#################################################

init.by.year<-dget(file="/Users/Natasha/Dropbox/Smoking Results/Smokecode-R/INPUT_smoke_pop_pull.R", keep.source=TRUE)
initial.pop<<-init.by.year(year)
year<-1965
###################################################
##Likelihood formula by year
##################################################

#load likelihood formular
llformula <- dget(file="FUNC_smoke_ll_year.R", keep.source=TRUE) 
llikelihood<-llformula()

####################################################
#Hessian matrix by year
####################################################
#to calculate
#hess.by.year<-dget(file="FUNC_init_hessian.R",keep.source=TRUE)
#new.hess<-hess.by.year(set.years, par=beta.start,init=initial.pop)

start.hess<-matrix(c(
-1.573330e-09,  3.475703e-10,  8.925389e-11,  9.794644e-10, -1.200455e-09,
3.475703e-10,  3.417890e-10,  5.981638e-10, -4.784591e-10,  1.026024e-10,
8.925389e-11,  5.981638e-10, -4.654548e-10,  7.994820e-10,  2.296867e-11,
9.794644e-10, -4.784591e-10,  7.994820e-10, -2.839998e-10,  1.553476e-10,
-1.200455e-09,  1.026024e-10,  2.296867e-11,  1.553476e-10, -1.362636e-09), nrow=length(beta.start), byrow=TRUE)
####################################################
#MCMC
#####################################################
#library(fluEvidenceSynthesis)
#require(adaptMCMC)

steps<-rep(1000,100);
sname<-c('MCMCresult1965', 'MCMCresult1970', 'MCMCresult1975', 'MCMCresult1980', 'MCMCresult1985', 'MCMCresult1990')

iterave.save<-function(u)
#for( u in 1:length(steps))
 {
  if(u==1)
      {
master.mcmc<-Metro_Hastings(li_func=llikelihood,
                            pars=beta.start,
                            prop_sigma=start.hess,
                            par_names = c('betacc','betajj','betayy','betaja','betaya'),
                            iterations = steps[u],
                            burn_in = 1000,
                            adapt_par = c(100, 20, 0.5, 0.75),
                            quiet = FALSE,
                            init=initial.pop, year=set.years[i])

save(master.mcmc,file=paste0(sname[i]))
    } else {
      master.mcmc<-Metro_Hastings(li_func=llikelihood,
                                  pars=tail(master.mcmc$trace, n=1),
                                  prop_sigma=master.mcmc$prop_sigma,
                                  par_names = c('betacc','betajj','betayy','betaja','betaya'),
                                  iterations = steps[u],
                                  burn_in = 0,
                                  adapt_par = c(100, 20, 0.5, 0.75),
                                  quiet = FALSE,
                                  init=initial.pop, year=set.years[i])
      master.mcmc<-rbind(master.mcmc,mcmc.result)
      save(master.mcmc,file=paste0(sname[i]))
            }
      } #iterative save loop

sapply(1:length(steps),iterave.save)

}

sapply(programme, mcmc.smoke)
#new.hess<-mcmc.result$prop_sigma
#new.hess==mcmc.result$prop_sigma
#plotMH(mcmc.result)
mcmc.result<-MCMCmetrop1R(fun=llikelihood, theta.init=c(0.3, 0.278, 0.7, 0.7, 0.7), mcmc = 2, burnin = 0, 
                         tune=1, logfun=TRUE, init=initial.pop)

#Hessian matrix
#[,1]          [,2]          [,3]          [,4]          [,5]
#[1,] -1.573330e-09  3.475703e-10  8.925389e-11  9.794644e-10 -1.200455e-09
#[2,]  3.475703e-10  3.417890e-10  5.981638e-10 -4.784591e-10  1.026024e-10
#[3,]  8.925389e-11  5.981638e-10 -4.654548e-10  7.994820e-10  2.296867e-11
#[4,]  9.794644e-10 -4.784591e-10  7.994820e-10 -2.839998e-10  1.553476e-10
#[5,] -1.200455e-09  1.026024e-10  2.296867e-11  1.553476e-10 -1.362636e-09
