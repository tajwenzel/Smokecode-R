
library(deSolve)
rm(list = ls())
setwd("/Users/Natasha/Dropbox/Smoking Results")
#1965, 1970, 1975, 1980, 1985, 1990
period<-c(1,2,3,4,5)

unknown.pars<-c(betacc=0.3, betajj=0.278, betayy=0.7, betaja=0.7, betaya=0.7)

known.pars<-function()
  {
  c(lc=0.125, lj=0.143, ly=0.25, sigma=0.5, q=(0.074), phica=runif(1,0.1,1), 
  phicj=runif(1,0.1,0.5), uj=0.000883, uy=0.00125, ua=0.00801, uc=0.000427)
  };

combo.pars<-c(unknown.pars,known.pars())

lc=0.125
lj=0.143
ly=0.25
sigma=0.5
q=(0.074)
phica=runif(1,0.25,1) 
phicj=runif(1,0.25,1)
uj=0.000883 
uy=0.00125 
ua=0.00801 
uc=0.000427
betacc=0.1
betajj=0.278
betayy=0.7
betaja=0.7
betaya=0.7
######################################################################################
#Differential Equations
####################################################################################

ode.time <- seq(0, 3000, by = 1)
#ode.time<-c(0,5000)


#from https://www.census.gov/popest/data/national/asrh/pre-1980/PE-11.html
pop.data<-read.csv(file="pop1965.csv",header=FALSE,sep=",",stringsAsFactors=FALSE)
pop.data<-na.omit(as.data.frame(as.numeric(gsub(",","",pop.data[-(1:7),2]))))


initial.pop <- c(B=sum(pop.data[1:7,]), #ages 0-6
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

#init=initial.pop
ode.storage<-as.data.frame(matrix(NaN,ncol=length(initial.pop),nrow=length(ode.time)))
### The function ###
solver.function <- function(t,init=initial.pop,pars=parms,...) {
    B<-init[1];
    Sc<-init[2];
    Ec<-init[3];
    Sj<-init[4];
    Ej<-init[5];
    Sy<-init[6];
    Ey<-init[7];
    Non<-init[8];
    Ha<-init[9];
    Na<-init[10];
    Qa<-init[11];
    
    Pop<-B+Sc+Ec+Sj+Ej+Sy+Ey+Non+Ha+Na+Qa;
    
    with(as.list(pars), {
      
      dB<-((Sc+Ec)*uc+(Sj+Ej)*uj+(Sy+Ey)*uy+(Non+Ha+Na+Qa)*ua+B*uc)-B*uc-B*0.14
  
      dSc <- B*0.14+(-Sc*betacc*Ec)/Pop-(Sc*betacc*phicj*Ej)/Pop-(Sc*betacc*phica*(Ey+Ha+Na))/Pop-Sc*lc-Sc*uc;
      
      dEc <- (Sc*betacc*Ec)/Pop+(Sc*betacc*phicj*Ej)/Pop+(Sc*betacc*phica*(Ey+Ha+Na))/Pop-Ec*lc-Ec*uc;
      
      dSj <-(-Sj*betajj*Ej)/Pop-(Sj*betaja*(Ey+Ha+Na))/Pop-Sj*lj+Sc*lc-Sj*uj;
      
      dEj <- (Sj*betajj*Ej)/Pop+(Sj*betaja*(Ey+Ha+Na))/Pop-Ej*lj+Ec*lc-Ej*uj;
      
      dSy <- (-Sy*betayy*Ey)/Pop-(Sy*betaya*(Na+Ha))/Pop-Sy*ly-Sy*uy+Sj*lj;
      
      dEy <- (Sy*betayy*Ey)/Pop+(Sy*betaya*(Na+Ha))/Pop-Ey*ly+Ej*lj-Ey*uy;
      
      dNon<- Sy*ly-Non*ua;
      dHa <- Ey*ly*0.1-Ha*ua;
      dNa <- Ey*ly*0.9-q*Na-Na*ua;
      dQa <- q*Na-Qa*ua;
      
      #####ODE list
      
      stored <-list(c(dB,dSc,dEc,dSj,dEj,dSy,dEy,dNon,dHa,dNa,dQa))
      
    })
  }


ode.output<-as.data.frame(rk(y=initial.pop, times= ode.time, func=solver.function, parms=combo.pars, rtol=1e-4,atol=1e-4,maxsteps = 1))

ode.solution<-ode.output[2000,]

par(mfrow=c(1,1))
matplot(ode.output$time,ode.output[ ,2:length(ode.output)], type='l')
legend('topleft',legend=2:length(ode.output),col=1:length(ode.output), pch=2)

###################################################
##Likelihood calculation
##################################################
llformula <- dget(file="FUNC_smoke_ll.R", keep.source=TRUE) 
llikelihood<-llformula()

try<-llikelihood.s(combo.pars,ode.solution)
try2<-llikelihood(combo.pars,ode.solution)
####################################################
#MCMC
#####################################################
library(fluEvidenceSynthesis)

burnin=100
nbatch=100
blen=1


ptm <- proc.time() 
mcmc.result <- adaptive.mcmc(lprior = 0, llikelihood, 
                             nburn=burnin,
                             initial = initial.pop,
                             nbatch = out, blen = saveiteration,
                             outfun= NULL, acceptfun= NULL,
                             pars=initial.pop, solution=ode.solution )
proc.time() - ptm


