#Source file for fluEvidenceSynthesis for changing likelihood. Options for varying age group sizes, risk groups, (and susceptibility, and ascertainment through the starting parameters). Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.

#DO NOT EDIT ANYMORE

#pars<-initial.parameters;
#polymod<-set

#riskratios<-risk.ratios.null
#pars<-unknown.pars
#solutions<-ode.solution
llikelihood.wrap<-function(pars)
{
llikelihood.s <- function(pars,init,...)
  {
    ode.time <- seq(0, 3000, by = 1)
    
ode.eqs<-dget(file="/Users/Natasha/Dropbox/Smoking Results/Smokecode-R/FUNC_smoke_ode_solver.R", keep.source=TRUE)

combo<<-all.pars(pars)

###############################################################
    ## REMEBER PARMS IS NOT THE SAME AS PARS. PARS=UNKNOWN PARAMETERS, PARMS IS ALL PARAMETERS AS GIVEN BY    COMBO
###############################################################
   #init=initial.pop
   #combo<<-all.pars(unknown.pars)
   
    ode.out<-as.data.frame(rk(y=init, times= ode.time, func=ode.eqs, parms=combo, rtol=1e-4,atol=1e-4,maxsteps = 1))
    
    solutions<-ode.out[2000,];
    
    with(as.list(c(combo,solutions)), {
    pop<-Sc+Ec+Sj+Ej+Sy+Ey+Non+Ha+Na+B;
    
    prevalenceC<- Ec/(Ec + Sc)
    prevalenceJ<- Ej/(Ej + Sj)
    prevalenceY<- Ey/(Ey + Sy)
    prevalenceA<- (Ha + Na)/(Ha + Na + Non + Qa) 
    truequit<- Qa/(Na + Ha + Qa) 
    fic <-(betacc*Ec+betacc*phicj*Ej+phica*(Ey+Ha+Na))/pop 
    fij <-(betajj*Ej+betaja*(Ey+Ha+Na))/pop 
    fiy <-(betayy*Ey+betaya*(Na+Ha))/pop
    
            
            
        if(
          (any(combo[1:length(combo)] < 0) || any(combo[1:length(combo)] > 1) || prevalenceA <= 0 
            || prevalenceA >= 1 || prevalenceY <= 0 || prevalenceY >= 1 || prevalenceJ <= 0 
            || prevalenceJ >= 1 || prevalenceC <=0 || prevalenceC >= 1 || truequit <= 0 
            || truequit >= 1 || fic <= 0 || fic>=1 || fij <= 0 || fij>=1 || fiy<=0 || fiy>=1)
          ==TRUE)
          
        {LL<<-(-Inf)}
       
     else  
            {LL<<-dbeta(prevalenceC,150, 8309, log=TRUE) + 
                  dbeta(prevalenceJ,1625, 5549,log=TRUE) + 
                  dbeta(prevalenceY,1281, 2428,log=TRUE) + 
                  dbeta(prevalenceA,7260, 17108,log=TRUE) + 
                  dbeta(truequit,93, 422,log=TRUE) + 
                  dbeta(fic,205, 20995,log=TRUE) + 
                  dbeta(fij,80, 7268,log=TRUE) + 
                  dbeta(fiy,612, 14638,log=TRUE);
            }
   
            })
   return(LL)
  }
  return(llikelihood.s)
}


#####MLE for starting parameters

