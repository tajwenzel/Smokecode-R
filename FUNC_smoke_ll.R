#Source file for fluEvidenceSynthesis for changing likelihood. Options for varying age group sizes, risk groups, (and susceptibility, and ascertainment through the starting parameters). Adapted from code by Baguelin 2013 and Edwin van Leevvan, by TajWenzel.

#pars<-initial.parameters;
#polymod<-set

#riskratios<-risk.ratios.null
#pars<-combo.pars
#solutions<-ode.solution

smoke.llikelihood<-function()
{
  llikelihood.s <- function(pars,solutions,...)
  {
    
    with(as.list(c(pars,solutions)), {
    pop<-Sc+Ec+Sj+Ej+Sy+Ey+Non+Ha+Na+B;
    
    prevalenceC<<- Ec/(Ec + Sc)
    prevalenceJ<<- Ej/(Ej + Sj)
    prevalenceY<<- Ey/(Ey + Sy)
    prevalenceA<<- (Ha + Na)/(Ha + Na + Non + Qa) 
    truequit<<- Qa/(Na + Ha + Qa) 
    fic <<-(betacc*Ec+betacc*phicj*Ej+phica*(Ey+Ha+Na))/pop 
    fij <<-(betajj*Ej+betaja*(Ey+Ha+Na))/pop 
    fiy <<-(betayy*Ey+betaya*(Na+Ha))/pop
    
            
            
        if(
          (any(pars[1:length(pars)] < 0) || any(pars[1:length(pars)] > 1) || prevalenceA <= 0 
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

